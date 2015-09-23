
#defines the library paths within which the libraries are found
.libPaths(c("C:/Users/rds/Documents/R/win-library/3.0", "C:/Program Files/R/R-3.0.2/library"))
#.libPaths("C:/Program Files/R/R-3.0.2/library")

library(plyr)
library(vegan)
library(stringr)
library(gdata)
library(scales)
# Set the filepath

setwd("/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/cancer_grade/freq_filtered/NEW/working_folder/")
  #setwd("C:/Users/rsutherland/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/cancer_grade/logistic/tumor_stage_cov")

topDirectory<-("/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/cancer_grade/freq_filtered/NEW/working_folder/")


models<-c("age_gender","age_gender_vars","age_gender_proteins","age_gender_vars_proteins")
#models<-c("age_gender","age_gender_vars","age_gender_proteins","age_gender_proteins_vars")



clin_params<- c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE)
clin_v_params<- c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE)
clin_p_params<-c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE)
clin_v_p_params<-c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE)


#dataframe of the parameters for each model type
params<-as.data.frame(cbind(clin=clin_params,clin_var= clin_v_params, clin_pro=clin_p_params, clin_var_pro_params=clin_v_p_params))


# choose to balance the outcome classes through downsampling of the majority class
balance<-TRUE

#use FalseDR correction for variable filtering or simply a rank threshold
FalseDR=TRUE


# are the SNV TV/TS frequencies included. For calculating FDR of protein variables
SNV_Types=TRUE
SNV_Types_Only=TRUE

## use clinicalData only?
clinical_only<-FALSE

#
covs<-c("age","stage","gender","tumourType")#covariates to include change this when necessary
#age,stage,gender, tumourType
##########################################################################################
##########################################################################################
# the filter to use

## true =frequency of mutation false= fisher exact test for enrichment in ns and functional mutations

freqfilter<- TRUE
##########################################################################################
# use proteins in the model?

proteins=TRUE
proteinsOnly=FALSE



# set the random seed

set.seed(10)

#
#
#
# The data file in vcf-like format.
#scores.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/sep_2013/sex_specific_removed/input"
scores.path <- "/Users/rds/Dropbox/PhD/tumour_classifier_data/sep_2013/input"
#scores.path <- "C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/sep_2013/input"

#scores.files<- "pancan12_cleaned.maf"


#Tissue source site file
TSS.path<-"/Users/rds/Dropbox/PhD/tumour_classifier_data/sep_2013/TSS"
#TSS.path<-"C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/sep_2013/TSS"
TSS.file<- "tissueSourceSite.csv"

######################################################################################################################
### metadata
#########

clinical.path<-"/Users/rds/Dropbox/PhD/tumour_classifier_data/SynapsePanCancer/PanCan12/tumour_grade/clinical"# load all clinical data at once
#clinical.path<-"C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/SynapsePanCancer/PanCan12/tumour_grade/clinical"# load all clinical data at once

##################################################################################################################################################
## load scores
scores.files<-unlist(list.files(scores.path))
dataFiles<-file.path(scores.path,scores.files)
scores<-createScoreTable(dataFiles)
center<-seqCentre(scores$Tumor_Sample_Barcode)
scores["center"]<- center


#scores<-labelCancer(scores,scores.files)
#scores<- rbind.data.frame(scores[[1]],scores[[2]])

##remove x y and mt chromosomes

scores<- scores[-which(scores$Chromosome=="X"|scores$Chromosome=="Y"|scores$Chromosome=="MT"),]


#silent mutations
scores.silent<-scores[which(scores$Variant_Classification=="Silent"),]
scores.function<-scores[which(scores$Variant_Classification!="Silent"),]


scores.functionSamples<- unique(scores.function$sampleID)# I need to remove the samples containing 0 silent mutations
scores.silentSamples<- unique(scores.silent$sampleID)# The samples containing at east 1 silent mutation
scores<- extractIntersectSamples(scores$sampleID,scores.silentSamples, scores)


#This line is here in case there are samples that have sequence information, but no phenotype data
#scores<-extractIntersectSamples(scores$sampleID,intersectSamples,scores)# this produces the same result as in the original code, but much faster
scores<- scores[order(scores$sampleID),]# final ordered scores dataframe...

##############################################################################################################################################
## load metadata and intersection between metadata and phenotype
##############################################################################################################################################

clinicalData.directories<-unlist(list.dirs(clinical.path))[-1]# load all clinical data at once

clinicalData.cancers<-lapply(seq(1:length(clinicalData.directories)), function(x) strsplit(clinicalData.directories[x],"/",)[[1]][11])



createClinicalData<- function(clinicalData.directories, clinicalData.cancers){
  
  clinicals<-lapply(seq(1:length(clinicalData.cancers)),function(x) clinicalDataTable(clinicalData.directories[x],strsplit(clinicalData.directories[x],"/",)[[1]][11]))
 names(clinicals)<- clinicalData.cancers
    return(clinicals)
}

allClinicalData<-createClinicalData(clinicalData.directories, clinicalData.cancers)# loads all of the clinical data

#change the tumour_grade variable name to neoplasm_histological_grade for UCEC to match the other two cancers
colnames(allClinicalData$UCEC)[which(colnames(allClinicalData$UCEC)=="tumor_grade")]<- "neoplasm_histologic_grade"
colnames(allClinicalData$UCEC)[which(colnames(allClinicalData$UCEC)=="gynecologic_tumor_grouping_figo_stage")]<- "tumor_stage"
# now the grade will be included in the reduced variable set



# list of all clinical data variables

allClinicalNames<-lapply(seq(1:length(allClinicalData )), function(x) colnames(allClinicalData[[x]]))

getReducedNameSet<-Reduce(intersect, allClinicalNames)# the clinical data names common to all clinical data files


#uses the list of clinical datatables as input together with the reduced nameset common to all clinical datatables
# and outputs the list of clinical datatables using only the reduced set of variable names.
getReducedClinicalData<- function(data,nameset){
  relevantData<-lapply(seq(1:length(data)), function(x) data[[x]][,nameset])
  names(relevantData)<-names(data)
  return(relevantData)
}

allClinicalDataReduced<- getReducedClinicalData(allClinicalData, getReducedNameSet)
allClinicalDataReduced<-do.call("rbind", allClinicalDataReduced)# all cancer clinical data in a single table! IMPORTANT


clinicalData<-allClinicalDataReduced
#

###################################################################################################################
# tumor grade phenotype work

# tumor grade frequencies

tumor_grade_freq<- table(clinicalData$neoplasm_histologic_grade, useNA="ifany")
write.table(tumor_grade_freq, file="tumour_grade_freq.txt", sep="\t", append=FALSE)
barplot(tumor_grade_freq, las=3, xlab=as.character(names(tumor_grade_freq)))
barplot(tumor_grade_freq, las=3)
# remove tumours without grade information
clinicalData<-clinicalData[-which(is.na(clinicalData$neoplasm_histologic_grade)),]
clinicalData<-clinicalData[-which(clinicalData$neoplasm_histologic_grade=="GB"|clinicalData$neoplasm_histologic_grade=="GX"),]# removing the tumours that could not be graded

###########################################################################################################
possible_grades<- unique(clinicalData$neoplasm_histologic_grade)# the gradings in my dataset
###defining high and low grade
# tumour grade high or low
grade<- rep(0, dim(clinicalData)[1])

high_low<-c("low","high","high","low","high", "low","high", "low")# grading assignments in the order of possible_grades shown below

#possible_grades   [1] "G2" "G4" "G3" "G1" "Grade 3" "Grade 2" "High Grade" "Grade 1"   

##assign high and low grade to the remaining samples
for (i in 1:(length(high_low))){
grade[which(clinicalData$neoplasm_histologic_grade==possible_grades[i])]<-high_low[i] 
}
#add grade to the clinical Data table. this is now ready for analysis
clinicalData<-cbind(clinicalData, grade)

############################################################################################################
## define stages
################################################################
# remove tumours without stage information
clinicalData<-clinicalData[-which(is.na(clinicalData$tumor_stage)),]


possible_stages<- unique(clinicalData$tumor_stage)

# group the tumour stages in to 1,2,3,4
stage_t<- str_replace(clinicalData$tumor_stage, "[A*,B*,C*]", "")
stage_u<-gsub("Stage ","",stage_t)
stage_v<- str_replace(stage_u, "[1,2]", "")

#stage<-as.factor(stage_v)# The stage variable to combine with the other clinical data. Stages have been standardised to I,II,III,IV

#recode stage to a numeric
stage_v[which(stage_v=="I")]<-1
stage_v[which(stage_v=="II")]<-2
stage_v[which(stage_v=="III")]<-3
stage_v[which(stage_v=="IV")]<-4
# recode to numeric
stage<- as.numeric(stage_v)# stage is now coded as a numeric ordinal variable For the purposes of logistic regression this is useful
                           # because stage is ordinal and you can only progress up the scale.

# dichotomise stage to be high/low the same way I did for tumour grade
stage[stage<3]<- "low_stage"
stage[stage!="low_stage"]<- "high_stage"
stage1<- factor(stage, levels=c("low_stage", "high_stage"))

clinicalData<- cbind(clinicalData, stage)# clinical data with the stage variable added
###################################################################################################################
#random sample of clinical data names

#sampleSubset<- sample(clinicalData$SampleIds, floor(length(clinicalData$SampleIds)/20), replace=FALSE)
#sampleSubsetIndex<-sapply(seq(1:length(sampleSubset)), function(x) which(clinicalData$SampleIds==sampleSubset[x]))
#clinicalData<- clinicalData[sampleSubsetIndex,]# random subset of Clinical Data
###################################################################################################################



##################################################################################################
#extract samples with phenotype and clinical data from clinical data

samplesAnalysis<- intersect(clinicalData$SampleIds, unique(scores$sampleID))
clinicalData2<- extractIntersectSamples(clinicalData$SampleIds,samplesAnalysis, clinicalData)
clinicalData2<-clinicalData2[order(clinicalData2$SampleIds),]

###################################################################################################
#

### classifier for each cancer type

dim(clinicalData2)

OVclinicalData2<-clinicalData2[which(clinicalData2$tumourType=="OV"),]
KIRCclinicalData2<-clinicalData2[which(clinicalData2$tumourType=="KIRC"),]
UCECclinicalData2<-clinicalData2[which(clinicalData2$tumourType=="UCEC"),]


#clinicalData2<- UCECclinicalData2

#important variables to include
#clinicalData2$stage##not available for BLCA
#clinicalData2$age_at_initial_pathologic_diagnosis
#clinicalData2$gender


if(length(which(is.na(clinicalData2$age_at_initial_pathologic_diagnosis)))>0){
  clinicalData2<- clinicalData2[-which(is.na(clinicalData2$age_at_initial_pathologic_diagnosis)),]
if(length(which(is.na(clinicalData2$gender)))>0){
  clinicalData2<- clinicalData2[-which(is.na(clinicalData2$gender)),]
  }
  else if(length(which(is.na(clinicalData2$gender)))>0){
      clinicalData2<- clinicalData2[-which(is.na(clinicalData2$gender)),]
}
}
#if(length(which(is.na(clinicalData2$days_to_last_followup)))>0){
#  clinicalData2<- clinicalData2[-which(is.na(clinicalData2$days_to_last_followup)),] 
#}

#remove samples that have clinical data missing for included clinical predictors
#clinicalData2<- clinicalData2[-which(is.na(clinicalData2$age_at_initial_pathologic_diagnosis)),]
# remove samples with missing gender data
#clinicalData2<- clinicalData2[-which(is.na(clinicalData2$gender)),]
samplesAnalysis<-intersect(clinicalData2$SampleIds, samplesAnalysis)


#extract samples with phenotype and clinical data from scores
scores<- extractIntersectSamples(scores$sampleID,samplesAnalysis, scores)

##################################################################################################


mutations<-scores[which(scores$Variant_Classification!="Silent"),]
colnames(mutations)<- colnames(scores)



mutationTable<- table(mutations$Hugo_Symbol, mutations$sampleID)
ranksGenesTest<- rownames(mutationTable)



mutationMatrix<-matrix(mutationTable,length(mutationTable[,1]),length(mutationTable[1,]))# This is my mutation matrix to be used to calculate sample distances

rownames(mutationMatrix)<- rownames(mutationTable)
colnames(mutationMatrix)<- colnames(mutationTable)
mutationMatrixBinary<- mutationMatrix

mutationMatrixBinary[which(mutationMatrix>1)]<-1# reassigns genes carrying more than one mutation to 1.
# all genes carrying more than one mutation in each individual ahve been set to one. This matrix
# now indicates if an individual carries at least one mutation per gene.

# logical mutation matrix all 1s converted to TRUE and all 0s converted to FALSE.
mutationMatrixLogical <- matrix(as.logical(mutationMatrixBinary),length(ranksGenesTest),dim(mutationMatrixBinary)[2])
rownames(mutationMatrixLogical)<- rownames(mutationTable)
colnames(mutationMatrixLogical)<- colnames(mutationTable)



scores.SNV<- scores[which(scores$Variant_Type=="SNP"),]



#swapStrand<- function(negativeVariant){
#
#  v<- negativeVariant
#  
#if(v=="A"){
#  v<-"T"
#}else if(v=="T"){
#  v<-"A"
#}else if(v=="C"){
#  v<- "G"
#}else if(v=="G"){
#  v<-"C"
#}
#return(v)
#}


#swapNegativeStrand<- function(scoreTable){
#
#  # negative strand index
#  indexNegative<-which(scoreTable$strand=="-1")
#  scoreTable.negative<-scoreTable[indexNegative,]
#  
#  
#  scoreTable.Negative.variant<-scoreTable.negative$variant
#
#swapped<-sapply(scoreTable.Negative.variant, swapStrand, USE.NAMES=FALSE)
#scoreTable$variant[indexNegative]<- swapped
#
#return(scoreTable)
#
#}
#scores.SNV<- swapNegativeStrand(scores.SNV)



#scores.SNV by sample
scores.SNV.samples<-lapply(samplesAnalysis, function(x) scores.SNV[which(scores.SNV$sampleID==x),])



####Function to get SNV type frequency across the entire dataset
getSNVType<- function(scoresDataFrame, b1,b2,b3,b4){
  substitution1<- length(which(scoresDataFrame$reference== b1 & scoresDataFrame$variant== b2))
  substitution2<- length(which(scoresDataFrame$reference== b3 & scoresDataFrame$variant== b4))
  total<- sum(substitution1,substitution2)
  return(total)
}

## scores$variant is the same as Tumour_seq_Allele2
## Both Match_Norm_Seq_Allele1 and Match_Norm_Seq_Allele2 are identical and identical to scores$reference
## i may need to change the strand of begative mutations to be the positive strand

SNV_Type_freqs<-function(scoresDataFrame){
  
  AC_TG<- getSNVType(scoresDataFrame, "A","C","T","G")
  AG_TC<- getSNVType(scoresDataFrame, "A","G","T","C")
  AT_TA<- getSNVType(scoresDataFrame, "A","T","T","A")
  CA_GT<- getSNVType(scoresDataFrame, "C","A","G","T")
  CG_GC<- getSNVType(scoresDataFrame, "C","G","G","C")
  CT_GA<- getSNVType(scoresDataFrame, "C","T","G","A")
  
  SNV_Types_freq<- c(AC_TG,AG_TC, AT_TA, CA_GT, CG_GC, CT_GA)
  names(SNV_Types_freq)<- c("A>C/T>G","A>G/T>C", "A>T/T>A", "C>A/G>T", "C>G/G>C", "C>T/G>A")
  
  return(SNV_Types_freq)
}


#substitution type frequency across dataset
# I can place this in within my loop to get SNV_type frequncies for each of my classes.
SNV_Types_freq_dataset<-SNV_Type_freqs(scores.SNV)

##SNV substitutions by sample

SNV_Types_freq_samples<-sapply(scores.SNV.samples, function(x) SNV_Type_freqs(x))#number of SNV_Type substitutions for each sample.Use this as a predictor in my random forest model
colnames(SNV_Types_freq_samples)<- sapply(1:length(scores.SNV.samples), function(x) unique(scores.SNV.samples[[x]]$sampleID))

SNV_Types_freq_samples<- t(SNV_Types_freq_samples)


####mutation types##### BEGIN here to create plots!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mutationTypes<-t(table(scores$Variant_Classification, scores$sampleID))
mutationTypes<- mutationTypes[,-9]# do not use the silent frequencies
#aggregate(mutationTypes, by=list(clinicalData2$tumourType), FUN=mean)

barplot(rowMeans(mutationTypes))

###########
# VariantFreqTable

variantFreqTable<- cbind(SNV_Types_freq_samples, mutationTypes)

###########################################################################################################

#######################################################################################################################################
#plots

#plots of the samples with clinical and mutation data
# tumor grade frequencies

tumor_grade_freq_clinical_and_mutation<- table(clinicalData2$neoplasm_histologic_grade, useNA="ifany")
write.table(tumor_grade_freq_clinical_and_mutation, file="tumour_grade_freq_samples_with_clin_and_mut_data.txt", sep="\t", append=FALSE)
barplot(tumor_grade_freq, las=3, main="grade frequencies for samples with both mutation and clinical data",xlab=as.character(names(tumor_grade_freq)))

######barplot of hig/low grade frequency
library(colorspace)
colors2<- rainbow_hcl(length(table(outcome)),c=50,l=70, start=0)
png("tumour_freq.png", width=1200, height=800, bg="transparent")
barplot(table(outcome), main="high and low grade tumour frequency", col=colors2, border=NA, ylim=c(0,800), cex.axis=1.5, cex.names=1.5)
dev.off()
#######barplot of tumor stage frequency

barplot(table(clinicalData2$stage), main ="tumor stage frequency")


####### grade by stage chi-square
# contingency table
as.matrix(table(clinicalData2$grade, clinicalData2$stage))

#chisquare test
chisq.test(as.matrix(table(clinicalData2$grade, clinicalData2$stage)))
#######################################################################################################################################







samples<- colnames(mutationMatrixLogical)
sample_TSS<- sapply(samples, function(x) substr(x, start=5, stop=6))
# in the TSS file the numeric categories belo 10 are single digit so I must convert my categories generated from the TCGA ID from "02" etc to just "2".
sample_TSSpos1<-substr(sample_TSS,start=1, stop=1)
sample_TSS[which(sample_TSSpos1==0)] <-substr(sample_TSS[which(sample_TSSpos1==0)], start=2, stop=2)

TSS<- read.table(file.path(TSS.path,TSS.file) , header=TRUE, sep = ",", quote = "", stringsAsFactors = FALSE)
cancer_Type1<- sapply(sample_TSS, function(x) TSS[which(TSS$TSS.Code==x),3])



groups<- unlist(unique(cancer_Type1))
#groupsClass<- unlist(sapply(cancer_Type1, function(x) which(groups==x)-1))
#groupsClass<- as.factor(groupsClass)# can be my outcome variable




##############################################################################
outcome<- clinicalData2$grade # my outcome is tumour grade in this analysis
outcome<- reorder(outcome,new.order=c("low","high"))
names(outcome)<- clinicalData2$SampleIds
write.table(file="grade_type_contingency", table(clinicalData2$neoplasm_histologic_grade,clinicalData2$grade), row.names=TRUE, col.names=TRUE, sep="\t")
#I can use this to extract the appropriate variables from the MutationMatrixLogical before I calculate F.tests or SNV_Type frequencies.

#gender table
genderTable<- table(outcome,clinicalData2$gender)
Ftestgender<-fisher.test(genderTable, alternative="two.sided")
sink("gender_stats.txt", append=TRUE)
print(list("gender Table",genderTable,"gender Fisher test",Ftestgender))
sink()

# tumourType by grade Chi-square
chisq.test(clinicalData$grade, clinicalData$tumourType)
table(clinicalData$grade, clinicalData$tumourType)

#age table
ageDataFrame <- as.data.frame(cbind(clinicalData2$age_at_initial_pathologic_diagnosis,outcome))
ageAOV<- aov(clinicalData2$age_at_initial_pathologic_diagnosis~outcome)
ageAOVANOVA<-anova(ageAOV)
pairwise.t.test(clinicalData2$age_at_initial_pathologic_diagnosis, outcome, p.adjust="none")
sink("age_stats", append=TRUE)

print(list(ageDataFrame,ageAOV, ageAOVANOVA))
sink()
# there are no statistically significant differences between the tumour types in terms of age

# fisher test to find genes associated with grade
fisherTables<-lapply(seq(1:dim(mutationMatrixLogical)[1]), function(x) table(outcome, mutationMatrixLogical[x,]))
names(fisherTables)<-rownames(mutationMatrixLogical)
fisherTestGrade<-lapply(seq(1:dim(mutationMatrixLogical)[1]), function(x) fisher.test(fisherTables[[x]], alternative="two.sided"))
names(fisherTestGrade)<-rownames(mutationMatrixLogical)
fisherTestPval<-sapply(seq(1:dim(mutationMatrixLogical)[1]), function(x) fisherTestGrade[[x]]$p.value)
significantTests<-which(fisherTestPval<0.05/dim(mutationMatrixLogical)[1])
significantPval<-fisherTestPval[significantTests]
#for output of significant differences
significantGenes<- rownames(mutationMatrixLogical)[significantTests]# the overall genes distinguishing between the tumour grades
significantFisherTestResults<- lapply(significantGenes,function(x) fisherTestGrade[x])
names(significantFisherTestResults)<-significantGenes
significantTables<-lapply(significantGenes, function(x) fisherTables[[x]])
names(significantTables)<-significantGenes

sigFisherResultsOutput<-data.frame(cbind(significantGenes,significantPval))
write.table(sigFisherResultsOutput, "significant_overall_fisher_test_grade_results.txt", sep="\t")

sink("significant_overall_fisher_test_grade_tables.txt", append=TRUE)
print(significantTables)
sink()


############

#distance1<-dist(mutationMatrixBinary, method="manhattan", diag=FALSE, upper=FALSE)
#distance2<- vegdist(mutationMatrixBinary, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE)
#mds1<-cmdscale(distance2, k=2, eig=TRUE, add=FALSE, x.ret=FALSE)# object too large for memory
#pcAnalysis<-prcomp(t(mutationMatrixBinary), scale=TRUE, center=TRUE)
#pcAnalysis2<-prcomp(mutationMatrixBinary, scale=TRUE, center=TRUE)
#plot(pcAnalysis$x[,1], pcAnalysis$x[,2])
#plot(pcAnalysis2$rotation[,1], pcAnalysis2$rotation[,2])


#put all the code in here
#twoSubTypes<- unlist(strsplit(comps[i],split="_"))
#subType1<- twoSubTypes[1]
#subType2<- twoSubTypes[2]
#dir.create(paste("C:/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/cancer_grade/", subType1,"_", subType2, sep=""))
#setwd(paste("C:/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/cancer_grade/", subType1,"_", subType2, sep=""))
##########################################################################################################################
#Random Forest starts Here
##########################################################################################################################

##Defining a training and test set for 
library(caret)
library(mlbench)
library(Hmisc)
library(randomForest)
library(e1071)

set.seed(10)
##############################################################################
##The two cancerClasses
class1<-outcome[which(outcome=="low")]
class2<-outcome[which(outcome=="high")]



class1_2<-c(class1,class2)
y<- outcome
ySamples<-names(y)
minorityClass<- names(which(table(y)==min(table(y))))
majorityClass<- names(which(table(y)==max(table(y))))



balanceClasses<-function(y){
minorityClass<- y[which(y==names(which(table(y)==min(table(y)))))]
majorityClass<- y[which(y==names(which(table(y)==max(table(y)))))]

mClass<-sample(majorityClass,length(minorityClass), replace=FALSE)

sampled<-c(names(mClass),names(minorityClass))# the names of the samples selected
sampledy<-y[intersect(names(y),sampled)]# the balanced outcome variable
return(sampledy)
}

if(balance==TRUE){
y<- balanceClasses(y)#balanced classes

  
ySamples<-names(y)}

set.seed(10)
TenFolds<- createFolds(y,k=3,list=FALSE, returnTrain=FALSE)

###############################################################################
#are the folds stratified?
table(y)/length(y)

table(TenFolds)

table(y[which(TenFolds==1)])/length(which(TenFolds==1))
###############################################################################
#test and train samples
testSamples<- ySamples[which(TenFolds==1)]
trainSamples<- ySamples[which(TenFolds!=1)]

#test and train outcomes
testOutcome<- y[which(TenFolds==1)]
trainOutcome<- y[which(TenFolds!=1)]

##############################################################################################################################

####################################################################################################################################
#####################################################################################################################################
######each analysis starts here


for(i in seq(1,length(models))){
  i<-1
 #model<- "age_gender_tumour"
  model<-models[i]
  balance<-params[1,i]
  FDR<-params[2,i]
  SNV_Types<-params[3,i]
  SNV_Types_Only<-params[4,i]
  clinical_only<-params[5,i]
  freqfilter<-params[6,i]
  proteins<-params[7,i]
  proteinsOnly<-params[8,i]

print(balance)
print(freqfilter)


dir.create(paste(topDirectory, model, sep=""))
setwd(paste(topDirectory, model, sep=""))




#test and train data matrices
testMutationMatrixLogical<-as.data.frame.matrix(t(mutationMatrixLogical[,testSamples]))
trainMutationMatrixLogical<-as.data.frame.matrix(t(mutationMatrixLogical[,trainSamples]))



#################################################################################################################################################
## filter out proteins according to the number of samples in which they carry mutations ################################

tumourType<-unique(clinicalData2$tumourType)

train_ClinicalData2<- clinicalData2[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds),]
#the names of the proteins carried through to the next stage, by cancer type

if(freqfilter){
filterMuts<- function(mutationMatrix){
  #filter out the genes that do not carry mutations in more than 5% of the samples
  #5% of samples is 970/20
  
  FivePercent<- ceiling(dim(mutationMatrix)[1]/20)
  mutationMatrix_proteins<-names(which(colSums(mutationMatrix)>FivePercent))
  
  return(mutationMatrix_proteins)
}
}else{
 # funcsTotal<- sum(rowSums(variantFreqTable_genes)[-c(8,9)])
#  nonfuncsTotal<- sum(rowSums(variantFreqTable_genes)[c(8,9)])
  
  
  ftestfilter<-function(variantFreqGene, funcsTotal, nonfuncsTotal){    
    funcs<-sum(variantFreqGene[-c(8,9)])
    nonfuncs<- sum(variantFreqGene[c(8,9)])
    
    # functional in gene
    A<- funcs
    #nonfunc in gene
    B<- nonfuncs
    # func not it gene
    C<- funcsTotal-funcs
    # non funcs not in gene
    D<- nonfuncsTotal-nonfuncs
    
    ftable<- matrix(data= c(A,B,C,D), nrow=2, ncol=2,byrow=TRUE)
    
    
    filtergene<-fisher.test(ftable, alternative="greater")$p.value
    return(filtergene)
}


# the train samples split according to cancer type
cancergroups<-split(trainSamples,train_ClinicalData2$tumourType)


# function to get the list of gene names that are enriched for functional mutations for each cancer type.
getgenes<- function(groups,scores,i){
  # extract the relevant mutations for the tumour type in question in the training set
  groupsScoresList<- lapply(seq(1:length(groups[[i]])), function(x) scores[which(scores$sampleID==groups[[i]][x]),])
  groupsScores<- do.call("rbind", groupsScoresList)
  
  
  #For each tumourType find the number of proteins significantly enriched for functional mutation regardless of tumour grade
  
  # mutation type count per gene across all samples coming from the "scores" table
  variantFreqTable_genes<-t(table(groupsScores$Hugo_Symbol, groupsScores$Variant_Classification))
  
  #TP53 fnctional variants i.e. not silent or RNA variants 16880
  
  ## silent and non silent frequencies across the entire sample
  
  funcsTotal<- sum(rowSums(variantFreqTable_genes)[-c(8,9)])
  nonfuncsTotal<- sum(rowSums(variantFreqTable_genes)[c(8,9)])
  
  
  filterPvalues<- sapply(seq(1:length(colnames(variantFreqTable_genes))), function(x) ftestfilter(variantFreqTable_genes[,x], funcsTotal, nonfuncsTotal))
  names(filterPvalues)<- colnames(variantFreqTable_genes)
  
  #useFDR
  retainedGenes<-cbind(which(p.adjust(filterPvalues,method="BH")<=0.1), filterPvalues[which(p.adjust(filterPvalues,method="BH")<=0.1)])

  
  genenames<-rownames(retainedGenes)
  return(genenames)
  return(retainedGenes)
}

}



if(freqfilter){
filtered_protein_Names<-lapply(tumourType, function(x) filterMuts(trainMutationMatrixLogical[train_ClinicalData2$tumourType==x,]))
names(filtered_protein_Names)<- tumourType

# the unique list of proteins carrying mutations in at least 5% of at least cancer type
proteins_analysis<- unique(unlist(filtered_protein_Names))
proteins_analysisN<- length(proteins_analysis)

}else{
  filtered_genes<- sapply(seq(1:length(cancergroups)), function(x) getgenes(cancergroups,scores,x))
  
  #the proteins to be included in the trainMutationMatrix
  proteins_analysis<-unique(unlist(filtered_genes))
  proteins_analysisN<- length(proteins_analysis)
}

write.table(proteins_analysis, file="filtered__proteins_for_log_regression.txt", append =FALSE, sep="\t")


###########################################################################################################################
######## if the above qc filter removes all of the proteins a whole section must be skipped out
##########################################################################################################################
if(proteins_analysisN>0){
  ##redefining the train and testMutationMatrices.
  
  trainMutationMatrixLogical<- trainMutationMatrixLogical[,proteins_analysis]
testMutationMatrixLogical<-testMutationMatrixLogical[,proteins_analysis]

trainSamples<- rownames(trainMutationMatrixLogical)
testSamples<- rownames(testMutationMatrixLogical)
}else{
  # if all of the proteins were removed using the protein filtering step keep 2 so that the code below doesn't have to change
 proteins_analysis<- colnames(trainMutationMatrixLogical)[c(seq(1:200))]
  
  trainMutationMatrixLogical<- trainMutationMatrixLogical[,proteins_analysis]#this is a hack
  testMutationMatrixLogical<-testMutationMatrixLogical[,proteins_analysis]# this is a hack
  
  trainSamples<- rownames(trainMutationMatrixLogical)
  testSamples<- rownames(testMutationMatrixLogical)
}

# may have to redefine the samples used if non of them have any proteins carrying mutations in the new matrix.
####################################################################################################################################################

#######################################################################################################################################


setdiff(class1,intersect(names(trainOutcome),class1))


###############################################################################################################################
##add extra test samples from the majority training class


if(unique(class1)==majorityClass){
  
  extraTestSamples<-setdiff(names(class1),intersect(names(trainOutcome),names(class1)))## the majority samples left out of the training set
  extraTest<-t(mutationMatrixLogical[proteins_analysis,extraTestSamples])
  extraTest<-extraTest[-which(extraTestSamples%in%rownames(testMutationMatrixLogical)),]#remove majority samples that are already present in the test set
}else{
  extraTestSamples<-setdiff(names(class2),intersect(names(trainOutcome),names(class2)))## the majority samples left out of the training set
  extraTest<-t(mutationMatrixLogical[proteins_analysis,extraTestSamples])
  extraTest<-extraTest[-which(extraTestSamples%in%rownames(testMutationMatrixLogical)),]#remove majority samples that are already present in the test set
  
}

##add extraTest samples to the testMutation matrix
testMutationMatrixLogical<-rbind(testMutationMatrixLogical, extraTest)

##add extra samples to the testSamples vector
testSamples<- c(testSamples, rownames(extraTest))
## add extra samples to the testOUtcome vector
testOutcome<- factor(outcome[c(names(testOutcome),rownames(extraTest))])
######################################################################################################################


######################################################################################################################################
######### filter out genes with mutations only in one class high or low because they are perfect predictors and log regression will complain

#trueRemove<-sapply(seq(1,length(colnames(trainMutationMatrixLogical))), function(x) 0%in%table(trainOutcome,trainMutationMatrixLogical[,x])[,TRUE])


#trainMutationMatrixLogical<-trainMutationMatrixLogical[,which(trueRemove==FALSE)]# remove the perfect predictors
#testMutationMatrixLogical<-testMutationMatrixLogical[,which(trueRemove==FALSE)]


######################################################################################################################
##Fisher's exact test variable selecton ######## to be used if you want to do variable selection manually
nonZeroGenes<- which(colSums(trainMutationMatrixLogical)!=0)
trainMutationMatrixLogical<- trainMutationMatrixLogical[,nonZeroGenes]
testMutationMatrixLogical<- testMutationMatrixLogical[,nonZeroGenes]



if("stage"%in% covs){
metadata_vars<-c("age_at_initial_pathologic_diagnosis","stage","gender","tumourType")
}else if("tumourType"%in% covs){
  metadata_vars<-c("age_at_initial_pathologic_diagnosis","gender","tumourType")
}else{
metadata_vars<-c("age_at_initial_pathologic_diagnosis","gender")
}

metadata_vars_index<-unlist(sapply(metadata_vars, function (x) which(colnames(clinicalData2)==x)))

metadata_vars_names<-names(metadata_vars_index)
metadata_vars_names[1]<- strsplit(metadata_vars_names[1], "_")[[1]][1]# if age is included

##filter metadata_var_names based on the covariates I want to inlude in the analysis (defined in "covs" at the top of the script)
#metadata_vars_names<-intersect(covs, metadata_vars_names)

##change this each time I wish to change the covariates in my analysis
covariates<-match(covs,metadata_vars_names)# index of the metadata variables I want to take to analysis


##names_string for files
namesString<-paste(metadata_vars_names[covariates], collapse="_")



##################################################################
# define the test and train metadata variables
trainMutationMetadataVars<- as.data.frame(sapply(seq(1,length(metadata_vars_names)), function(x) clinicalData2[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds),metadata_vars_index[x], drop=FALSE]))
colnames(trainMutationMetadataVars)<- metadata_vars_names
trainMutationMetadataVars<-trainMutationMetadataVars[covariates]

testMutationMetadataVars<- as.data.frame(sapply(seq(1,length(metadata_vars_names)), function(x) clinicalData2[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds),metadata_vars_index[x], drop = FALSE]))
colnames(testMutationMetadataVars)<- metadata_vars_names
testMutationMetadataVars<-testMutationMetadataVars[covariates]
#############################################################################################################################
## add SNV_Type measure
## add SNV_TYPE measures to the train/testmutationMatrixLogical

if(SNV_Types){

  if(SNV_Types_Only){
    testMutationMatrixLogical<-as.data.frame(variantFreqTable[match(testSamples,rownames(variantFreqTable)),], stringsAsFactors=TRUE)
    trainMutationMatrixLogical<-as.data.frame(variantFreqTable[match(trainSamples,rownames(variantFreqTable)),], stringsAsFactors=TRUE)
    extraVars<-dim(variantFreqTable)[2]+1
    
  }
else if (proteins){
testMutationMatrixLogical<-as.data.frame(cbind(variantFreqTable[match(testSamples,rownames(variantFreqTable)),],testMutationMatrixLogical), stringsAsFactors=TRUE)
trainMutationMatrixLogical<-as.data.frame(cbind(variantFreqTable[match(trainSamples,rownames(variantFreqTable)),],trainMutationMatrixLogical), stringsAsFactors=TRUE)
}## The extra variables other than proteins in the mutationMatrixlogical
extraVars<-dim(variantFreqTable)[2]+1
}else if(proteinsOnly){
  #i no protins made it through the filterr then I need to make this just the metadatavariables
  if(proteins_analysisN==0){
    testMutationMatrixLogical<-testMutationMetadataVars
    testMutationMatrixLogical<-as.data.frame(cbind(outcome=testOutcome,testMutationMatrixLogical))
    trainMutationMatrixLogical<- trainMutationMetadataVars
    trainMutationMatrixLogical<-as.data.frame(cbind(outcome=trainOutcome,trainMutationMatrixLogical))
    
    testMutationMatrixLogical2<-testMutationMatrixLogical
    trainMutationMatrixLogical2<-trainMutationMatrixLogical 
    extraVars<-1
  }else{
    testMutationMatrixLogical<- testMutationMatrixLogical
    trainMutationMatrixLogical<-trainMutationMatrixLogical
    extraVars<-1
    }
}else if (clinical_only){#just the metadataVariables
  testMutationMatrixLogical<-testMutationMetadataVars
  testMutationMatrixLogical<-as.data.frame(cbind(outcome=testOutcome,testMutationMatrixLogical))
  trainMutationMatrixLogical<- trainMutationMetadataVars
  trainMutationMatrixLogical<-as.data.frame(cbind(outcome=trainOutcome,trainMutationMatrixLogical))
  
  testMutationMatrixLogical2<-testMutationMatrixLogical
  trainMutationMatrixLogical2<-trainMutationMatrixLogical  
  extraVars<-1}


###############START HERE ######################################3
#####offset the covariates depending on if the tumourType covariate is included
if(is.na(match("tumourType", covs))){
  offset_TTypes<- 2
}else{offset_TTypes<-3}# 






if(proteinsOnly|SNV_Types_Only){
  ##for just proteins
  if(proteinsOnly){
  namesString<- paste(namesString,"_proteins", sep="")
  }else{
    namesString<-paste(namesString,"_variantTypes", sep="")
  }
  extraVars<- length(covariates)+1
  
  testMutationMatrixLogical<-as.data.frame(cbind(outcome=testOutcome,testMutationMetadataVars[covariates],testMutationMatrixLogical), stringsAsFactors=TRUE)
  trainMutationMatrixLogical<-as.data.frame(cbind(outcome=trainOutcome,trainMutationMetadataVars[covariates],trainMutationMatrixLogical), stringsAsFactors=TRUE)
  genes<- colnames(trainMutationMatrixLogical)[-seq(1:(length(covariates)+1))]
  
  
  allGLM<-lapply(seq(length(covariates)+2, (length(genes)+(length(covariates)+1))), function(x) summary(glm(outcome~., data = trainMutationMatrixLogical[,c(seq(1:(length(covariates)+1)),x)], family="binomial")))
  names(allGLM)<- genes

  
  ############################################
  #finding out the variant variables with only on value across the two groups
  oneValue<-sapply(seq(length(covariates)+2, dim(trainMutationMatrixLogical)[2]),function(x) length(colnames(table(trainMutationMatrixLogical$outcome,trainMutationMatrixLogical[,x]))))
  fixedVariables<-which(oneValue==1)
  if(length(fixedVariables)>0){
    #reassign the fixed variable to be a pvlue of 1
    allGLM[[fixedVariables]]$coefficients<-rbind(allGLM[[fixedVariables]]$coefficients,c(NA,NA,NA,1))
  }
  ###########################################
  
  
  allGLMPvals<- sapply(seq(1,length(genes)), function(x) allGLM[[x]]$coefficients[length(covariates)+offset_TTypes,4])
  names(allGLMPvals)<- genes
  
  
  
  if(FalseDR==TRUE){
    
    significant_allGLMPvals<- which(p.adjust(allGLMPvals,method="BH")<=0.1)
    significant_allGLMPvals_genes<- names(significant_allGLMPvals)
  }else{
    #allGLMPvals<-allGLMPvals[order(allGLMPvals, decreasing=FALSE)]# order the genes by logistic regression P-value
    significant_allGLMPvals<- which(allGLMPvals<=0.05)
    significant_allGLMPvals_genes<-names(significant_allGLMPvals)
  }
  
  # genes that predict classification as high/low grade signifcantly well adjusted for age and gender
  significant_GLM_output<-data.frame(genes=significant_allGLMPvals_genes, Pval= allGLMPvals[significant_allGLMPvals_genes])
  # write to file
  write.csv(significant_GLM_output, file=(paste("high","_","low","_2_thirds_train_logistic_significant_",namesString,".csv", sep="")), row.names=FALSE)
  bonferroni_pvalue<-paste(c("bonferroni p-value = 0.05/",length(genes),",",0.05/length(genes), sep=","))
  cat(bonferroni_pvalue, file=(paste("high","_","low","_2_thirds_train_logistic_significant_",namesString,".csv", sep="")), append=TRUE)
  
  
  colIndex_sig_genes<-significant_allGLMPvals# col index of the genes in the mutation matrices
  trainMutationMatrixLogical2<- trainMutationMatrixLogical[,c(seq(1, (length(covariates)+1)),(length(covariates)+1+colIndex_sig_genes))]
  testMutationMatrixLogical2<-testMutationMatrixLogical[,c(seq(1, (length(covariates)+1)),(length(covariates)+1+colIndex_sig_genes))]
  
  
  
  
  
  #proteins and variant types analysis
}else if(proteins&SNV_Types){
namesString<- paste(namesString,"_proteins_variants", sep="")

if(length(covariates==1)){
  testMutationMatrixLogical<-as.data.frame(cbind(testMutationMetadataVars[covariates],testMutationMatrixLogical), stringsAsFactors=TRUE)
  trainMutationMatrixLogical<- as.data.frame(cbind(trainMutationMetadataVars[covariates],trainMutationMatrixLogical), stringsAsFactors=TRUE)  
}else{
  testMutationMatrixLogical<-as.data.frame(cbind(testMutationMetadataVars[,covariates],testMutationMatrixLogical), stringsAsFactors=TRUE)
  trainMutationMatrixLogical<- as.data.frame(cbind(trainMutationMetadataVars[,covariates],trainMutationMatrixLogical), stringsAsFactors=TRUE)  
}
testMutationMatrixLogical<-as.data.frame(cbind(outcome=testOutcome,testMutationMatrixLogical))
trainMutationMatrixLogical<-as.data.frame(cbind(outcome=trainOutcome,trainMutationMatrixLogical))
genes<- colnames(trainMutationMatrixLogical)[-seq(1,length(covariates)+1)]# gets the gene names from the mutation matrix without the clinical data names


#glmTrainMutationMatrixLogical<-trainMutationMatrixLogical[,c(1,2,3,i)]

allGLM<-lapply(seq(length(covariates)+2, length(genes)+(length(covariates)+1)), function(x) summary(glm(outcome~., data = trainMutationMatrixLogical[,c(seq(1:(length(covariates)+1)),x)], family="binomial")))

#Remove the variables that have a fixed value


names(allGLM)<- c(genes)


############################################
#finding out the variant variables with only on value across the two groups
oneValue<-sapply(seq(length(covariates)+2, (length(covariates)+extraVars)),function(x) length(colnames(table(trainMutationMatrixLogical$outcome,trainMutationMatrixLogical[,x]))))
fixedVariables<-which(oneValue==1)
if(length(fixedVariables)>0){
  #reassign the fixed variable to be a pvlue of 1
allGLM[[fixedVariables]]$coefficients<-rbind(allGLM[[fixedVariables]]$coefficients,c(NA,NA,NA,1))
}
###########################################

allGLMPvals<- sapply(seq(1,length(genes)), function(x) allGLM[[x]]$coefficients[length(covariates)+offset_TTypes,4])
names(allGLMPvals)<- c(genes)




if(FalseDR==TRUE){
  
  significant_allGLMPvals<- which(p.adjust(allGLMPvals,method="BH")<=0.1)# the tests where at that p-value, less than 10* of positives are expected to be false positives. 
  significant_allGLMPvals_genes<- names(significant_allGLMPvals)
}else{
  #allGLMPvals<-allGLMPvals[order(allGLMPvals, decreasing=FALSE)]# order the genes by logistic regression P-value
  significant_allGLMPvals<- which(allGLMPvals<=0.05)
  significant_allGLMPvals_genes<-names(allGLMPvals)[significant_allGLMPvals]
}

# genes that predict classification as high/low grade signifcantly well adjusted for age and gender
significant_GLM_output<-data.frame(genes=significant_allGLMPvals_genes, Pval= allGLMPvals[significant_allGLMPvals_genes])
# write to file
write.csv(significant_GLM_output, file=(paste("high","_","low","_2_thirds_train_logistic_sig_",namesString,".csv", sep="")), row.names=FALSE)
bonferroni_pvalue<-paste(c("bonferroni p-value = 0.05/",length(genes),",",0.05/length(genes), sep=","))
cat(bonferroni_pvalue, file=(paste("high","_","low","_2_thirds_train_logistic_sig_",namesString,".csv", sep="")), append=TRUE)



#fisherTestsGenes<-sapply(seq(1:dim(trainMutationMatrixLogical)[2]), function(x) fisher.test(table(trainOutcome, trainMutationMatrixLogical[,x]), alternative="two.sided")$p.value)
#significantTests<-which(fisherTestsGenes<0.05/dim(trainMutationMatrixLogical)[2])
#trainMutationMatrixLogical<-trainMutationMatrixLogical[,significantTests]# train dataset for significant F.tests oin training set
#testMutationMatrixLogical<- testMutationMatrixLogical[,significantTests]# test dataset for significant F tests

colIndex_sig_genes<-significant_allGLMPvals+(length(covariates)+1)# col index of the genes in the mutation matrices
trainMutationMatrixLogical2<- trainMutationMatrixLogical[,c(seq(1,length(covariates)+1), colIndex_sig_genes)]
testMutationMatrixLogical2<-testMutationMatrixLogical[,c(seq(1,length(covariates)+1),colIndex_sig_genes)]


}
  


#####################################################################################################################
## add age and gender covariates to the train/testmutationMatrixLogical
#testMutationMatrixLogical<- as.data.frame(cbind(age=clinicalData2$age_at_initial_pathologic_diagnosis[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],testMutationMatrixLogical), stringsAsFactors=TRUE)
#testMutationMatrixLogical<- as.data.frame(cbind(gender=clinicalData2$gender[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],testMutationMatrixLogical), stringsAsFactors=TRUE)
#testMutationMatrixLogical<-as.data.frame(cbind(stage=clinicalData2$stage[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],testMutationMatrixLogical), stringsAsFactors=TRUE)
#testMutationMatrixLogical<- as.data.frame(cbind(outcome=testOutcome,testMutationMatrixLogical))



#trainMutationMatrixLogical<- as.data.frame(cbind(age=clinicalData2$age_at_initial_pathologic_diagnosis[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))
#trainMutationMatrixLogical<- as.data.frame(cbind(gender=clinicalData2$gender[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))
#trainMutationMatrixLogical<- as.data.frame(cbind(stage=clinicalData2$stage[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))
#trainMutationMatrixLogical<- as.data.frame(cbind(outcome=trainOutcome,trainMutationMatrixLogical))

###no stage information
#trainMutationMatrixLogical2<- trainMutationMatrixLogical[,c(1,2,3,4,colIndex_sig_genes)]
#testMutationMatrixLogical2<-testMutationMatrixLogical[,c(1,2,3,4,colIndex_sig_genes)]

###just the covariates
#trainMutationMatrixLogical2<- trainMutationMatrixLogical[,c(1,2,3,4)]
#testMutationMatrixLogical2<-testMutationMatrixLogical[,c(1,2,3,4)]

###just genetics
#trainMutationMatrixLogical2<- trainMutationMatrixLogical[,c(1,colIndex_sig_genes)]
#testMutationMatrixLogical2<-testMutationMatrixLogical[,c(1,colIndex_sig_genes)]

###just stage
#trainMutationMatrixLogical2<- trainMutationMatrixLogical[,c(1,2)]
#testMutationMatrixLogical2<-testMutationMatrixLogical[,c(1,2)]

#####################################################################################################################
## Train accuracy using my significant predictors


##recode stagelevels to be low high

if("stage" %in% covs && "stage" %in% colnames(trainMutationMatrixLogical2) ){
trainMutationMatrixLogical2$stage<-factor(trainMutationMatrixLogical2$stage, levels=c("low_stage", "high_stage"))
testMutationMatrixLogical2$stage<-factor(testMutationMatrixLogical2$stage, levels=c("low_stage", "high_stage"))
}

#recode grade outcome levels to be low high
trainMutationMatrixLogical2$outcome<-factor(trainMutationMatrixLogical2$outcome, levels=c("low", "high"))
testMutationMatrixLogical2$outcome<-factor(testMutationMatrixLogical2$outcome, levels=c("low", "high"))



train_GLM<- glm(outcome~., data=trainMutationMatrixLogical2, family="binomial")
train_GLM_summary<- summary(train_GLM)
train_GLM_OR<-exp(cbind(OR = coef(train_GLM), confint(train_GLM)))

sink(paste("train_GLM_stats_", namesString,".txt"))
print(train_GLM_summary)
print(train_GLM_OR)
sink()


#modelfit

dev<-with(train_GLM, null.deviance-deviance)# difference in deviance between the null model and our model
dfs<-with(train_GLM, df.null - df.residual)# df for the difference between the two models = difference between the number of predictors

train_GLM_modelfit<-with(train_GLM, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))# pvalue for a chisquare comparing difference between the deviance statistic for the null modela nd our model.
loglik<-logLik(train_GLM)

modelFit<- list(deviance_diff=dev,degrees_of_freedom=dfs,chisquare=train_GLM_modelfit, loglikelihood_of_model=loglik)
sink(paste("train_GLM_model_",namesString,".txt"), append =TRUE)
print(train_GLM_summary)
print(modelFit)
sink()

####################################################################################################################
  ###############model fit assessment using AIC ######################################################################
##The model after AIC backwards elimination should be better

train_GLM_step<- step(train_GLM,direction="backward", k=2)
train_GLM_step_OR<-exp(cbind(OR = coef(train_GLM_step), confint(train_GLM_step)))

sink(paste("train_GLM_step_stats_",namesString,".txt"))
print(summary(train_GLM_step))
print(train_GLM_step_OR)
sink()


##backwards elimination and model fitting###################################

dev2<-with(train_GLM_step, null.deviance-deviance)# difference in deviance between the null model and our model
dfs2<-with(train_GLM_step, df.null - df.residual)# df for the difference between the two models = difference between the number of predictors

train_GLM_step_modelfit<-with(train_GLM_step, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))# pvalue for a chisquare comparing difference between the deviance statistic for the null modela nd our model.
loglik2<-logLik(train_GLM_step)

modelFit2<- list(deviance_diff=dev2,degrees_of_freedom=dfs2,chisquare=train_GLM_modelfit, loglikelihood_of_model=loglik2)
sink(paste("train_glm_step_model_", namesString,".txt"), append =TRUE)
print(summary(train_GLM_step))
print(modelFit2)
sink()






###################################################################
##

#test_GLM2<- predict(train_GLM, newdata=new.data, type="response")
# TP53 is predictive of tumour grade when adjusted for age, gender, adn tumor stage

#with(train_GLM, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
#####################################################################

#ROC
##############################################################
## ROC curves


library(pROC)
library(scales)

#true adeno class
trueClass<-testOutcome
trueClass<-factor(trueClass, levels=c("low","high"))


#trueClass<-revalue(trueClass,c("squamous"="other", "BLCA"="other","leukemia"="other", "adenocarcinoma"="other"))# re-assigning the factor labels to be "adenocarinoma" or "other"
# now trueClass is the true gold standard class. This should be used for building ROC curves
test_GLM<- predict(train_GLM, newdata=testMutationMatrixLogical2, type="response")

clssvec<- test_GLM



test<-roc(trueClass,clssvec)
sink(paste("AUC_train_GLM",namesString,".txt", sep=""), append=TRUE)# append ROC AUC results
print("ROC RESULTS")
print(test)
sink()


png(paste("ROC_","high","_vs_","low","_logistic_train_GLM_",namesString,".png", sep=""), width=1200, height=1200)
par(mar=c(10, 10, 10, 10))
plot(test, col="red", bg="transparent", main = paste("ROC curve ", "high"," vs ", "low" ), cex.axis=2.5, ann=FALSE, lwd=2)


dev.off()



#################for the stepwise model###########################################################

test_GLM_step<- predict(train_GLM_step, newdata=testMutationMatrixLogical2, type="response")# reassign test_GLM predictions to be based on the glm_step model
write.table(test_GLM_step, paste("ROC_data_",namesString,".txt"), sep="\t")


clssvec<- test_GLM_step



####################################
####sensitivity and specificity for the logistic regression model

true_pred_table<-table(trueClass, test_GLM>0.5)# contingency table for true class and predicted class
true_pred_table2<-true_pred_table# The rows are the true class and the columns are the predicted
true_pred_table2[1,1]<-true_pred_table[2,2]# in order for the confusion matrix function to work, the TP and FP indeces needed to be swapped
true_pred_table2[2,2]<- true_pred_table[1,1]
colnames(true_pred_table2)<-c("high","low")
rownames(true_pred_table2)<-c("high","low")


colnames(true_pred_table)<- rownames(true_pred_table)

confusionMatrix(true_pred_table2)

sink(paste("Test_GLM_step_classifier_stats",namesString,".txt", sep=""), append=TRUE)

print("High class is defined as a positive and Low class is a negative. index [1,1] is TP and index[2,2] is TN")
print(confusionMatrix(true_pred_table2))
sink()


################################################################

test2<-roc(trueClass,clssvec)
sink(paste("AUC_GLM_step",namesString,".txt", sep=""), append=TRUE)# append ROC AUC results
print("ROC RESULTS")
print(test2)
sink()





sink(paste("Test_","GLM_step_ROC_predictions_",namesString,".txt", sep=""), append=TRUE)# append ROC AUC results
print("ROC RESULTS")
print(test2$predictor)
sink()


png(paste("ROC_","high","_vs_","low","_logistic_train_GLM_step_",namesString,".png", sep=""), width=1200, height=1200)
par(mar=c(10, 10, 10, 10))
plot(test2, col="blue", bg="transparent", main = paste("ROC curve ", "high"," vs ", "low GLM step" ),cex.axis=2.5,ann=FALSE, lwd=2)


dev.off()


#######combined ROC curve###################################


AUC<-c(test$auc,test2$auc)
AUC<-round(AUC,3)

textForLegend<-cbind(c(paste(gsub("_", " ",namesString)),paste(gsub("_", " ",namesString)," AIC")),AUC)
legendText<-sapply(1:dim(textForLegend)[1], function(x) paste(textForLegend[x,1], textForLegend[x,2]))
legend(0.65,0.25, legendText, col=c(1,2), cex=2.5, lty=1, lwd=2)

png(paste("ROC_","high","_vs_","low_","initial_and_AIC_",namesString,".png", sep=""), width=1200, height=1200)
par(mar=c(10, 10, 10, 10))
plot(test, col=alpha("red", alpha=0.6), bg="transparent", main = paste("ROC curve ", "high"," vs ", "low GLM step" ),cex.axis=2.5,ann=FALSE, lwd=2)
lines(test2, col=alpha("blue", alpha=0.6),lwd=2)
AUC<-c(test$auc,test2$auc)
AUC<-round(AUC,3)


legendText<-sapply(1:dim(textForLegend)[1], function(x) paste(textForLegend[x,1], textForLegend[x,2]))
legend(0.65,0.25, legendText, col=c("red","blue"), cex=2, lty=1, lwd=2)


dev.off()

}


########Post hoc tests for bias in predictors according to cancer type

##get tumour type

trainTumourType<- clinicalData2[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds),"tumourType", drop=FALSE]
colnames(trainTumourType)<- "tumourType"
rownames(trainTumourType)<- clinicalData2[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds),"SampleIds"]


##chi-square to test if predictive protein mutations are biased among tumour types.


predictors<-colnames(trainMutationMatrixLogical2[-c(1:4)])

#chisquares for the protein predictors
predictorTumourTypeChisquares<-apply(trainMutationMatrixLogical2[,-c(1:4)], 2, function(x) chisq.test(x,trainTumourType[,1]))


#tables for the protein predictors and tumour type
predictorTumourTypeTables<- lapply(seq(5:(dim(trainMutationMatrixLogical2[,-c(1:4)])[2]+4)), function(x) table(trainTumourType[,1],trainMutationMatrixLogical2[,x+4]))

#####################################################################################################################
#check for multicollinearity

 #all pairwise comparisons between cancer classes
#pairs<- outer(4:7,4:7,paste, sep="_")
#comps<-pairs[upper.tri(pairs)]



#nComps<- length(comps)# iterate over this when I want to automate the comparisons I do using lapply or a loop on the cluster

# use chi-square to check for multicollinearity significant result indicates no
#multicollinearity_chisquares<-

#for (i in 1:length(comps)){

  
  #put all the code in here
 # twogenes<- as.numeric(unlist(strsplit(comps[6],split="_")))
#  gene1<- twogenes[1]
 # gene2<- twogenes[2]
  
  #table(trainMutationMatrixLogical2[,gene1], trainMutationMatrixLogical2[,gene2])
  #chisq.test(x=trainMutationMatrixLogical2[,gene1], y=trainMutationMatrixLogical2[,gene2], correct=TRUE)

#  table(testMutationMatrixLogical2[,gene1], testMutationMatrixLogical2[,gene2])
 # chisq.test(x=testMutationMatrixLogical2[,gene1], y=testMutationMatrixLogical2[,gene2], correct=TRUE)
  
#}



















#######RandomForest analysis########################################################################################


#####################################################################################################################
## add SNV_TYPE measures to the train/testmutationMatrixLogical
#testMutationMatrixLogical<-as.data.frame(cbind(SNV_Types_freq_samples[match(testSamples,rownames(SNV_Types_freq_samples)),],testMutationMatrixLogical), stringsAsFactors=TRUE)
#trainMutationMatrixLogical<-as.data.frame(cbind(SNV_Types_freq_samples[match(trainSamples,rownames(SNV_Types_freq_samples)),],trainMutationMatrixLogical), stringsAsFactors=TRUE)





#rfe control fuction
#set.seed(10)
#ctrl <- rfeControl(functions = rfFuncs, method = "repeatedcv", repeats = 5, verbose = FALSE, returnResamp = "final")

#subsets<-c(seq(1:50),60,70,80,90,100)

#############################################################################################################################
#Model Building
#trainModel<-rfe(trainMutationMatrixLogical2, trainOutcome, subsets, rfeControl=ctrl)


#  sink(paste("Train_","high","_","low","_2_thirds_randomForest_results.txt", sep=""))
#print(trainModel)
#print(trainModel$optVariables)
#sink()

#testClasses <- predict(trainModel, testMutationMatrixLogical)# the test run classes
#############  confusionMatrix(testClasses$pred,testOutcome)# the test run accuracy results
#sink(paste("ROC_curve_data_", "high","_","low",".txt", sep=""))
#print(testClasses)
#sink()


#sink(paste("Test_balanced_","high","_", "low","_2_thirds_randomForest_results.txt", sep=""), append=TRUE)
#print(testClasses)
#print(confusionMatrix(testClasses$pred,testOutcome))
#sink()



##############################################################
## ROC curves


#library(pROC)

##true adeno class
#trueClass<-testOutcome
###################trueClass<-revalue(trueClass,c("squamous"="other", "BLCA"="other","leukemia"="other", "adenocarcinoma"="other"))# re-assigning the factor labels to be "adenocarinoma" or "other"
# now trueClass is the true gold standard class. This should be used for building ROC curves


#clssvec<- testClasses[,"high"]


#test<-roc(trueClass,clssvec)
#sink(paste("Test_balanced_","high","_", "low","_2_thirds_randomForest_results.txt", sep=""), append=TRUE)# append ROC AUC results
#print("ROC RESULTS")
#print(test)
#sink()


#png(paste("ROC_","high","_vs_","low",".png", sep=""), width=800, height=800)
#par(mar=c(10, 10, 10, 10))
#plot(test, col="red", bg="transparent", main = paste("ROC curve ", "high"," vs ", "low" ))


#dev.off()

####ROCR ROC curve#########################################


#tst<-mapvalues(trueClass, c(toString("high"), toString("low")), c(1,0))
#library(ROCR)
#pred<- prediction(clssvec,tst)
#perf<- performance(pred,"tpr", "fpr")
#png(paste("ROCR_ROC_","high","_vs_","low",".png", sep=""), width=800, height=800)
#par(mar=c(10, 10, 10, 10))
#plot(perf, colorize=TRUE, lwd=3, main = paste("ROC curve ","high"," vs ","low", sep=""), bg="transparent")
#dev.off()
#}






