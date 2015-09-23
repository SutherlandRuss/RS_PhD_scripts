#defines the library paths within which the libraries are found
.libPaths(c("/Users/rds/Documents/R/win-library/3.0", "/Program Files/R/R-3.0.2/library"))

library(plyr)
library(vegan)
library(sets)
#library(sampling)
library(scales)
library(DMwR)


#install.packages("sampling")
#
# Set the filepath
dir.create("/Users/rds/Dropbox/PhD/CUP_classifier/outputs/11_class_balanced/10_class_test")
setwd("/Users/rds/Dropbox/PhD/CUP_classifier/outputs/11_class_balanced/10_class_test")

#setwd("/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/5_class/sex_gene_filtered/170214/RRF")
#setwd("/Users/rsutherland/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/5_class/sex_gene_filtered")




########Parameters to set before starting to run the program#######################

# choose to balance the outcome classes through downsampling of the majority class
balance<-TRUE

#bonferroni corrected or ranked

falseDR<-TRUE

# Use SNV_type predictors

SNV_Types<-TRUE

# do I just want to use SNV_Type predictors

SNV_Types_And_clinical<- FALSE

# do I just want clinical data only?

clinicalOnly<- FALSE


# set the random seed

set.seed(10)

#
#
#
# The data file in vcf-like format.
#scores.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/sep_2013/input"
scores.path <- "/Users/rds/Dropbox/PhD/tumour_classifier_data/sep_2013/input"
#scores.path <- "/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/sep_2013/input"

#scores.files<- "pancan12_cleaned.maf"


#Tissue source site file
TSS.path<-"/Users/rds/Dropbox/PhD/tumour_classifier_data/sep_2013/TSS"
#TSS.path<-"/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/sep_2013/TSS"
TSS.file<- "tissueSourceSite.csv"

######################################################################################################################
### metadata
#########

clinical.path<-"/Users/rds/Dropbox/PhD/tumour_classifier_data/SynapsePanCancer/PanCan12/clinical" # load all clinical data at once
#clinical.path<-"/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/SynapsePanCancer/PanCan12/sex_specific_removed/clinical"# load all clinical data at once


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

clinicalData.cancers<-lapply(seq(1:length(clinicalData.directories)), function(x) strsplit(clinicalData.directories[x],"/",)[[1]][10])



createClinicalData<- function(clinicalData.directories, clinicalData.cancers){
  
  clinicals<-lapply(seq(1:length(clinicalData.cancers)),function(x) clinicalDataTable(clinicalData.directories[x],strsplit(clinicalData.directories[x],"/",)[[1]][10]))
  names(clinicals)<- clinicalData.cancers
  return(clinicals)
}

allClinicalData<-createClinicalData(clinicalData.directories, clinicalData.cancers)# loads all of the clinical data

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


####################################################################################
## AML is not a solid tumour remove it

clinicalData<- clinicalData[-which(clinicalData$tumourType=="LAML"),]# is not a solid tumour



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


#important variabls to include
#clinicalData2$stage##not available for BLCA
clinicalData2$age_at_initial_pathologic_diagnosis
clinicalData2$gender


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


#################################################################################################################################################

####################################################################################################################################################


#######################################################################################################################################
scores.SNV<- scores[which(scores$Variant_Type=="SNP"),]

#swapStrand<- function(negativeVariant){
#  
#  v<- negativeVariant
#  
#  if(v=="A"){
#    v<-"T"
#  }else if(v=="T"){
#    v<-"A"
#  }else if(v=="C"){
#    v<- "G"
#  }else if(v=="G"){
#    v<-"C"
#  }
#  return(v)
#}#
#
#
#swapNegativeStrand<- function(scoreTable){
#  
#  # negative strand index
#  indexNegative<-which(scoreTable$strand=="-1")
#  scoreTable.negative<-scoreTable[indexNegative,]
#  
#  
#  scoreTable.Negative.variant<-scoreTable.negative$variant
#  
#  swapped<-sapply(scoreTable.Negative.variant, swapStrand, USE.NAMES=FALSE)
#  scoreTable$variant[indexNegative]<- swapped
#  
#  return(scoreTable)
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



SNV_Type_freqs<-function(scoresDataFrame){
  
  AC_TG<- getSNVType(scoresDataFrame, "A","C","T","G")
  AG_TC<- getSNVType(scoresDataFrame, "A","G","T","C")
  AT_TA<- getSNVType(scoresDataFrame, "A","T","T","A")
  CA_GT<- getSNVType(scoresDataFrame, "C","A","G","T")
  CG_GC<- getSNVType(scoresDataFrame, "C","G","G","C")
  CT_GA<- getSNVType(scoresDataFrame, "C","T","G","A")
  
  SNV_Types_freq<- c(AC_TG,AG_TC, AT_TA, CA_GT, CG_GC, CT_GA)
  names(SNV_Types_freq)<- c("A>C_T>G","A>G_T>C", "A>T_T>A", "C>A_G>T", "C>G_G>C", "C>T_G>A")
  
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

barplot(colMeans(mutationTypes))

###########
# VariantFreqTable

variantFreqTable<- cbind(SNV_Types_freq_samples, mutationTypes)


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

outcome<-clinicalData2$tumourType
outcome[which(outcome=="COAD"| outcome=="READ")]<-"COADREAD"

outcome<- as.factor(outcome)#The classes to be used for the randomforest I will run. 
names(outcome)<- colnames(mutationMatrixLogical)
write.table(file="tumourclass_type_contingency", table(cancer_Type1,outcome), row.names=TRUE, col.names=TRUE, sep="\t")

#I can use this to extract the appropriate variables from the MutationMatrixLogical before I calculate F.tests or SNV_Type frequencies.


#gender table
genderTable<- table(outcome,clinicalData2$gender)
XsqTest<-chisq.test(genderTable)


#age table
ageDataFrame <- as.data.frame(cbind(clinicalData2$age_at_initial_pathologic_diagnosis,outcome))
ageAOV<- aov(clinicalData2$age_at_initial_pathologic_diagnosis~outcome)
ageAOVANOVA<-anova(ageAOV)
pairwise.t.test(clinicalData2$age_at_initial_pathologic_diagnosis, outcome, p.adjust="none")
# there are differences between the tumour types in terms of age

#all pairwise comparisons between cancer classes
pairs<- outer(subTypes,subTypes,paste, sep="_")
comps<-pairs[upper.tri(pairs)]

############

#distance1<-dist(mutationMatrixBinary, method="manhattan", diag=FALSE, upper=FALSE)
#distance2<- vegdist(mutationMatrixBinary, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE)
#mds1<-cmdscale(distance2, k=2, eig=TRUE, add=FALSE, x.ret=FALSE)# object too large for memory
#pcAnalysis<-prcomp(t(mutationMatrixBinary), scale=TRUE, center=TRUE)
#pcAnalysis2<-prcomp(mutationMatrixBinary, scale=TRUE, center=TRUE)
#plot(pcAnalysis$x[,1], pcAnalysis$x[,2])
#plot(pcAnalysis2$rotation[,1], pcAnalysis2$rotation[,2])

nComps<- length(comps)# iterate over this when I want to automate the comparisons I do using lapply or a loop on the cluster
AUCs<- rep(0, length(comps))


##############################################################################
###identify minority class


mutationMatrixLogicalA<-t(mutationMatrixLogical)
#########################################################################################################################
#########################################################################################################################
## add SNV_TYPE measures to the train/testmutationMatrixLogical
if(SNV_Types_And_clinical==TRUE){
  mutationMatrixLogicalB<-variantFreqTable
}else if (SNV_Types==TRUE){
  mutationMatrixLogicalB<-as.data.frame(cbind(variantFreqTable, mutationMatrixLogicalA), stringsAsFactors=TRUE)
}


#####################################################################################################################
## add age and gender covariates to the train/testmutationMatrixLogical
mutationMatrixLogicalB<- as.data.frame(cbind(age=clinicalData2$age_at_initial_pathologic_diagnosis,mutationMatrixLogicalB), stringsAsFactors=TRUE)
mutationMatrixLogicalB<- as.data.frame(cbind(gender=clinicalData2$gender,mutationMatrixLogicalB), stringsAsFactors=TRUE)
#mutationMatrixLogicalB<- as.data.frame(cbind(tumourType=clinicalData2$tumourType,mutationMatrixLogicalB), stringsAsFactors=TRUE)
#mutationMatrixLogicalB<- as.data.frame(cbind(stage=clinicalData2$stage,mutationMatrixLogicalB), stringsAsFactors=TRUE)
#mutationMatrixLogicalB<- as.data.frame(cbind(grade=clinicalData2$grade,mutationMatrixLogicalB), stringsAsFactors=TRUE)
mutationMatrixLogicalB<-cbind(outcome, mutationMatrixLogicalB)

# number of non protein variables if including outcome, gender, age, tumourtype, stage, grade
nNonProteinVars<- 3+(dim(variantFreqTable)[2])



## if I only want to use clinical Data this statement changes the train and test data to just use age and gender information
#if(clinicalOnly==TRUE){
#  trainMutationMatrixLogical2<- trainMutationMatrixLogical[,1:2]
#  testMutationMatrixLogical2<- testMutationMatrixLogical[,1:2]
#}else if (SNV_Types_And_clinical==TRUE){
#  trainMutationMatrixLogical2<-trainMutationMatrixLogical[,-1]
#  testMutationMatrixLogical2<-trainMutationMatrixLogical[,-1]
#}


#library(DMwR)

#trainsmote<-SMOTE(outcome~., data=mutationMatrixLogicalB)





if(balance){
  #identify the minority class and balance
  
  # minority class training set size
  minclass23rds<- floor(((min(table(outcome))/3)*2))
  
  
}

# ---------------------------------------------------------------------------
random_sample = function (dt, sample_size) {
  # determine the number of records in the data frame
  length_dt = length (dt[,1])
  # extract the sample
  dt [sample (1:length_dt, size = sample_size),]
}
# ---------------------------------------------------------------------------



dt<-mutationMatrixLogicalB
strat_by<-"outcome"
sample_intensity<-0.66
get_stratified_sample = function (dt, strat_by, sample_intensity=0.1) {
  # dt contains the data
  # strat_by contains the name of the variable according to which the sample is to be stratified
  # sample intensity is proportion of data to be sampled, default = 10%
  # deterine the number of levels of factor 'strat_by' and store the number of rows relevant to the factor in table tmp
  f = factor(dt[,strat_by])
  tmp = aggregate (dt[,strat_by], by=list(f), FUN=length)
  #
  # determine the number of rows to be included within the sample for each level of the factor
  #the matching gets the index of the minority class
  tmp$sz = round(tmp[match(min(tmp[,2]), tmp[,2]),2] * sample_intensity, 0)
  #
  f.number = length (tmp [,1])      # determine the number of levels of the factor
  #
  # obtain the sample for the first level and store in strat_sample
  strat_sample = random_sample (subset (dt, f == tmp[1,1]), tmp[1,3])
  # obtain subsample for the remaining levels of the factor.
  if (f.number > 1) {
    for (i in 2:f.number) { 
      rs = random_sample (subset (dt, f == tmp[i,1]), tmp[i,3])
      strat_sample = rbind (strat_sample, rs)
    }
  }
  return(strat_sample)      
}




mutationMatrixLogicalTrain<-get_stratified_sample(mutationMatrixLogicalB,"outcome", 0.66)

#sample IDS that are not in the Train Set.
testSamples<-setdiff(rownames(mutationMatrixLogicalB), rownames(mutationMatrixLogicalTrain))
trainSamples<- intersect(rownames(mutationMatrixLogicalTrain), rownames(mutationMatrixLogicalB))  

#test MutationMatrixLogical samples
mutationMatrixLogicalTest<-mutationMatrixLogicalB[testSamples,]


#test and train outcome

trainOutcome<-outcome[trainSamples]
testOutcome<-outcome[testSamples]




###################################BEGINE here tomorrow

## filter out proteins according to the number of samples in which they carry mutations ################################
tumourType<-unique(outcome)

filterMuts<- function(mutationMatrix){
  #filter out the genes that do not carry mutations in more than 5% of the samples
  #5% of samples is 970/20
  
  FivePercent<- ceiling(dim(mutationMatrix)[1]/20)
  mutationMatrix_proteins<-names(which(colSums(mutationMatrix)>FivePercent))
  
  return(mutationMatrix_proteins)
}

##########

train_ClinicalData2<- clinicalData2[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds),]
#the names of the proteins carried through to the next stage, by cancer type
filtered_protein_Names<-lapply(tumourType, function(x) filterMuts(mutationMatrixLogicalTrain[train_ClinicalData2$tumourType==x,-c(1:nNonProteinVars)]))
names(filtered_protein_Names)<- tumourType

# the unique list of proteins carrying mutations in at least 5% of at least cancer type
proteins_analysis<- unique(unlist(filtered_protein_Names))




##redefining the train and testMutationMatrices.
# clinical metadata
trainMutationMatrixLogical<-mutationMatrixLogicalTrain[,c(1:nNonProteinVars)]
testMutationMatrixLogical<-mutationMatrixLogicalTest[,c(1:nNonProteinVars)]


trainMutationMatrixLogical<- as.data.frame(cbind(trainMutationMatrixLogical,mutationMatrixLogicalTrain[,proteins_analysis]))
testMutationMatrixLogical<-as.data.frame(cbind(testMutationMatrixLogical,mutationMatrixLogicalTest[,proteins_analysis]))

#################################################################################################################################
################### CARET random forest

##Defining a training and test set for 
library(caret)
library(mlbench)
library(Hmisc)
library(randomForest)
library(e1071)



## define the sizetolerance function
rfFuncs2<-rfFuncs

# selectModeSize function
rfFuncs2$selectSize<-function (x, metric="Accuracy", tol=1, maximize) 
{
  if (!maximize) {
    best <- min(x[, metric])
    perf <- (x[, metric] - best)/best * 100
    flag <- perf <= tol  
  }
  else {
    best <- max(x[, metric])
    perf <- (best - x[, metric])/best * 100
    flag <- perf <= tol
  }
  min(x[flag, "Variables"])
}


rfFuncs2$selectVar<-pickVars


set.seed(10)

#rfe control fuction
set.seed(10)
ctrl <- rfeControl(functions = rfFuncs2, method = "repeatedcv", repeats = 5, verbose = FALSE, returnResamp = "final", allowParallel=TRUE, rerank=FALSE)

subsets<-c(seq(1:50),60,70,80,90,100,110,120,130,140,150,160,170,180,190,200)
#subsets<-c(10,20,30,40,50)


#############################################################################################################################
#Model Building
trainModel<-rfe(trainMutationMatrixLogical[,-1], trainOutcome, subsets, rfeControl=ctrl, proximity=TRUE, localImp=TRUE)


# local importance and proximity are used for creating heatmaps

#####################################################################################################################
##### model selection and accuracy plots



#plot how each of the size tolerant models was arrived at
picksizeResults<- data.frame(Accuracy=trainModel$results$Accuracy, Variables=seq(1:dim(trainModel$results)[1]))


# write the acuracies and tolerances to the results folder
write.table(picksizeResults, file=paste("tolerance_results.txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

#in order to calculate tolerance I need the maximum accuracy from the list of models
maxAcc <- max(picksizeResults$Accuracy)
picksizeResults$Tolerance <- 100-(abs((picksizeResults$Accuracy / maxAcc))*100)


within1pcnt<- pickSizeTolerance(picksizeResults, metric="Accuracy", tol=1, maximize=TRUE)
smallest<- pickSizeBest(picksizeResults, metric="Accuracy", maximize=TRUE)
#png(paste("size_tolerance_plots_",subType1,"_",subType2,".png"), 1200,1000, antialias="cleartype")
pdf(paste("size_tolerance_plots.pdf"), 12,10, pointsize=12)
par( mar = c(7, 7, 7, 15), mgp=c(4,1,0))
plot(picksizeResults$Variables[-c(smallest, within1pcnt)],
     picksizeResults$Accuracy[-c(smallest, within1pcnt)],
     ylim = extendrange(picksizeResults$Accuracy),
     ylab = "Accuracy", xlab = "Variables",cex=1.8,cex.axis=2, cex.lab=2)
points(picksizeResults$Variables[smallest],
       picksizeResults$Accuracy[smallest], pch = 16, cex= 2.3)
points(picksizeResults$Variables[within1pcnt],
       picksizeResults$Accuracy[within1pcnt], pch = 17, cex= 2.3)
par(xpd=TRUE)
legend(dim(trainModel$results)[1]+2, max(trainModel$results$Accuracy), c("Maximum\naccuracy",paste("\nSelected\nModel\n",within1pcnt, " features",sep="")), pch=c(16,17), cex=2.0, bty="n",y.intersp=2)
# now I do not need the second plot
#with(picksizeResults, plot(Variables, Tolerance, cex=1.8, cex.lab=2,cex.axis=2))
#points(picksizeResults$Variables[smallest],
#       picksizeResults$Tolerance[smallest], pch = 16, cex= 2.3)
#points(picksizeResults$Variables[optimalModel],
#       picksizeResults$Tolerance[optimalModel], pch = 17, cex= 2.3)
#legend(60, 12, c("Best\nmodel",paste("\nSelected\nModel\n ",optimalModel, " features", sep="")), pch=c(16,17), cex=2.0, bty="n",y.intersp=2)
#   par(xpd=FALSE)
# abline(h = sizetol[i], lty = 2, col = "darkgrey")
dev.off()

# plot the distribution of the accuracy statistic

pdf(paste("accuracy_hist_distribution.pdf",sep=""),12,10, pointsize=12)
hist(picksizeResults[,1],xlab="Accuracy", main="", border=NA, col="lightgray", cex=2)
#curve(dnorm(x, mean=mean(picksizeResults[,1]), sd=sd(picksizeResults[,1])), add=TRUE, col=col.plot[2], lwd=2)
lines(density(picksizeResults[,1],adjust=2), col=col.plot[2],lwd=2)
dev.off()




##### variable importance parameters

trainModelImportances<- importance(trainModel$fit, type="1")
trainModelImportances<- trainModelImportances[order(trainModelImportances, decreasing=TRUE),]

#print the model parameters     
sink(paste("Train_10class_all_genes_randomForest_results.txt", sep=""))
print(trainModel$optsize)
print("trainmodel classes, the first factor is the positive class")
print(trainModel$fit$classes)
print(trainModel)
print(trainModel$optVariables)
print(cbind(trainModelImportances))
sink()

# optimal model variables for kegg pathway enrichment
write.table(trainModel$optVariables,paste("kegg_protein_list.txt", sep=""), row.names=FALSE, col.names=FALSE )



#sink(paste("Train_2_thirds_randomForest_results.txt", sep=""))
#print(trainModel)
#print(trainModel$optVariables)
#print(varImp(trainModel))
#sink()







########plots and model selection

testClasses <- predict(trainModel, testMutationMatrixLogical[,-1])# the test run classes
#  confusionMatrix(testClasses$pred,testOutcome)# the test run accuracy results
sink(paste("ROC_curve_data",".txt", sep=""))
print(testClasses)
sink()


classificationStats<- confusionMatrix(testClasses$pred, testOutcome)
sink(paste("classification_stats_confusion_matrix.txt", sep="\t"))
print(classificationStats)
sink()



#True positive classifications in 5 class random forest
#TP<-sapply(seq(1:length(classStatsSubTypes)), function(x) classificationStats$table[x,x])

#FP rate
#FP<-



###############################################################################
library(pROC)
library(colorspace)


# vector of classes
classStatsSubTypes<-colnames(classificationStats$table)

#True classes for ROC
trueClass<- testOutcome

AUCs<- rep(0, length(classStatsSubTypes))

for(i in seq(1:length(classStatsSubTypes))){
  
  #predicted probability of classification as the cancer subtype under investigation
  predprobabilities<- testClasses[,i+1]
  
  trueClassBinary<-as.character.factor(trueClass)
  trueClassBinary[which(trueClass!=classStatsSubTypes[i])]<-"other"
  trueClassBinary<- as.factor(trueClassBinary)
  
  test<-roc(trueClassBinary,predprobabilities)
  AUCs[i]<-test$auc
  
  sink(paste("Test_2_thirds_5_class_randomForest_results.txt", sep=""), append=TRUE)# append ROC AUC results
  print(paste("ROC RESULTS",classStatsSubTypes[i]))
  print(test)
  sink()
  
  
  #col.rainbow<- col.rainbow<- rainbow_hcl(10, c=70, l=50, start = 30, end= 300)
  col.plot<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99")
  palette(col.plot)
  if(i==1){
    
    png("10_classifier_ROC.png", width=1200, height=1200)
    par(mar=c(5,5,5,25))
    par(xpd=T)
    plot(test, xlim=c(1,0), ylim=c(0,1), col=alpha(i,alpha=0.8),lwd=8, ann=FALSE, cex.axis=2, bg="transparent",bty="L" )
  }else{
    
    lines(test, col=alpha(i,0.5), lwd=8)
  }
}
AUCs<-round(AUCs,2)
textForLegend<-cbind(classStatsSubTypes,AUCs)
legendText<-sapply(1:dim(textForLegend)[1], function(x) paste(textForLegend[x,1],"=", textForLegend[x,2]))

legend(0.40,0.45,legendText,title="AUC",col=c(1:10), cex=2.4, lty=1, lwd=8, bty="n")
dev.off()

## save the random forest object

tenClassRF_all_features2<-trainModel
save(tenClassRF_all_features2, file="tenClassRF_all_features.RData")
print(tenClassRF_all_features2)

