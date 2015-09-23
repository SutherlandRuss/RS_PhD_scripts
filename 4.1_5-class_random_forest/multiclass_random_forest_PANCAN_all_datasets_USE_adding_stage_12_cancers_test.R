#defines the library paths within which the libraries are found
.libPaths(c("C:/Users/rds/Documents/R/win-library/3.0", "C:/Program Files/R/R-3.0.2/library"))

library(plyr)
library(vegan)
library(sets)
library(scales)
library(fmsb)
library(changepoint)
library(stringr)
library(gtools)
#
# Set the filepath
dir.create("C:/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/5_class/5_class_pairwise/test/")
setwd("C:/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/5_class/5_class_pairwise/test/")
#setwd("C:/Users/rsutherland/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/5_class/sex_gene_filtered")



########Parameters to set before starting to run the program#######################

# choose to balance the outcome classes through downsampling of the majority class
balance<-TRUE

#bonferroni corrected or ranked

falseDR<-TRUE

# Use SNV_type predictors

SNV_Types<-FALSE

# do I just want to use SNV_Type predictors

SNV_Types_And_clinical<- FALSE

# do I just want clinical data only?

clinicalOnly<- FALSE

## do not use gender specific cancers?

removeGenderSpec<- FALSE


#covariates
covs<-c("age","gender","stage","outcome")

########size tolerant random forest

sizeTolerant1pcnt<- TRUE

# set the random seed

set.seed(10)

#
#
#
# The data file in vcf-like format.
#scores.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/sep_2013/input"
scores.path <- "C:/Users/rds/Dropbox/PhD/tumour_classifier_data/sep_2013/input"
#scores.path <- "C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/sep_2013/input"

#scores.files<- "pancan12_cleaned.maf"


#Tissue source site file
TSS.path<-"C:/Users/rds/Dropbox/PhD/tumour_classifier_data/sep_2013/TSS"
#TSS.path<-"C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/sep_2013/TSS"
TSS.file<- "tissueSourceSite.csv"

######################################################################################################################
### metadata
#########
if (removeGenderSpec){
  clinical.path<-"C:/Users/rds/Dropbox/PhD/tumour_classifier_data/SynapsePanCancer/PanCan12/sex_specific_removed/clinical"# load all clinical data at once
}else{
clinical.path<-"C:/Users/rds/Dropbox/PhD/tumour_classifier_data/SynapsePanCancer/PanCan12/clinical"# load all clinical data at once
#clinical.path<-"C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/SynapsePanCancer/PanCan12/sex_specific_removed/clinical"# load all clinical data at once
}

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
clinicalData.cancers.index<-length(unlist(strsplit(clinicalData.directories[1],"/")))# the index of the cancer names in the clinicalData.directories

clinicalData.cancers<-lapply(seq(1:length(clinicalData.directories)), function(x) strsplit(clinicalData.directories[x],"/",)[[1]][clinicalData.cancers.index])



createClinicalData<- function(clinicalData.directories, clinicalData.cancers){
  
  clinicals<-lapply(seq(1:length(clinicalData.cancers)),function(x) clinicalDataTable(clinicalData.directories[x],strsplit(clinicalData.directories[x],"/",)[[1]][clinicalData.cancers.index]))
  names(clinicals)<- clinicalData.cancers
  return(clinicals)
}

allClinicalData<-createClinicalData(clinicalData.directories, clinicalData.cancers)# loads all of the clinical data



#########trying to add the stage information here
if(removeGenderSpec==FALSE){
colnames(allClinicalData$BLCA)[which(colnames(allClinicalData$BLCA)=="ajcc_neoplasm_disease_stage")]<- "tumor_stage"
colnames(allClinicalData$BRCA)[which(colnames(allClinicalData$BRCA)=="ajcc_neoplasm_disease_stage")]<- "tumor_stage"
colnames(allClinicalData$UCEC)[which(colnames(allClinicalData$UCEC)=="gynecologic_tumor_grouping_figo_stage")]<- "tumor_stage"
}else{
  colnames(allClinicalData$BLCA)[which(colnames(allClinicalData$BLCA)=="ajcc_neoplasm_disease_stage")]<- "tumor_stage"
}

#create dummy variables for GBM and LAML cancers
allClinicalData$LAML<-cbind(allClinicalData$LAML,tumor_stage=rep(0,dim(allClinicalData$LAML)[1]))
allClinicalData$GBM<-cbind(allClinicalData$GBM,tumor_stage=rep(0,dim(allClinicalData$GBM)[1]))


# list of all clinical data variables

allClinicalNames<-lapply(seq(1:length(allClinicalData )), function(x) colnames(allClinicalData[[x]]))

#getStage<-sapply(allClinicalNames, function(x) match("tumor_stage",x))




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


########recode stage
#remove stage classes that cannot be grouped

# remove tumours without grade information
clinicalData<-clinicalData[-which(is.na(clinicalData$tumor_stage)),]
if(removeGenderSpec==FALSE){
clinicalData<-clinicalData[-which(clinicalData$tumor_stage=="Stage X"|clinicalData$tumor_stage=="Stage Tis"),]# removing the tumours that could not be graded
}


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
stage<- as.numeric(stage_v)
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


#important variabls to include
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


#################################################################################################################################################

####################################################################################################################################################


#######################################################################################################################################
scores.SNV<- scores[which(scores$Variant_Type=="SNP"),]

#scores.SNV by sample
scores.SNV.samples<-lapply(samplesAnalysis, function(x) scores.SNV[which(scores.SNV$sampleID==x),])



####Function to get SNV type frequency across the entire dataset
getSNVType<- function(scoresDataFrame, b1,b2,b3,b4){
  substitution1<- length(which(scoresDataFrame$variant== b1 & scoresDataFrame$reference== b2))
  substitution2<- length(which(scoresDataFrame$variant== b3 & scoresDataFrame$reference== b4))
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
outcome<- rep(0,length(cancer_Type1))


subTypes<- c("adenocarcinoma","squamous","urothelial", "GBM","leukemia")
if(removeGenderSpec){
groups.subTypes<- list(adenocarcinoma=groups[c(2,4,5,7)],squamous=groups[c(3,8)], urothelial=groups[c(9)], GBM=groups[1], leukemia=groups[6])
#classes<-c("adenocarcinoma", "BLCA")# block out for all class
}else{
groups.subTypes<- list(adenocarcinoma=groups[c(2,3,5,6,7,8,10)], squamous=groups[c(4,11)], urothelial=groups[12], GBM=groups[1], leukemia=groups[9])# block out for all class
}
classes<- list(classes=subTypes, subtype_names=groups.subTypes)# a list of the cancer cell type groups and the cancer subtypes that belong to them

# re-assigns a 0 class label to 1 depedent on a cancer_Type label matching one ine the groups.adenocarcinomas vector.

outcomeClass<- function(outcome,classes){
  for(j in seq_along(classes[[1]])){
    for(i in seq_along(classes[[2]][[j]])){
      
      outcome[which(cancer_Type1==classes[[2]][[j]][i])]<- classes[[1]][j]
      names(outcome)<- names(cancer_Type1)
    }
  }
  return(outcome)
}

outcome<-outcomeClass(outcome,classes)
write.table(file="tumourclass_type_contingency", table(cancer_Type1,outcome), row.names=TRUE, col.names=TRUE, sep="\t")
outcome<- as.factor(outcome)#The classes to be used for the randomforest I will run. 
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

##Leave out the stage variable from comparisons 4-10 because they include GBM or Leukeamia, both of which do not have stage information
leaveStage<- seq(4,10)


############

#distance1<-dist(mutationMatrixBinary, method="manhattan", diag=FALSE, upper=FALSE)
#distance2<- vegdist(mutationMatrixBinary, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE)
#mds1<-cmdscale(distance2, k=2, eig=TRUE, add=FALSE, x.ret=FALSE)# object too large for memory
#pcAnalysis<-prcomp(t(mutationMatrixBinary), scale=TRUE, center=TRUE)
#pcAnalysis2<-prcomp(mutationMatrixBinary, scale=TRUE, center=TRUE)
#plot(pcAnalysis$x[,1], pcAnalysis$x[,2])
#plot(pcAnalysis2$rotation[,1], pcAnalysis2$rotation[,2])

nComps<- length(comps)# iterate over this when I want to automate the comparisons I do using lapply or a loop on the cluster
AUCs<- rep(0, nComps)# vector of AUC values to fill
cancer_subtypes<- matrix(data=NA, nrow=nComps, ncol=2)# cancer subtype matrix for legend labels in ROC curves
ROC_objects<- list()# list to hold ROC objects for use in the ROC curve
RegModels<- list()# list to hold the regression model outputs
classStats<-list()# list to hold the models classificfation stats
subsetmodelsize<-list()# output from change analysis
subsetmodelsizeselect<-list()# vector of model sizes from the breakpoint analysis
for (i in 1:nComps){
i<-9
  #put all the code in here
  twoSubTypes<- unlist(strsplit(comps[i],split="_"))
  subType1<- twoSubTypes[1]
  subType2<- twoSubTypes[2]
#assign the aubtypes to the appropriate place in the cancer_subTypes matrix
cancer_subtypes[i,1]<-subType1
cancer_subtypes[i,2]<-subType2
  dir.create(paste("C:/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/5_class/5_class_pairwise/test/", subType1,"_", subType2, sep=""))
  setwd(paste("C:/Users/rds/Dropbox/PhD/RandomForest_Project/Fisher_Test_filtered_results/5_class/5_class_pairwise/test/", subType1,"_", subType2, sep=""))
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
  class1<-names(outcome[which(outcome==subType1)])
  class2<-names(outcome[which(outcome==subType2)])
  
  
  
  class1_2<-c(class1,class2)
  y<- factor(outcome[class1_2])#outcome class for random forest prediction
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
  #testSamples<- ySamples[which(TenFolds==1)]
  #trainSamples<- ySamples[which(TenFolds!=1)]
majLabs<-table(table(TenFolds))# the number of folds of a particular size, used to identify the majory fold size
majFoldSize<-names(majLabs[match(max(majLabs),majLabs)]) #which of the TenFolds class size is the majority out of the three classes?
majFoldNumbers<-which(table(TenFolds)==majFoldSize)# the class numbers which are of the majority class size

# this section is important for making sure that the training set is balanced for the outcome classes, by making sure the training set is composed of the two third folds that are of equal size. see table(TenFolds).

# majLabels<-which(table(TenFolds)==max(table(TenFolds)))
 trainSamples<-ySamples[which(TenFolds==majFoldNumbers[1]| TenFolds==majFoldNumbers[2])]
 testSamples<-ySamples[which(TenFolds!=majFoldNumbers[1]& TenFolds!=majFoldNumbers[2])]






  #test and train outcomes
  testOutcome<- y[which(TenFolds!=majFoldNumbers[1]& TenFolds!=majFoldNumbers[2])]
  trainOutcome<- y[which(TenFolds==majFoldNumbers[1]| TenFolds==majFoldNumbers[2])]
  
  ##############################################################################################################################
  
  #######################################
  
  #test and train data matrices
  testMutationMatrixLogical<-as.data.frame.matrix(t(mutationMatrixLogical[,testSamples]))
  trainMutationMatrixLogical<-as.data.frame.matrix(t(mutationMatrixLogical[,trainSamples]))
  
  # write number of protein coding genes carrying at least one mutation
  sink(paste("trainMutationLogical_proteins_caryying_one_mutation_",subType1,"_",subType2,".txt",sep=""))
  print(c("proteins carrying at least one protein across training set",length(which(colSums(trainMutationMatrixLogical)!=0))))
  sink()
  
  # write background kegg list
  write.table(rownames(mutationMatrixLogical), "mutationMatrixLogical_kegg_background.txt", row.names=FALSE,col.names=FALSE)
  
  setdiff(class1,intersect(names(trainOutcome),class1))
  
  #################################################################################################################################################
  ## filter out proteins according to the number of samples in which they carry mutations ################################
  
  tumourType<-unique(clinicalData2$tumourType)
  
  filterMuts<- function(mutationMatrix){
    #filter out the genes that do not carry mutations in more than 5% of the samples
    #5% of samples is 970/20
    
    FivePercent<- ceiling(dim(mutationMatrix)[1]/20)
    mutationMatrix_proteins<-names(which(colSums(mutationMatrix)>FivePercent))
    
    return(mutationMatrix_proteins)
  }
  
  
  train_ClinicalData2<- clinicalData2[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds),]
  #the names of the proteins carried through to the next stage, by cancer type
  filtered_protein_Names<-lapply(tumourType, function(x) filterMuts(trainMutationMatrixLogical[train_ClinicalData2$tumourType==x,]))
  names(filtered_protein_Names)<- tumourType
  
  # the unique list of proteins carrying mutations in at least 5% of at least cancer type
  proteins_analysis<- unique(unlist(filtered_protein_Names))
  
  # writing the list of proteins that made it through the filter
  write.table(proteins_analysis, paste("proteins_after_filter_",subType1,"_",subType2,".txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  ##redefining the train and testMutationMatrices.
  
  trainMutationMatrixLogical<- trainMutationMatrixLogical[,proteins_analysis]
  testMutationMatrixLogical<-testMutationMatrixLogical[,proteins_analysis]
  
  trainSamples<- rownames(trainMutationMatrixLogical)
  testSamples<- rownames(testMutationMatrixLogical)
  
  # may have to redefine the samples used if non of them have any proteins carrying mutations in the new matrix.
  ####################################################################################################################################################
  
  
  
  
  ###############################################################################################################################
  ##add extra test samples from the majority training class
  
  
  if(subType1==majorityClass){
    
    extraTestSamples<-setdiff(class1,intersect(names(trainOutcome),class1))## the majority samples left out of the training set
    extraTest<-t(mutationMatrixLogical[proteins_analysis,extraTestSamples])
    extraTest<-extraTest[-which(extraTestSamples%in%rownames(testMutationMatrixLogical)),]#remove majority samples that are already present in the test set
  }else{
    extraTestSamples<-setdiff(class2,intersect(names(trainOutcome),class2))## the majority samples left out of the training set
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
  #######################################################################################################################
  ######################################################################################################################
  ## filter out genes which carry too few mutations
  
  #colsumtrain<-table(colSums(trainMutationMatrixLogical))
  #hist(table(colSums(trainMutationMatrixLogical)), breaks=50)
  
  ##cumulative frequency
  ########plot(cumsum(colsumtrain)/dim(trainMutationMatrixLogical)[2])
  #top_quartile_proteins<-which(cumsum(colsumtrain)/dim(trainMutationMatrixLogical)[2]>0.75)[1]## tells me how many mutated samples there must be for a gene to be in the top 20% of genes according to number of samples carrying mutations.
  
  ##filter
  #trainMutationMatrixLogical<-trainMutationMatrixLogical[,names(which(colSums(trainMutationMatrixLogical)>=top_quartile_proteins))]
  #testMutationMatrixLogical<- testMutationMatrixLogical[,names(which(colSums(trainMutationMatrixLogical)>=top_quartile_proteins))]
  #####d<- density(colsumtrain, bw="SJ")
  #####plot(d)
  
  #####xfit<-seq(min(colsumtrain), max(colsumtrain), length=50)
  #####plot(dgeom(xfit, prob=0.5))
  
  #######################################################################################################################
  #######################################################################################################################
  
  
  
  fisherTestsGenes<-sapply(seq(1:dim(trainMutationMatrixLogical)[2]), function(x) fisher.test(table(trainOutcome, trainMutationMatrixLogical[,x]), alternative="two.sided")$p.value)
  significantTests<-which(fisherTestsGenes<0.05/dim(trainMutationMatrixLogical)[2])
  ##trainMutationMatrixLogical<-trainMutationMatrixLogical[,significantTests]# train dataset for significant F.tests oin training set
  ##testMutationMatrixLogical<- testMutationMatrixLogical[,significantTests]# test dataset for significant F tests
  
  
  ###output Fisher's test results to file
  ##significantTestsOutput<-data.frame(genes=colnames(trainMutationMatrixLogical),pvalues=as.numeric(fisherTestsGenes[significantTests]), stringsAsFactors=FALSE)
  
  ##significantTestsOutput<-significantTestsOutput[order(significantTestsOutput[,2], decreasing =TRUE),]
  ##write.csv(significantTestsOutput, file=(paste(subType1,"_",subType2,"_2_thirds_train_f.testsignificant.csv", sep="")), row.names=FALSE)
  ##bonferroni_pvalue<-paste(c("bonferroni p-value = 0.05/",length(fisherTestsGenes),",",0.05/length(fisherTestsGenes), sep=","))
  ##cat(bonferroni_pvalue, file=(paste(subType1,"_",subType2,"_2_thirds_train_f.testsignificant.csv", sep="")), append=TRUE)
  
  
  #########################################################################################################################
  #########################################################################################################################
  ## add SNV_TYPE measures to the train/testmutationMatrixLogical
  if(SNV_Types_And_clinical==TRUE){
    testMutationMatrixLogical<-as.data.frame(cbind(SNV_Types_freq_samples[match(testSamples,rownames(SNV_Types_freq_samples)),]), stringsAsFactors=TRUE)
    trainMutationMatrixLogical<-as.data.frame(cbind(SNV_Types_freq_samples[match(trainSamples,rownames(SNV_Types_freq_samples)),]), stringsAsFactors=TRUE)
  }else if (SNV_Types==TRUE){
    testMutationMatrixLogical<-as.data.frame(cbind(SNV_Types_freq_samples[match(testSamples,rownames(SNV_Types_freq_samples)),],testMutationMatrixLogical), stringsAsFactors=TRUE)
    trainMutationMatrixLogical<-as.data.frame(cbind(SNV_Types_freq_samples[match(trainSamples,rownames(SNV_Types_freq_samples)),],trainMutationMatrixLogical), stringsAsFactors=TRUE)
  }
  
 
  #####################################################################################################################
  ## add age and gender covariates to the train/testmutationMatrixLogical
  #test if the comparison inlcudes LAML or GBM. if so do not include stage in formation
  if(i %in% leaveStage){
  testMutationMatrixLogical<- as.data.frame(cbind(age=clinicalData2$age_at_initial_pathologic_diagnosis[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],testMutationMatrixLogical), stringsAsFactors=TRUE)
  testMutationMatrixLogical<- as.data.frame(cbind(gender=clinicalData2$gender[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],testMutationMatrixLogical), stringsAsFactors=TRUE)
  testMutationMatrixLogical<- as.data.frame(cbind(outcome=testOutcome,testMutationMatrixLogical))
 
  trainMutationMatrixLogical<- as.data.frame(cbind(age=clinicalData2$age_at_initial_pathologic_diagnosis[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))#  trainMutationMatrixLogical<- as.data.frame(cbind(gender=clinicalData2$gender[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))
  trainMutationMatrixLogical<- as.data.frame(cbind(gender=clinicalData2$gender[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))
  trainMutationMatrixLogical<- as.data.frame(cbind(outcome=trainOutcome,trainMutationMatrixLogical))
  
  }else{
    testMutationMatrixLogical<- as.data.frame(cbind(stage=clinicalData2$stage[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],testMutationMatrixLogical), stringsAsFactors=TRUE)
    testMutationMatrixLogical<- as.data.frame(cbind(age=clinicalData2$age_at_initial_pathologic_diagnosis[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],testMutationMatrixLogical), stringsAsFactors=TRUE)
    testMutationMatrixLogical<- as.data.frame(cbind(gender=clinicalData2$gender[match(intersect(testSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],testMutationMatrixLogical), stringsAsFactors=TRUE)
    testMutationMatrixLogical<- as.data.frame(cbind(outcome=testOutcome,testMutationMatrixLogical))
    
    trainMutationMatrixLogical<- as.data.frame(cbind(stage=clinicalData2$stage[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))
    trainMutationMatrixLogical<- as.data.frame(cbind(age=clinicalData2$age_at_initial_pathologic_diagnosis[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))
    trainMutationMatrixLogical<- as.data.frame(cbind(gender=clinicalData2$gender[match(intersect(trainSamples,clinicalData2$SampleIds),clinicalData2$SampleIds)],trainMutationMatrixLogical))
    trainMutationMatrixLogical<- as.data.frame(cbind(outcome=trainOutcome,trainMutationMatrixLogical))
  }
   
    testOutcome<-testMutationMatrixLogical$outcome
    

###############################age, gender and stage effectcs in the training set##########################
#gender table
genderTable<- table(trainOutcome,trainMutationMatrixLogical$gender)
XsqTest<-chisq.test(genderTable)


#age table
ageDataFrame <- as.data.frame(cbind(trainMutationMatrixLogical$age,trainOutcome))
ageAOV<- aov(trainMutationMatrixLogical$age~trainOutcome)
ageAOVANOVA<-anova(ageAOV)

pairwise.t.test(trainMutationMatrixLogical$age, trainOutcome, p.adjust="none")

if(i %in% leaveStage==FALSE){


#test for stage biases
Kwallis_stage<-kruskal.test(stage~ outcome, trainMutationMatrixLogical)
#stage bias is present
}
sink(file="age_stage_gender_training_set_bias_tests.txt",split=TRUE)
print("test for age bias in training set")
print(ageAOVANOVA)
print("test for gender bias in training set")
print(XsqTest)
if(i %in% leaveStage==FALSE){
print("test for stage bias in training set")
print(Kwallis_stage)
}
sink() 
#############################################################################################################################  
  
  
  
  
  #########################################################################################################################
  ### filter proteins by jaccard similarity
  #jaccardSim<- 1-vegdist(t(trainMutationMatrixLogical[,-c(1,2,3,4,5,6,7,8,9)]), method="jaccard", diag=FALSE, upper=TRUE)
  #jaccardSim1<- as.matrix(jaccardSim)
  #jaccardSim1Table<-table(jaccardSim1)
  #correlated<-findCorrelation(jaccardSim1, cutoff=0.90)
  #trainMutationMatrixLogical1<- trainMutationMatrixLogical[,-correlated]
  
  
  #mdsPoints<-cmdscale(jaccardSim1, k=2, eig=FALSE)
  #########################################################################################################################
  ############univariate filter random forest##########################log regression#####################################
  if(SNV_Types_And_clinical== FALSE){
    if(i %in% leaveStage){
    extraVars<-length(c("outcome","gender","age"))
    }else{
      extraVars<-length(c("outcome","gender","age","stage"))
    }
    #if SNV_Types are to be included in the model
    if(SNV_Types==TRUE){
      extraVars<- dim(SNV_Types_freq_samples)[2]+length(covs)
    }else{
      if(i %in% leaveStage){
        extraVars<-length(c("outcome","gender","age"))
      }else{
        extraVars<-length(c("outcome","gender","age","stage"))  
    }
    }
  }
    genes <- colnames(trainMutationMatrixLogical)[-seq(1:extraVars)]
    
   
    
    ## if I only want to use clinical Data this statement changes the train and test data to just use age and gender information
    if(clinicalOnly==TRUE){
      trainMutationMatrixLogical2<- trainMutationMatrixLogical[,1:extraVars]
      testMutationMatrixLogical2<- testMutationMatrixLogical[,1:extraVars]
    
  }else if(SNV_Types_And_clinical==TRUE){
    trainMutationMatrixLogical2<-trainMutationMatrixLogical[,-1]
    testMutationMatrixLogical2<-trainMutationMatrixLogical[,-1]
  }else if(SNV_Types==FALSE & clinicalOnly==FALSE){
    trainMutationMatrixLogical2<- trainMutationMatrixLogical[,-1]
    testMutationMatrixLogical2<- testMutationMatrixLogical[,-1]
  }
  
  

###################################################################################################

##################################################################################################



  #output training and test set sizes
  sink(paste("training_and_test_set_sizes",subType1,"_",subType2,".txt",sep=""), append=FALSE)
  print(c("training set size ",length(trainSamples)))
  print(c("training set outcome table", table(trainOutcome)))
  print(c("test set size ",length(testSamples)))
  sink()
  
#############################################################################################################################  
#model building
if(clinicalOnly==TRUE){
  # need to get ROC objects and AUCs for plotting and AUC for models just using the clinical variables
  
 # subtype1 must be the positive logistic regression class. it must be the second outcome factor
  
  trainOutcome<- factor(trainOutcome, levels=c(subType2, subType1))
  testOutcome<- factor(testOutcome, levels=c(subType2,subType1))
  

  trueClass<-factor(testOutcome, levels=c(subType2,subType1))
  
  
  #recode grade outcome levels to be low high
  trainMutationMatrixLogical2$outcome<-factor(trainMutationMatrixLogical2$outcome, levels=c(subType2, subType1))# squamous is the positive class
  testMutationMatrixLogical2$outcome<-factor(testMutationMatrixLogical2$outcome, levels=c(subType2, subType1))# squamous is the positive class
  
  
  
  train_GLM<- glm(outcome~., data=trainMutationMatrixLogical2, family="binomial")
  train_GLM_summary<- summary(train_GLM)
  train_GLM_OR<-exp(cbind(OR = coef(train_GLM), confint(train_GLM)))
  
  sink(paste("train_GLM_stats_",subType1,"_",subType2,".txt"))
  print(train_GLM_summary)
  print(train_GLM_OR)
  sink()
  
  
  #classification models list
  regOR<- round(train_GLM_OR[,1],3)
  regCI<-cbind(paste("(",round(train_GLM_OR[,2],3),", ",round(train_GLM_OR[,3],3),")", sep=""))
  regORCI<- cbind(regOR,regCI)
  regPval1<- format(round(train_GLM_summary[[12]][,4],4), scientific=TRUE)# rounded scientific notation
  regPval2<- train_GLM_summary[[12]][,4]
  
  RegModels[[i]]<- data.frame(regORCI,regPval1,regPval2, stringsAsFactors=FALSE)
  colnames(RegModels[[i]])<-c("OR","ORCI","Pvalround","rawPval")
  # output the RegModels object to file at the end of the anlaysis
  
  
  #modelfit
  
  dev<-with(train_GLM, null.deviance-deviance)# difference in deviance between the null model and our model
  dfs<-with(train_GLM, df.null - df.residual)# df for the difference between the two models = difference between the number of predictors
  
  train_GLM_modelfit<-with(train_GLM, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))# pvalue for a chisquare comparing difference between the deviance statistic for the null modela nd our model.
  loglik<-logLik(train_GLM)
  
  modelFit<- list(deviance_diff=dev,degrees_of_freedom=dfs,chisquare=train_GLM_modelfit, loglikelihood_of_model=loglik)
  sink(paste("train_GLM_model_",subType1,"_",subType2,".txt"), append =TRUE)
  print(train_GLM_summary)
  print(modelFit)
  sink()
  
  test_GLM<- predict(train_GLM, newdata=testMutationMatrixLogical2, type="response")
  
  
  #####################################################################################################
  #####sensitivity and specificity of the GLM model
  
  true_pred_table<-table(trueClass, test_GLM>0.5)# contingency table for true class and predicted class
  true_pred_table2<-true_pred_table# The rows are the true class and the columns are the predicted
  true_pred_table2[1,1]<-true_pred_table[2,2]# in order for the confusion matrix function to work, the TP and FP indeces needed to be swapped
  true_pred_table2[2,2]<- true_pred_table[1,1]
  colnames(true_pred_table2)<-c(subType1,subType2)
  rownames(true_pred_table2)<-c(subType1,subType2)
  
  
  colnames(true_pred_table)<- rownames(true_pred_table)
  
  stats<-confusionMatrix(true_pred_table2)
  
  sink(paste("Test_GLM_classifier_stats",subType1,"_",subType2,".txt", sep=""), append=TRUE)
  
  print("High class is defined as a positive and Low class is a negative. index [1,1] is TP and index[2,2] is TN")
  print(confusionMatrix(true_pred_table2))
  sink()
  
  classStats[[i]]<-c(stats$byClass,positive=stats$positive,stats$overall)
    
}else{
#############################################################################################################################
  # random forest Model Building

#rfe control fuction
set.seed(10)

# subtype1 must be the positive logistic regression class. it must be the FIRST outcome factor in Random Forest

trainOutcome<- factor(trainOutcome, levels=c(subType1, subType2))
testOutcome<- factor(testOutcome, levels=c(subType1,subType2))


#sizetol<-c(rep(5,10))
sizetol<-c(rep(2,10))
# set the number of changepoints
#chpts<- rep(1,nComps)
#chpts[c(2,6,8)]<-2# for the comparisons with the slightly different distributions of accuracy statistic- no gender
#chpts[c(4,6,8)]<-2# for the comparisons with the slightly different distributions of accuracy statistic-all


if(sizeTolerant1pcnt){
#define the picksize tolerance
rfFuncs2<-rfFuncs
#rfFuncs2$selectSize<-pickSizeTolerance


rfFuncs2$selectSize<-function (x, metric="Accuracy", tol=sizetol[i], maximize) 
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





#selectModelsize<- function (x, metric="Accuracy", maximize, chpt)
#{
#  #picksizeResults<- data.frame(Accuracy=x$results$Accuracy, Variables=seq(1:dim(x$results)[1]))
#  picksizeResults<- x[,metric]
#  changepnt<- cpt.meanvar(picksizeResults,penalty="AIC",pen.value="n", method="AMOC",Q=chpt)
#  #changepnt<- cpt.mean(picksizeResults,penalty="AIC",pen.value="n",test.stat="CUSUM", method="BinSeg",Q=2)
#  modelsize<-changepnt@cpts[1]
#  as.integer(modelsize)
#  return(modelsize)
#  
#}
#environment(selectModelsize)<- as.environment("package:caret")# add my function to the caret environment so it can pass aruments to other functions within caret
#rfFuncs2$selectSize<-selectModelsize# assign the select model size criterea in caret to be th selectModelsize function


rfFuncs2$selectVar<-pickVars

subsets<-c(seq(1:100))

ctrl <- rfeControl(functions = rfFuncs2, method = "repeatedcv", repeats = 5, verbose = FALSE, returnResamp = "final", allowParallel=TRUE, rerank=TRUE)
trainModel<-rfe(trainMutationMatrixLogical2, trainOutcome, subsets, rfeControl=ctrl,proximity=TRUE, localImp=TRUE)

}else{
  subsets<-c(seq(1:100))
  
  ctrl <- rfeControl(functions = rfFuncs, method = "repeatedcv", repeats = 5, verbose = FALSE, returnResamp = "final", allowParallel=TRUE, rerank=TRUE)
  trainModel<-rfe(trainMutationMatrixLogical2, trainOutcome, subsets, rfeControl=ctrl,proximity=TRUE, localImp=TRUE)

}



#plot how each of the size tolerant models was arrived at
picksizeResults<- data.frame(Accuracy=trainModel$results$Accuracy, Variables=seq(1:dim(trainModel$results)[1]))


#picking the sizetolerance
changepnt<- cpt.meanvar(picksizeResults[,1],penalty="AIC",pen.value="n",method="AMOC",Q=chpts[i])
#changepnt<- cpt.mean(picksizeResults[,1],penalty="AIC",pen.value="n", method="AMOC",Q=1)


pdf(paste("changepoint_analysis",subType1,"_",subType2,".pdf"), 12,10, pointsize=12)
plot( changepnt, ylab="tolerance", xlab="subset size")
dev.off()
subsetmodelsize[[i]]<- changepnt
subsetmodelsizeselect[[i]]<-changepnt@cpts# forthe list after all analyses are run

optimalModel<- as.integer(changepnt@cpts[1])


# write the acuracies and tolerances to the results folder
write.table(picksizeResults, file=paste("tolerance_results_",subType1,"_", subType2,".txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

#in order to calculate tolerance I need the maximum accuracy from the list of models
maxAcc <- max(picksizeResults$Accuracy)
picksizeResults$Tolerance <- 100-(abs((picksizeResults$Accuracy / maxAcc))*100)


within1pcnt<- pickSizeTolerance(picksizeResults, metric="Accuracy", tol=sizetol[i], maximize=TRUE)
smallest<- pickSizeBest(picksizeResults, metric="Accuracy", maximize=TRUE)
#png(paste("size_tolerance_plots_",subType1,"_",subType2,".png"), 1200,1000, antialias="cleartype")
pdf(paste("size_tolerance_plots_",subType1,"_",subType2,".pdf"), 12,10, pointsize=12)
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

pdf(paste("accuracy_hist_distribution",subType1,"_",subType2,".pdf",sep=""),12,10, pointsize=12)
hist(picksizeResults[,1],xlab="Accuracy", main="", border=NA, col="lightgray", cex=2)
#curve(dnorm(x, mean=mean(picksizeResults[,1]), sd=sd(picksizeResults[,1])), add=TRUE, col=col.plot[2], lwd=2)
lines(density(picksizeResults[,1],adjust=2), col=col.plot[2],lwd=2)
dev.off()



          


     
# variable importance parameters

trainModelImportances<- importance(trainModel$fit, type="1")
trainModelImportances<- trainModelImportances[order(trainModelImportances, decreasing=TRUE),]

  #print the model parameters     
  sink(paste("Train_",subType1,"_",subType2,"_2_thirds_randomForest_results.txt", sep=""))
  print(trainModel$optsize)
  print("trainmodel classes, the first factor is the positive class")
  print(trainModel$fit$classes)
  print(trainModel)
  print(trainModel$optVariables)
  print(cbind(trainModelImportances))
  sink()
  
  # optimal model variables for kegg pathway enrichment
  write.table(trainModel$optVariables,paste("kegg_protein_list",subType1,"_",subType2,".txt", sep=""), row.names=FALSE, col.names=FALSE )
  
######plots and model selection


  
  testClasses <- predict(trainModel, testMutationMatrixLogical2)# the test run classes
  #  confusionMatrix(testClasses$pred,testOutcome)# the test run accuracy results
  sink(paste("ROC_curve_data_", subType1,"_",subType2,".txt", sep=""))
  print(testClasses)
  sink()
  
  
  sink(paste("Test_balanced_",subType1,"_", subType2,"_2_thirds_randomForest_results.txt", sep=""), append=TRUE)
  print(testClasses)
  print(confusionMatrix(testClasses$pred,testOutcome))
  sink()

stats<-confusionMatrix(testClasses$pred,testOutcome)
classStats[[i]]<-c(stats$byClass,positive=stats$positive,stats$overall)

#################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
# clustered heatmaps to represent what is happening over the forest
if(trainModel$optsize!=dim(trainMutationMatrixLogical2)[2]& trainModel$optsize!=2){
library(gplots)
library(colorspace)


#vars4Heatmaps<-c(colnames(trainMutationMatrixLogical2)[seq(1,(extraVars-1))], trainModel$optVariables)# the list of variables to include in my heatmap. clinicl variables should update according to the number of extraVars.

vars4Heatmaps<-trainModel$optVariables
trainMutationMatrixLogical2_model<- trainMutationMatrixLogical2[,vars4Heatmaps]
trainMutationMatrixLogical2_model_Logical<- trainMutationMatrixLogical2_model

#####################################################################################################
#fisher's exact test for overlap of mutations
#all pairwise comparisons between cancer classes

#indeces of clinicalVariables do I can remove them from the possible fisher's tests I want to do
stageIndex<-match("stage",colnames(trainMutationMatrixLogical2_model))
ageIndex<-match("age",colnames(trainMutationMatrixLogical2_model))
genderIndex<-match("gender", colnames(trainMutationMatrixLogical2_model))
clinicalIndeces<- c(stageIndex,ageIndex,genderIndex)
clinicalIndeces<- clinicalIndeces[!is.na(clinicalIndeces)]# final index to use once missing clinical data indeces are removed

if(length(clinicalIndeces)>0){
# gene names for comparisons
pairwiseGeneCompsNames<-combn(colnames(trainMutationMatrixLogical2_model_Logical)[-clinicalIndeces],2)

indexList<- seq(1:trainModel$optsize)[-clinicalIndeces]
}else{
  pairwiseGeneCompsNames<-combn(colnames(trainMutationMatrixLogical2_model_Logical),2)
  indexList<- seq(1:trainModel$optsize)
  
}

pairwiseIndexComps<- combn(indexList,2)# use this one for defining the fisher's tests



       #carry out fisher tests for association between model protein variables
nftests<- dim(pairwiseIndexComps)[2]
if(nftests>1){
       modelVarsFTests<-lapply(seq(1:dim(pairwiseIndexComps)[2]), function(x) fisher.test(table(trainMutationMatrixLogical2_model_Logical[,pairwiseIndexComps[1,x]],trainMutationMatrixLogical2_model_Logical[,pairwiseIndexComps[2,x]],dnn=(c(pairwiseGeneCompsNames[1,x], pairwiseGeneCompsNames[2,x])))))
       #get p values
       modelVarsFTestsPvals<- sapply(seq(1:dim(pairwiseIndexComps)[2]), function(x) modelVarsFTests[[x]]$p.value)
       #adjust pvalues BH
       modelVarsFTestsPvalsAdj<-p.adjust(modelVarsFTestsPvals, method="BH")
       #get significant values at FDR<0.1
       modelVarsFTestsPvalsSig<-which(modelVarsFTestsPvalsAdj<0.1)
       # get the comparisons that are significant
       modelVarsAssocComps<- sapply(seq(1:length(modelVarsFTestsPvalsSig)), function(x) pairwiseGeneCompsNames[,x])
       modelVarsAssocCompsPvals<- lapply(modelVarsFTestsPvalsSig,function(x) modelVarsFTestsPvals[[x]])
       modelVarsAssocCompsFDR<- lapply(modelVarsFTestsPvalsSig,function(x) modelVarsFTestsPvalsAdj[[x]])
       
       modelVarsAssocCompsoutput<-cbind(protein1=modelVarsAssocComps[1,], protein2=modelVarsAssocComps[2,],pvalue=modelVarsAssocCompsPvals,fdr=modelVarsAssocCompsFDR)
       write.table(modelVarsAssocCompsoutput, file=paste("associatedModelVariables_proteins",subType1,"_",subType2,".txt"), sep="\t", col.names=TRUE, row.names=FALSE)
}else{
  
  modelVarsFTests<-fisher.test(table(trainMutationMatrixLogical2_model_Logical[,pairwiseIndexComps[1,1]],trainMutationMatrixLogical2_model_Logical[,pairwiseIndexComps[2,1]],dnn=(c(pairwiseGeneCompsNames[1,1], pairwiseGeneCompsNames[2,1]))))
  sink(file=paste("associatedModelVariables_proteins",subType1,"_",subType2,".txt"), append=FALSE)
sink()    
}
       
       #####################################################################################################
       
       
       
       
       #######################################use the random forest variables ###################################################
       
       ##heatmap using the vegdist functions to create dendrograms and use those dendrograms in the reslting heatmap.
       #clusterRows
       
       proximity<- trainModel$fit$proximity
       rd<-1-proximity# dissimilarity 1-proximity
       rd<- as.dist(rd)
       rc<-hclust(rd, method="ward.D")
       
       
       
       v.importance<-t(trainModel$fit$localImportance)
       cd<-1-abs(cor(v.importance, method="pearson"))# correlation distance 1- absolute correlation coefficient. dissimilarity.
       cd<-as.dist(cd)
       cc<- hclust(cd, method="ward.D")
       
       
       #recode FEMALE to 0 and Male to 1# need to modify this so that male and female can be visualised, try 0.5 and 0.7
       if("gender" %in% colnames(trainMutationMatrixLogical2_model)){
         trainMutationMatrixLogical2_model$gender<-as.character(trainMutationMatrixLogical2_model$gender)
         trainMutationMatrixLogical2_model$gender[which(trainMutationMatrixLogical2_model$gender=="FEMALE")]<-0.5
         trainMutationMatrixLogical2_model$gender[which(trainMutationMatrixLogical2_model$gender=="MALE")]<-0.6
         trainMutationMatrixLogical2_model$gender<- as.numeric(trainMutationMatrixLogical2_model$gender)
       }
       ###convert trainMutationMatrixLogical2_Model to a numeric matrix
       
       
       trainMutationMatrixLogical2_model[trainMutationMatrixLogical2_model==TRUE]<-0.4
       
       trainMutationMatrixLogical2_model<-as.matrix.data.frame(trainMutationMatrixLogical2_model)
       
       
       ##heatmap using the vegdist functions to create dendrograms and use those dendrograms in the resulting heatmap. by clustering raw binary mutation matrix
       #clusterRows
       #rd<-dist(trainMutationMatrixLogical2_model, method="binary")
       #rc<-hclust(rd, method="ward.D")
       
       #clustercolumns
       #cd <- vegdist(t(trainMutationMatrixLogical2_model), method="jaccard")
       #cc<- hclust(cd, method="ward.D")
       
       
       
       # define the palette for the 5 different cancer types
       mycol<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99")
       
       # define the vector of color names of each sample in the training set
       classCols<- rep(0,length(trainOutcome))
       subTypeIndex<-c(match(subType1,classes$classes),match(subType2,classes$classes))
       classCols[which(trainOutcome==names(table(trainOutcome))[1])]<- mycol[ subTypeIndex[1]]
       classCols[which(trainOutcome==names(table(trainOutcome))[2])]<- mycol[ subTypeIndex[2]]
       
       
       colage<-heat_hcl(length(table(clinicalData2$age_at_initial_pathologic_diagnosis)), h = 120, c. = c(0,80), l = c(90, 30), power = 1.3)
       colstage<-heat_hcl(4,h=60, c.=c(0,80), l=c(90,30), power=1.3)
       #colMut<-rainbow_hcl(2, start = 100, end = 210)
       colMut<-c("#FFFFFF","#BEBEBE")
       
       colgender<-c("#e31a1c","#fdbf6f")
       
       my_palette<-c(colMut,colgender,colstage,colage)
       par(xpd=FALSE)
       colorbreaks<-(c(0,0.4,0.5,0.6,1,2,3,4,sort(unique(clinicalData2$age_at_initial_pathologic_diagnosis), decreasing=FALSE),91)-0.01)#breaks based on the number of colors in my_palette
       png(paste(subType1,"_",subType2,"heatmap.png",sep=""),1200,1200, antialias="cleartype")
       #png(paste(subType1,"_",subType2,"heatmap.png",sep=""),11*300,11*300, res=300, pointsize=5)

       par(mar=c(10,10,10,10))
       par(pty="m")
       test<-heatmap.2(trainMutationMatrixLogical2_model, Rowv=as.dendrogram(rc), Colv=as.dendrogram(cc),col=my_palette,breaks=colorbreaks, RowSideColors=classCols, trace="none",labRow=FALSE, density.info="none", margins=c(12,6), key=FALSE, cexCol=2.4)
       par(xpd=TRUE)
       legend(0.000001,1.1, legend=(c(levels(trainOutcome),"not mutated", "mutated")), fill=c(mycol[subTypeIndex],my_palette[c(1,2)]), pch, bty="n", cex=2.4)
       #legend(0.000001,1, legend=(c(levels(trainOutcome),"not mutated", "mutated")), fill=c(mycol[subTypeIndex],my_palette[c(1,2)]), pch, bty="n", cex=2.4)

dev.off()
       
       ##########################################################################################################################################################
       ###########################################################################################################################################################
}
}

  #################################################################################################################################
  ## ROC curves
    
  library(pROC)

if(clinicalOnly==FALSE){
  
  #true class. for ROC calculation the second factor is the successful classification
  trueClass<-factor(testOutcome, levels=c(subType2,subType1))
  #trueClass<-reorder(trueClass,new.order=c(2,1))
}
  #trueClass<-revalue(trueClass,c("squamous"="other", "BLCA"="other","leukemia"="other", "adenocarcinoma"="other"))# re-assigning the factor labels to be "adenocarinoma" or "other"
  # now trueClass is the true gold standard class. This should be used for building ROC curves
  
  # changing the clas used to calculate AUC
  #if(SNV_Types_And_clinical){
  #  
  #  clssvec<-testClasses[,subType1]
  #} else if(i==4| i==7){
  #  
  #  clssvec<-testClasses[,subType2]
  #}else{clssvec<- testClasses[,subType1]
  #}

if(clinicalOnly==TRUE){
  trueClass<-factor(testOutcome, levels=c(subType2,subType1))
  clssvec<-test_GLM# probability of being classified as subType1.
}else{
  clssvec<-testClasses[,subType1]# probability of bbeing classified as the positive class in random forest
  

}
  test<-roc(trueClass,clssvec)
  AUCs[i]<-test$auc# assign AUC to AUC object for ROC curve creation later
  ROC_objects[[i]]<-test# assign ROC object to the list of ROC_objects for ROC curve creation later

  
  sink(paste("Test_balanced_",subType1,"_", subType2,"_2_thirds_randomForest_results.txt", sep=""), append=TRUE)# append ROC AUC results
  print("ROC RESULTS")
  print(test)
  sink()


}
  

# output classification statistics
classStatsOutput<-as.data.frame.list(classStats, stringsAsFactors=FALSE)
colnames(classStatsOutput)<- comps
classStatsOutput<-as.data.frame(t(classStatsOutput), stringsAsFactors=FALSE)
positiveclass<- classStatsOutput[,8]#positive classes
classStatsOutput<-classStatsOutput[,-8]#all the other classification statistics

classStatsOutput<-as.matrix(sapply(classStatsOutput, as.numeric))# convert to numeric
classStatsOutput1<-as.data.frame(cbind(classStatsOutput,positive=positiveclass))
classStatsOutput1<-cbind(classStatsOutput,positive=positiveclass)# raw scores
rownames(classStatsOutput1)<-comps

classStatsOutputRound<-as.data.frame(cbind(round(classStatsOutput,3),positive=positiveclass))# rounded scores
rownames(classStatsOutputRound)<- comps

# write the raw classification stats
write.table(classStatsOutput1, file="classification_stats_raw.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

# write the rounded classification stats
write.table(classStatsOutputRound, file="classification_stats_rounded.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)


#######################
# output the changpoint analysis selected model sizes
if(clinicalOnly==FALSE){
if(sizeTolerant1pcnt==TRUE){
 sink("model_selection_changepoints.txt", append=FALSE)
print(subsetmodelsizeselect)
sink()

# output the model selection data

sink("changepoint_analysis.txt",append=FALSE)
print(subsetmodelsize)
print(subsetmodelsizeselect)
sink()


}
}else{
  names(RegModels)<-comps
  sink("RegressionModels.txt", append=FALSE)
  print(RegModels)
  sink()
  }

# final ROC-plot
col.plot<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")
  palette(col.plot)
  #col.rainbow<- col.rainbow<- rainbow_hcl(10, c=70, l=50, start = 30, end= 300)
  #palette(col.rainbow)

for(i in 1:nComps){
roctest<-ROC_objects[[i]]# the ROC object for the right test

if(i==1){
    
    png("allcomparisons_ROC.png", width=1200, height=1200)
    par(mar=c(5,5,5,25))
    par(pty="s")
    par(xpd=FALSE)
    plot(roctest, xlim=c(1,0), ylim=c(0,1), col=alpha(i,alpha=0.8),lwd=8, ann=FALSE, cex.axis=2, bg="transparent",bty="L" )
  }else{
    
    lines(roctest, col=alpha(i,0.8), lwd=8)
  }
  
}
  

AUCs<-round(AUCs,2)
comps2<-gsub("_"," ",comps)
textForLegend<-cbind(comps2,AUCs)
legendText<-sapply(1:dim(textForLegend)[1], function(x) paste(textForLegend[x,1],"=",textForLegend[x,2]))

legend(0.55,0.54,legendText,title="AUC",col=alpha(c(1:10), alpha=0.8),lwd=8, cex=2.4, lty=1, bty="n")
dev.off()




#old legend values 0.47,0.47


#############################################################################################################################################
#############################################################################################################################################
# clustered heatmaps to represent what is happening over the forest

library(gplots)
library(colorspace)


#vars4Heatmaps<-c(colnames(trainMutationMatrixLogical2)[seq(1,(extraVars-1))], trainModel$optVariables)# the list of variables to include in my heatmap. clinicl variables should update according to the number of extraVars.

vars4Heatmaps<-trainModel$optVariables
trainMutationMatrixLogical2_model<- trainMutationMatrixLogical2[,vars4Heatmaps]
trainMutationMatrixLogical2_model_Logical<- trainMutationMatrixLogical2_model

#####################################################################################################
#fisher's exact test for overlap of mutations
#all pairwise comparisons between cancer classes

#indeces of clinicalVariables do I can remove them from the possible fisher's tests I want to do
stageIndex<-match("stage",colnames(trainMutationMatrixLogical2_model))
ageIndex<-match("age",colnames(trainMutationMatrixLogical2_model))
genderIndex<-match("gender", colnames(trainMutationMatrixLogical2_model))
clinicalIndeces<- c(stageIndex,ageIndex,genderIndex)
clinicalIndeces<- clinicalIndeces[!is.na(clinicalIndeces)]# final index to use once missing clinical data indeces are removed

# gene names for comparisons
pairwiseGeneCompsNames<-combn(colnames(trainMutationMatrixLogical2_model_Logical)[-clinicalIndeces],2)

indexList<- seq(1:trainModel$optsize)[-clinicalIndeces]
pairwiseIndexComps<- combn(indexList,2)# use this one for defining the fisher's tests



sapply(fisher.test(table(trainMutationMatrixLogical2_model_Logical[,pairwiseIndexComps[1,1]],trainMutationMatrixLogical2_model_Logical[,pairwiseIndexComps[2,1]],dnn=(c(pairwiseGeneCompsNames[1,1], pairwiseGeneCompsNames[2,1]))))
#carry out fisher tests for association between model protin variables
#modelVarsFTests<-lapply(seq(1:dim(pairwiseIndexComps)[2]), function(x) fisher.test(table(trainMutationMatrixLogical2_model_Logical[,pairwiseIndexComps[1,x]],trainMutationMatrixLogical2_model_Logical[,pairwiseIndexComps[2,x]],dnn=(c(pairwiseGeneCompsNames[1,x], pairwiseGeneCompsNames[2,x])))))
#get p values
#modelVarsFTestsPvals<-sapply(seq(1:dim(pairwiseIndexComps)[2]), function(x) modelVarsFTests[[x]]$p.value)
#adjust pvalues BH
#modelVarsFTestsPvalsAdj<-p.adjust(modelVarsFTestsPvals, method="BH")
#get significant values at FDR<0.1
#modelVarsFTestsPvalsSig<-which(modelVarsFTestsPvalsAdj<0.1)
# get the comparisons that are significant
#modelVarsAssocComps<- sapply(seq(1:length(modelVarsFTestsPvalsSig)), function(x) pairwiseGeneCompsNames[,x])
#modelVarsAssocCompsPvals<- lapply(modelVarsFTestsPvalsSig,function(x) modelVarsFTestsPvals[[x]])
#modelVarsAssocCompsFDR<- lapply(modelVarsFTestsPvalsSig,function(x) modelVarsFTestsPvalsAdj[[x]])

#modelVarsAssocCompsoutput<-cbind(protein1=modelVarsAssocComps[1,], protein2=modelVarsAssocComps[2,],pvalue=modelVarsAssocCompsPvals,fdr=modelVarsAssocCompsFDR)
#write.table(modelVarsAssocCompsoutput, file=paste("associatedModelVariables_proteins",subType1,"_",subType2,".txt"), sep="\t", col.names=TRUE, row.names=FALSE)


#####################################################################################################
q()