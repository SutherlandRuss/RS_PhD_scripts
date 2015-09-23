

###########################################################################################
## FUNCTION TO LOAD CLINICAL DATA

clinicalDataTable<- function(rectumClinical.path, tumourLabel){
  #the names of the rectumClinical files
  rectumClinical.files<- list(list.files(rectumClinical.path))[[1]][1]
  rectumClinicalTables<- lapply(seq(1:length(rectumClinical.files[[1]])), function(x) cbind(read.table(file.path(rectumClinical.path,rectumClinical.files[[1]][x]) , header=TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("NA", "[Not Available]", "<NA>", ""))))
  names(rectumClinicalTables)<- unlist(rectumClinical.files)
  
  rectumClinicalTablesSamples<- lapply(rectumClinicalTables, function(x) strsplit(x[,1], split="-", fixed =TRUE))
  #variables for the table
  TablesForTable<- which(sapply(rectumClinicalTablesSamples, function(x)length(x[[1]]))==3)
  #modified table list based on the data in tables with duplicate SampleIDs. ommitted tables are not needed, drug table may be included at a later date.
  
  if(tumourLabel!="breast"){
    TablesForTable<- TablesForTable[-c(3,6,8)]# This needs to be changed manually
  }else{
    TablesForTable<-TablesForTable[-c(2,5,7)]# for breast cancer
  }
  #Identify the tables which use individual IDs as opposed to tissue sample ids which can come from the same patient.
  rectumClinicalTablesSamples<-lapply(seq(1:length(TablesForTable)),function(x) rectumClinicalTablesSamples[[TablesForTable[x]]])
  #get the appropriate Clinical Tables
  #no factor types here yet
  rectumClinicalTables<-lapply(seq(1:length(TablesForTable)),function(x)rectumClinicalTables[[TablesForTable[x]]])
  #get the sample names for each row of each metadata table
  rectumClinicalTablesSampleIDs1<- lapply(rectumClinicalTablesSamples, function(x) unlist(lapply(x, function(y) paste(y[1:3],collapse="")), use.names=FALSE))
  # assign sample ids to rectumClinicalTables. Now each table should have the same sampleIDS and I can look for uniques etc...
  #no factor tyes here yet
  rectumClinicalTables<-lapply(seq(1:length(rectumClinicalTables)), function(x) data.frame(rectumClinicalTables[[x]],rectumClinicalTablesSampleIDs1[[x]], stringsAsFactors=FALSE))
  
  #rename the newly added sample IDs column
  colnames(rectumClinicalTables[[1]])[length(rectumClinicalTables[[1]][1,])]<-"SampleIds"
  #assign new sampleid column names to all tables of rectumClinical Samples
  for( i in seq_along(rectumClinicalTables)){
    colnames(rectumClinicalTables[[i]])[length(rectumClinicalTables[[i]][1,])]<- "SampleIds"
  }
  #get the list of unique samples across all files
  uniqueSamples<-as.list(as.character(unique(unlist(lapply(rectumClinicalTables, function(x) x$SampleIds), use.names=FALSE))))
  
  
  # find the number of variables across all of my tables
  vnumber<- function(tabList){
    vnum<-0
    for(i in seq_along(tabList)){
      vnum<-vnum+((dim(tabList[[i]])[2])-2)# the -2 is because there are two id columns in each table that I don't want
    }
    return(vnum)
  }
  
  vn<-vnumber(rectumClinicalTables)
  # find the variable names
  #vars<-lapply(rectumClinicalTables, function(x) colnames(x)[c(-1,-(length(x[1,])))])
  #varNames<-unlist(vars)
  
  # create a table of all data from all tables based on unique sampleIDs
  
  #The join_all function acts like a mysql join, to combine all of the listed tables based on matching "sampleIds" variables.
  clinicalData<- join_all(rectumClinicalTables, by="SampleIds", type="full", match="first")
  #disease type matrix
  tumourType<- matrix(tumourLabel, ncol=1, nrow=length(clinicalData[,1]))
  # Add disease identifier to the clinicalData
  #no factors at this point
  clinicalData<-data.frame(clinicalData,tumourType, stringsAsFactors=FALSE)
  #change the SampleIds names to sampleID
  colnames(clinicalData)[1]<-"sampleID"
  return(clinicalData)
}


###########################################################################################################################



###########################################################################################################################
# function to combine metadata tables and ensures that variables remain in the correct order so that columns match netween the two joined tables
combineMetadata<- function(table1,table2){
  commonVariables<-intersect(colnames(table1), colnames(table2))
  table1<- table1[,commonVariables]
  table2<- table2[,commonVariables]
  table1<-table1[,order(match(colnames(table1), colnames(table2)))]
  combinedT<-as.data.frame(rbind(table1,table2), stringsAsFactors=FALSE)
  rownames(combinedT)<-combinedT[,1]
  return(combinedT)
}

###########################################################################################################################
#####################################################################################################################################################################
##loading seqdata
#####################################################################################################################################################################
## FUNCTION TO LOAD SEQUENCE DATA

createScoreTable<- function(dataFile){
  scores<- read.table(dataFiles , header=TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  scores<- unique(scores)
  #Splitting the TCGA barcode in to the minimum  number of fields to uniquely identify individuals
  samples<- strsplit(scores$Tumor_Sample_Barcode, split = "-", fixed = TRUE)
  # Unique sample IDs
  #sampleID<- matrix(sapply(1:length(samples),function(x) paste0(samples[[x]][1], samples[[x]][2], samples[[x]][3])))
  sampleID<- unlist(lapply(samples,function(x) paste(x[1:3],collapse="")), use.names=FALSE)
  
  #adding the sampleID to the geneScore dataframe
  scores <-cbind(scores, sampleID)
  scores$sampleID<- as.character(scores$sampleID)
  
  return(scores)
}

#####################################################################################################################################################################
##

##################################################################################################################
# function to get sequence center label
seqCentre<- function(scoresVariable){
  longID<- lapply(scoresVariable, function(x) strsplit(x, split="-", fixed =TRUE))
  seq_center_ID<-lapply(longID, function(x) unlist(lapply(x, function(y) y[7]), use.names=FALSE))
  seq_center_ID[seq_center_ID=="09"]<-"WashU"
  seq_center_ID[seq_center_ID=="10"]<-"BCM"
  seq_center_ID<-unlist(seq_center_ID)
  return(seq_center_ID)
}
####################################################################################################################


#######################################################################################################################
#function to label specific cancer samples accordign to their particular type of cancer
labelCancer<- function(scores,files){
  
  cancers<- sapply(seq(1:length(files)), function(x)unlist(strsplit(files[[x]],"[.]"))[1])
  cancer_type <- lapply(seq(1:length(files)), function(x) rep(cancers[x], length(scores[[x]][,1])))
  lScores<- lapply(seq(1:length(files)), function(x) cbind.data.frame(cancer_type=cancer_type[[x]],scores[[x]]))# when i leave out the cancer label this table is identical to scores. can I put it in to a for loop instead of a series of lapply and sapply?
  
  return(lScores)
}
#########################################################################################################################
#########################################################################################################################
# Function to find the intersection between two samples lists and return a dataset with only the intersecting samples

extractIntersectSamples<- function(dataNames,intersectList,data){
  if(class(dataNames)!= "character"||class(intersectList)!="character"){
    stop("Both SeqInfoNames and MetadataNames must be of character type")}else{
      SampleMatchI<- which(dataNames%in%intersectList)
      NewSeqInfo <- data[SampleMatchI,]
      #NewSeqInfo<- NewSeqInfo[order(NewSeqInfo$sampleID),]# it is better if I order the sample names after uniqueing
      return(NewSeqInfo)
    }
}
#########################################################################################################################
#########################################################################################################################
# Function to get non-silent mutations from the initial "scores" input from the somatic.maf files

getFuncMutations<-function(scoresT){
  mutations<-scoresT[which(scores$Variant_Classification!="Silent"),]
  colnames(mutations)<- colnames(scoresT)
  return (mutations)
}

#########################################################################################################################
#########################################################################################################################
##Function to get the hypermutated samples and an indicator of hypermutation status for TCGA data prior to PANCAN12.
getHypermutated<-function(rawdata,mutTable,mutationMatrix,metadata,ratio){
  
  silentMutations<-rawdata[which(rawdata$Variant_Classification=="Silent"),]
  silentAndNonSilentMutations <-cbind(table(silentMutations$sampleID), table(mutTable$sampleID))
  silentAndNonSilentMutations<-silentAndNonSilentMutations[order(silentAndNonSilentMutations[,2], decreasing=TRUE),]
  hypermutatedSampleNames<-rownames(silentAndNonSilentMutations[which(ratio>.18),])# I need this because the names are ordered differently in 
  #hypermutatedSampleNames<-rownames(silentAndNonSilentMutations[which(ratio>20),])# I need this because the names are ordered differently in 
  hyperIndex<-match(hypermutatedSampleNames,colnames(mutationMatrix))
  #hypermutated metadata variable
  hypermutatedMetaData<-rownames(metadata)%in%hypermutatedSampleNames# hypermutated samples are identified in the metadata2 table using this variable
  return(list(hyperIndex,hypermutatedMetaData))
  #return(hyperIndex)
}

##########################################################################################################################
##########################################################################################################################
## A function to assign samples a class and number ( for use in a color function) based on a metadata variable. It is used in an lapply fashion
## in order to create a list of samples by variable value and color value.

sampstr<-function(metadata, variableName){
  met<-metadata[,which(colnames(metadata)==variableName)]# the metadata variable
  groups<-unique(metadata[,which(colnames(metadata)==variableName)])
  colrs<-sapply(met, function(x) which(groups==x))# the colour indeces to be used to color plots points
  return(cbind.data.frame(colrs,met,stringsAsFactors=FALSE))
}

##########################################################################################################################
##########################################################################################################################


##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
## A plotting function to show functional mutation frequency, non-functional mutation frequency and an indelSNV ratio/ 
## indel percentage point.


mutFreqPlot<-function(freqData,indSNVratio, mutsPerIndiv){
  #plot
  par(mar=c(3,5,5,5))
  par(xpd=TRUE)
  barplot(t(freqData), col=c("blue", "red"),border=c("blue","red"), beside=FALSE,axisnames=FALSE, main= "non-silent mutation frequency and Indel:SNV ratio")# red = silent mutations and blue=nonsilent mutations
  par(new=T)
  #testing to see that the ordering of my stratification metadata for coloring plot points is correct.
  #identical(rownames(freqData), rownames(metadata2)[order(mutationsPerIndiv,decreasing=TRUE)])
  plot(indSNVratio, axes=F, pch=20, col=(stratification$hyperMutatedINDSNVRatio[order(mutationsPerIndiv, decreasing =TRUE),1]),cex=0.8,xlab="",ylab="", cex.axis=1.2) # a plot of the ratio of indels amongst all non-silent mutations
  #plot(indSNVratio, axes=F, pch=20,col=(stratification$vital_status[order(mutationsPerIndiv, decreasing =TRUE),1]), cex=0.8,xlab="",ylab="", cex.axis=1.2) # a plot of the percentage of indels amongst all non-silent mutations
  axis(4, pretty(c(0, max(indSNVratio))), pos= length(indSNVratio)+2)
  mtext("mutation frequency per sample", side=2, line=0, adj=0.5, padj=-3.5, cex=1.2)
  mtext("ratio of indels to SNV mutations",side=4,line=0, adj=0.5, padj=2.5, cex=1.2)
  mtext("samples ordered by mutation frequency",side=1, adj=0.5,padj=2.0, cex= 1.2)
  legend(0.8*length(indSNVratio),max(indSNVratio),bty="n", c("SNV","indel", "non-hyper", "hyper"), lty=c(1,1,NA,NA),lwd=c(2.5,2.5,NA,NA),pch=c(NA,NA,20,20),col=c("blue","red","black","red"), text.font=1.5)
  
}

##########################################################################################################################
##########################################################################################################################
# function to convert a table Object to a matrix object. Used when measuring mutation frequency per gene
tableToMatrix<- function(tableData){
  data<-as.matrix(tableData)
  class(data)<- "matrix"
  return(data)
}

##########################################################################################################################
##########################################################################################################################
# Function to return the mean number of mutations carried in each Individual for each gene.
getPerGenePerSamp<- function(PerIndpGene){
  #perIndivPerGene<-as.matrix(table(mutDataFrame$sampleID,mutDataFrame$Hugo_Symbol))
  meanPerIndivPerGene<- sapply(seq(1:dim(PerIndpGene)[2]), function(x) mean(PerIndpGene[,x]))
  return(meanPerIndivPerGene)
}

##########################################################################################################################
##########################################################################################################################
#
# Function to plot MDS results


MDSplotting<- function(mdsPoints,smplclrs,hyperMutIndex){
  par(mar=c(5,5,5,25))
  par(xpd=TRUE)
  plot(mdsPoints$points[,1],mdsPoints$points[,2], pch = 19, col=alpha(smplclrs[,1],0.5), cex=1.2, cex.lab= 1.5, cex.main = 1.5, cex.axis=1.5, xlab="dimension 1", ylab="dimension 2", main ="mds plot of colorectal cancer samples \n before network processing of mutation matrix")
  points(mdsPoints$points[hyperIndex,1], mdsPoints$points[hyperIndex,2],pch=1, cex = 2.0)
  #legend(0.5, 0.11, c("Illumina", "Solid","hypermutated"), cex=1., pch=c(19,19,1),col=c("red","blue","black"))
  legend(max(v$points[,1])+0.05, max(v$points[,2]),bty="n", c(unique(smplclrs[,2]),"hypermutated"), cex=1.5,pt.cex=1.2, pch=c(rep(19, length(unique(smplclrs[,2]))),1),col=c(seq_along(1:length(unique(smplclrs[,1]))),1))
}

##########################################################################################################################
##########################################################################################################################
#


