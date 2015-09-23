#defines the library paths within which the libraries are found
#.libPaths(c("/Users/rds/Documents/R/win-library/3.0", "/Program Files/R/R-3.0.2/library"))

################################################################################
#### Thisscript can be used on the cluster to find the

library(igraph)
library(hash)
library(colorspace)
args<- commandArgs(TRUE)
#args<- as.numeric(args)

huppi<-TRUE

NPerms<-10000



#densities<- c("0.95","0.90","0.85", "0.8", "0.75","0.70")

density<- args[1] # the DME density
# load the networks
network.path  <- "/home/rsutherlandbrc/networks/huppi2_noSL/"
network.name  <- "huppi2_noSL"
network.file  <- paste0(network.name,".simple",collapse=NULL)

DME.density   <- paste0("a",density)
DME.file      <- paste0(network.name,"_DME_results_",DME.density,".txt",collapse=NULL)
data.dir <- dirname(paste0("/home/rsutherlandbrc/networks/",network.name,"/",network.name,"_DME_results/",DME.file, collapse = NULL))
data.file<- basename (paste0("/home/rsutherlandbrc/networks/",network.name,"/",network.name,"_DME_results/",DME.file, collapse = NULL))


  DME<- read.table(file.path(data.dir,data.file), header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE)
  DME<- DME[which(sapply(1:length(DME[,1]),function(x) length(unlist(strsplit(DME[x,4], ","))))>3),]# removing modules with less than 3 nodes

  # the integer module table I map this back to the permuted gene names
  integerModuleTable<-strsplit(DME$V4, split=",", fixed = TRUE)

dir.create(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density,"/concatenated"))
setwd(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density,"/concatenated"))


  ########################################################################

  ########################################################################
  ###loading the list of genes which are mutated in any one sample
  ##make sure the path is correct for this
load("/home/rsutherlandbrc/networks/huppi2_noSL/mutationData/RanksGenesTest.RData")
load("/home/rsutherlandbrc/networks/huppi2_noSL/mutationData/mutationTableTest.RData")


###############
samplesTest <- unique(unlist(mutationTableTest,use.names=FALSE))
#samplesTest <- samples[sort.list(as.integer(sub("Apollo3_", "", samplesTest, fixed = TRUE)))]


#########################################################################################
#########################################################################################
# load the networks

###DME.density   <- "a0.95"
##DME.file	<- paste0(network.name,"_DME_results_",DME.density,".txt",collapse=NULL)
##data.dir <- dirname(paste0("/Users/rds/Google Drive/PhD/chapters/DME/Ranking_project/Gene_Ranking/DME/RESULTS/",network.name,"/",DME.file, collapse = NULL))
##data.file<- basename (paste0("/Users/rds/Google Drive/PhD/chapters/DME/Ranking_project/Gene_Ranking/DME/RESULTS/",network.name,"/",DME.file, collapse = NULL))
##enrich.name<-paste0(network.name,"_enrichment_",DME.density, collapse=NULL)
##enrich.path<-dirname(paste0("/Users/rds/Google Drive/PhD/chapters/DME/Ranking_project/Gene_Ranking/DME/RESULTS/",enrich.name,"/",collapse=NULL))
##enrich.file<-basename(paste0("/Users/rds/Google Drive/PhD/chapters/DME/Ranking_project/Gene_Ranking/DME/RESULTS/",enrich.name,"/",collapse=NULL))



##DME<- read.table(file.path(data.dir,data.file), header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE)
##DME<- DME[which(sapply(1:length(DME[,1]),function(x) length(unlist(strsplit(DME[x,4], ","))))>3),]# removing modules with less than 3 nodes
# as these are automatically a clique and canot be anything else
network<- scan(file = paste0(network.path, network.file, collapse=NULL), skip = 1L, what = list(nameA = character(), nameB = character()), sep = "\t")

#all genes in the network
networkGeneList<-(sort.int(unique(c(network$nameA,network$nameB))))
integerModuleTable<-strsplit(DME$V4, split=",", fixed = TRUE)



########################################################################
###loading the list of genes which are mutated in any one sample
##make sure the path is correct for this
scores.file<-"/home/rsutherlandbrc/networks/huppi2_noSL/mutationData/colorectal_non_hyper_mutsig_perm_data.txt"

scores<-read.table(file=scores.file, sep="\t",header=TRUE, stringsAsFactors = FALSE)



##################This finds the genes which are in the network but are not scored by MutsigCV

notInMutSig<-setdiff(networkGeneList,scores[,1])
## to each of the proteins not in the mutsig file assign them to the mean p-value of the mutsigcv file
notInMutSig<-data.frame(gene=notInMutSig,p=rep(mean(scores[,2])))
## then combine the scores object and the notInMutsig objects


##This is the scores object with the genes in the network, but the mutsig data added in
scores<-as.data.frame.matrix(rbind(scores,notInMutSig))





## dme index starts at 0!!!
# It is vitally important to use "subnetworks" variable when creating
# subgraphs as igraph will search for the node names within the huppi igraph
# object.

subnetworks<-lapply(lapply(integerModuleTable, as.integer),
                    function(i) networkGeneList[i + 1L])


 ##########the permutation test #################################
###############################################################################
####### FINDS THE MEAN MUTSIGCV P-VALUE FOR EACH SUBNETWORK ###################
###############################################################################



  #the mean mutsig score for the genes in the pseudoclique
#  meanMutsig<-sapply(subnetworks, function(x) mean(scores[match(x, scores[,1]),2]))

   #the median mutsig score for the genes in the pseudoclique
  medianMutsig<-sapply(subnetworks, function(x) median(scores[match(x, scores[,1]),2]))


###################################################################################################################
######load the permutations
###################################################################################################################
permutations<- read.table(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density,"/concatenated/","permutedMutSigCV_",network.name,"_",DME.density,"_",NPerms), header=FALSE, sep= "\t", stringsAsFactors = FALSE)



##order each column
permutationsSorted <- apply(permutations,2, sort, decreasing = F)

##permutations median
permutationsmedian<- round(apply(permutations,2, median),digits = 2)

##empirical pvalue
empirical_value<-sapply(seq(1:length(subnetworks)), function(x) length(which(permutationsSorted[,x]<=medianMutsig[x])))/NPerms# N permutations

##finding out the index of the subnetwork which is mutated in more samples than you expect given the network
sigIndex<- which(empirical_value<0.05)

########################################################################################################
## Data Ouput
########################################################################################################


# the index position of the elements of the empirical p values arranged by empirical p-values
i= sort.list(empirical_value, decreasing = FALSE)

alldataTest<- cbind(empirical_value,medianMutsig,permutationsmedian,sapply(subnetworks, paste, collapse = "|"))[i,]

write.table(alldataTest,file = paste0(network.name,DME.density,"perm_test.txt",collapse =NULL), append = FALSE, quote= FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

#dir.create(paste("/Users/rds/Google Drive/PhD/chapters/DME/RESULTS/",enrich.file, sep=""))
#setwd(paste("/Users/rds/Google Drive/PhD/chapters/DME/RESULTS/",enrich.file, sep=""))

getGenelistOutput<- function(sigIndex, subnetworks){
  #x<-(seq_along(ll$subnetworks[[77]]))
  x<-lapply(sigIndex, function(z) seq_along(subnetworks[[z]]))
  #significantGenelist<- split(ll$subnetworks[[77]], ceiling(x/5))
  significantGenelist<-lapply(1:length(sigIndex), function(z) split(subnetworks[[sigIndex[z]]], ceiling(x[[z]]/10)))# change the 10 here in order to change the number of genes on each line

  return(significantGenelist)
}






#the percentage of samples that are mutated in each gene
samplemutationpcntage<-sapply(seq(1:length(mutationTableTest)), function(x) length(mutationTableTest[[x]])/length(samplesTest)*100)
names(samplemutationpcntage)<-names(mutationTableTest)

#the break to which each sample belongs. This plus one is the index of the color the gene should be coloured in a network plot
#I have breaks between 0,1,2,3,4,5,6,7,8,9,20,30,40,50,60,70,80,90,100.
#There are 19 levels so I need to create a palette with 19 colors
samplemutationpcntageBreaks<-cut(samplemutationpcntage,breaks=c(0,1,2,3,4,5,seq(6,10),seq(20,100,10)), labels=c(0,1,2,3,4,5,seq(6,10),seq(20,90,10)))
#define the palette
palette1<-rev(terrain_hcl(length(levels(samplemutationpcntageBreaks)),h=c(3,180), c = c(80, 20), l = c(60,90), power = c(1/2,0.2)))



samplemutationpcntageBreaks<-as.numeric(samplemutationpcntageBreaks)
samplemutationpcntageBreaks+1# the index for the color in my palette1
names(samplemutationpcntageBreaks)<-names(mutationTableTest)

proteinColorMapping<-palette1[samplemutationpcntageBreaks]
names(proteinColorMapping)<-names(mutationTableTest)# I now simply match the node names to this vector to find the node's color


####calculate the five percent mutated samples
#Do this if I need to filter the genes by 5% sample mutation percentage
FivePercent<-ceiling(length(samplesTest)/20)

samplemutationcount<-sapply(seq(1:length(mutationTableTest)), function(x) length(mutationTableTest[[x]]))
names(samplemutationcount)<- names(mutationTableTest)
mutationMatrixProteins<-names(which(samplemutationcount>FivePercent))





##plot the subnetwork which an adjusted pvalue below 0.1.
if(huppi){
  edgeList<-do.call(cbind,network)
}else{
  edgeList<-lapply(seq(1:dim(network)[1]), function(x) unlist(network[x,]))
  edgeList<-do.call(rbind,edgeList)
}
netgraph<-graph.edgelist(edgeList, directed=FALSE)





if (length(sigIndex)<1){}else{


  geneListOutput<-getGenelistOutput(sigIndex, subnetworks)
  sink("significantGeneList", append=FALSE)
  print("Each significant gene list gets it's own highest order list, the length of which is equal to the number of significant gene set.")
  print(geneListOutput)
  sink()

  for(i in 1:length(sigIndex)){
    j<-sigIndex[i]

    graph1<-induced.subgraph(netgraph,subnetworks[[j]])
    V(graph1)
 ##change the density at line 573
    maxColorValue <- length(samplesTest)

    palette <- colorRampPalette(c("light blue","red"))(maxColorValue)
    ##plotcolors
    mappedcolrs<-proteinColorMapping[V(graph1)$name]
    names(mappedcolrs)<-V(graph1)$name
    mappedcolrs[which(is.na(mappedcolrs))]<-"#FFFFFF"#NA values are mapped to white
  
    V(graph1)$color<-mappedcolrs  #assigning the vertex colors to the vertices
    pdf(file=paste("significantnetwork",i,length(mappedcolrs),"nodes.pdf", sep="_"),width=7, height=7)
    par(mai=c(0,0,1,0))
    plot.igraph(graph1, vertex.frame.color="light gray",vertex.label.family="sans", vertex.label.color= "#545454", layout=layout.fruchterman.reingold,vertex.size=13, vertex.label.cex=0.1, 
edge.color="light gray")
    dev.off()
  }#loop for the number of significant subnetworks
}#the if else huppi condition







#create pseudoclique - pseudoclique overlap matrix and group the cliques in order to find a clustering of the cliques
#this should reduce the number of pseudo cliques
Nsubnetworks<-length(subnetworks)
seqsubnetworks<- seq(1,Nsubnetworks)

################################################################################################
#################################################################################################
####k-pseudocliques


psC<-lapply(seqsubnetworks, function(y) sapply(seqsubnetworks, function(x) c(y,x)))


#check if this is the same as the overlapmatrix
overlapMatrix<-sapply(seq(1:Nsubnetworks), function(y) sapply(seq(1:Nsubnetworks), function(x) length(intersect(subnetworks[[psC[[y]][1,x]]], subnetworks[[psC[[y]][2,x]]]))))
#for(i in 1:Nsubnetworks){
#  overlapMatrix[i,i]<-0
#}




pseudocgraph<- overlapMatrix
pseudocgraph<-as.matrix(pseudocgraph)

##this is the graph of pseudocliques
pcgraph<- graph.adjacency(pseudocgraph,mode="undirected",weighted=TRUE)
V(pcgraph)$label <- seq(vcount(pcgraph))
pcgraphClusters<-clusters(pcgraph, mode="strong")#identifies the components of the graph
#each separate component should be considered as an independent entity because non of the nodes overlap.
#in larger components we should decompose them in to either k-pseudo cliques using a graph decomposition, or k clusters using a distance matrx

#plot(pc)
pcgraphs<-decompose.graph(pcgraph)#decompose pc graph in to separate components. components composed of less than 2 should be considered independent groups
#size of each of the components in test graph
largest <- which.max(sapply(pcgraphs, vcount))





ks<- 3:7
##k  to define the k-communities
for(k in ks){
  print(k)

  ## I can reduce the pseudocliquegraph to that of only the pseudocliques I am interested in, those of size k and larger
  ## This means I must set to zero, any overlaps that are smaller than k.

  kpseudocgraph<-pseudocgraph
  kpseudocgraph[which(kpseudocgraph<k-1)]<-0

  for(i in 1:Nsubnetworks){

    #if(kpseudocgraph[i,i]<k){
    kpseudocgraph[i,i]<-0
    #}
  }



  #graph object for the k pseudoclique graph
  kpcgraph<- graph.adjacency(kpseudocgraph,mode="undirected",weighted=TRUE)
  V(kpcgraph)$label <- seq(vcount(kpcgraph))

  kpcgraphClusters<-clusters(kpcgraph, mode="strong")#identifies the components of the graph


  kpcgraphs<-decompose.graph(kpcgraph)#decompose pc graph in to separate components. components composed of less than 2 should be considered independent groups


  ####assign the pseudocliques to the pseudoclique communities####
  #the number of communities
  ncommunities<-kpcgraphClusters$no
  kcommunities<- which(kpcgraphClusters$csize>1)


 #The genes in each of the k communities
  kpcCommunities<-lapply(kcommunities, function(x) unique(unlist(subnetworks[c(which(kpcgraphClusters$membership==x))])))


  
  ###############################################################################
####### FINDS THE MEAN MUTSIGCV P-VALUE FOR EACH SUBNETWORK ###################
##############################################################################

  #the mean mutsig score for the genes in the pseudoclique
  medianMutsig<-sapply(kpcCommunities, function(x) median(scores[match(x, scores[,1]),2]))


  #add the above mutsigcV pvalues permutations test to the permutation


  #write out to a peruted sample diverity file which can be concatenated with the other permutation files
  # make sure the 1_10 11_20 permutation numbers are correct.

  ###################################################################################################################
  ######load the permutations
  ###################################################################################################################
permutations<- read.table(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density,"/concatenated/","permutedMutSigCV_huppi2noSL_",DME.density,"_pseudo_",k,"_",NPerms),
header=FALSE, sep= "\t", stringsAsFactors = FALSE)

 ##order each column
  permutationsSorted <- apply( permutations,2, sort, decreasing = F)

  ##permutations median
  permutationsmedian<- round(apply(permutations,2, median),digits = 2)

  ##empirical pvalue
  empirical_value<-sapply(seq(1:length(kcommunities)), function(x) length(which(permutationsSorted[,x]<=medianMutsig[x])))/NPerms# N permutations

  ##reassign empirical pvalues that are 0 to be 0.000099, just below what can be measured using 10000 permutations, but a positive value
  ##so that the FDR calculation will return a useful qvalue.
  empirical_value[empirical_value==0]<-0.000099

#empirical FDR of the empirical	p-values

empirical_FDR<-p.adjust(empirical_value,method="BH")


  ##finding out the index of the subnetwork which is mutated in more samples than you expect given the network
  sigIndex<- which(empirical_FDR<=0.1)

 ########################################################################################################
  ## Data Ouput
  ########################################################################################################

  # the index position of the elements of the empirical p values arranged by empirical p-values
  i= sort.list(empirical_value, decreasing = FALSE)

  alldataTest<- cbind(empirical_value,empirical_FDR,medianMutsig,permutationsmedian,sapply(kpcCommunities, paste, collapse = "|"))[i,]




#create the results folders for the k cliques

dir.create(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density,"/concatenated/",k,"_clique_comms", sep=""))
setwd(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density,"/concatenated/",k,"_clique_comms", sep=""))


  write.table(alldataTest,file = paste0(network.name,DME.density,"_",k,"_pseudo_perm.txt",collapse =NULL), append = FALSE, quote= FALSE, sep = "\t",
              row.names = FALSE, col.names = TRUE)
  #write the community gene sets
  #lapply(kpcCommunities, cat, "\n", file=paste0(enrich.file,"_",k,"_communities.txt", collapse=NULL), append=TRUE)


#write the graph clusters
  #lapply(kpcgraphClusters, cat, "\n", file=paste0(enrich.file,"_",k,"_graph_clustes.txt", collapse=NULL), append=TRUE)
  save(kpcgraphClusters,file=paste0(network.name, DME.density,"_",k,"_graph_clusters.Rdata", collapse=NULL))
  ###loadObject(file=paste0(enrich.file,"_",k,"_graph_clustes.Rbin", collapse=NULL))


  library(R.utils)
  library(ggplot2)
  #### The results do not show any significant enrichment of any subnetworks. Due to the large sample size most of the genes carry mutations. Therefore networks are not enriched for mutations.
  #### Does it make sense to look for dense subnetworks which are covered by the largest number of samples?
  #### Can I do this by permuting the gene labels within degree and look for subnetworks which are covered by the largest number of samples?
  #### Or is there a way to find the significantly mutated genes (nature paper) and find subnetworks which are enriched for significantly mutated genes and then look at the p-value of the non-significantly mut$


  ##plot the subnetwork which an adjusted pvalue below 0.1.
  #edgeList<-do.call(cbind,network)
  #netgraph<-graph.edgelist(edgeList, directed=FALSE)





######### create a new circle object for igraph ##########


mycircle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {        vertex.color <- vertex.color[v]
    }
    vertex.size  <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
        vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
        vertex.frame.width <- vertex.frame.width[v]
    }

    mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
           vertex.size, vertex.frame.width,
           FUN=function(x, y, bg, fg, size, lwd) {
               symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                       circles=size, add=TRUE, inches=FALSE)
           })
}
add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                 plot=mycircle, parameters=list(vertex.frame.color=1,
                                  vertex.frame.width=1))



##############################################################


mutsigGenes<-c("KRAS","NRAS","TP53","APC","SMAD4","FBXW7","SMAD2","BRAF","TNFRSF10C","TCF7L2","C17orf97","CTNNB1")




  if (length(sigIndex)<1){}else{

    #print the significant gene lists
    geneListOutput<-getGenelistOutput(sigIndex, kpcCommunities)
    sink(paste0("significantGeneList_",k,"txt", collapse = NULL), append=FALSE)
    print("Each significant gene list gets it's own highest order list, the length of which is equal to the number of significant gene lists. Lower order lists of length 10 indicate a partial list of a significant gene set")
    print(geneListOutput)
    sink()


    for(i in 1:length(sigIndex)){
      j<-sigIndex[i]

      graph1<-induced.subgraph(netgraph,kpcCommunities[[j]])
      V(graph1)

      ##change the density at line 573
      maxColorValue <- length(samplesTest)

      palette <- colorRampPalette(c("light blue","red"))(maxColorValue)
      ##plotcolors
      mappedcolrs<-proteinColorMapping[V(graph1)$name]
      names(mappedcolrs)<-V(graph1)$name
      mappedcolrs[which(is.na(mappedcolrs))]<-"#FFFFFF"#NA values are mapped to white
      V(graph1)$color<-mappedcolrs  #assigning the vertex colors to the vertices
  
      mutatedinMutsig<-match(V(graph1)$name,mutsigGenes)#find the proteins mutated in more than 5pcnt of samples
      mutatedinMutsigCol<-mutatedinMutsig
      mutatedinMutsig[which(!is.na(mutatedinMutsig))]<-5#assign the vertex border width
      mutatedinMutsig[which(is.na(mutatedinMutsig))]<-1
      V(graph1)$vertex.frame.width<-mutatedinMutsig

      mutatedinMutsigCol[which(!is.na(mutatedinMutsigCol))]<-"#feb24c"
      mutatedinMutsigCol[which(is.na(mutatedinMutsigCol))]<-"light gray"

 pdf(file=paste("significantnetwork",i,length(mappedcolrs),"nodes.pdf", sep="_"),width=7, height=7)
      par(mai=c(0,0,1,0))
plot.igraph(graph1, vertex.frame.color= mutatedinMutsigCol,vertex.shape="fcircle",vertex.frame.width= mutatedinMutsig,
vertex.label.family="sans", vertex.label.color= "#545454", layout=layout.fruchterman.reingold, vertex.size=13,
 vertex.label.cex=1, edge.color="light gray")
      dev.off()
    }# loop over all significant genesets and plot the graph
  }# if there are any significant genesets loop #k =ks[j] loop
}




