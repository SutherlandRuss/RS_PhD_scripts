#################################################################################
###sampleDiversity permutation
# a cluster script to identify if the sampleDiversity ( the number of samples mapping to a subnetwork) is higher in the original network than expected.

library(igraph)
library(hash)
args<- commandArgs(TRUE)
#args<-c(1,3,1000)
args<- as.numeric(args)

NPerms<-args[3]

# load the networks
network.path  <- "/home/rsutherlandbrc/networks/huppi2_noSL/"
network.name  <- "huppi2_noSL"
network.file  <- paste0(network.name,".simple",collapse=NULL)

network<- scan(file = paste0(network.path, network.file, collapse=NULL), skip = 1L, what = list(nameA = character(), nameB = character()), sep = "\t")

#all genes in the network
networkGeneList<-(sort.int(unique(c(network$nameA,network$nameB))))

##create the graph
edgeList<-do.call(cbind,network)
netgraph<-graph.edgelist(edgeList, directed=FALSE)

netgraph2<-netgraph#this is the graph in which I reassign the node names


#############################################################################################
##the Dense module enumeration load

densities<- c("0.95","0.90","0.85","0.8","0.75")
for(z in 1:length(densities)){
DME.density   <- paste0("a",densities[z])
DME.file      <- paste0(network.name,"_DME_results_",DME.density,".txt",collapse=NULL)
data.dir <- dirname(paste0("/home/rsutherlandbrc/networks/",network.name,"/",network.name,"_DME_results/",DME.file, collapse = NULL))
data.file<- basename (paste0("/home/rsutherlandbrc/networks/",network.name,"/",network.name,"_DME_results/",DME.file, collapse = NULL))

DME<- read.table(file.path(data.dir,data.file), header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE)
DME<- DME[which(sapply(1:length(DME[,1]),function(x) length(unlist(strsplit(DME[x,4], ","))))>3),]# removing modules with less than 3 nodes

# the integer module table I map this back to the permuted gene names
integerModuleTable<-strsplit(DME$V4, split=",", fixed = TRUE)

#set the results path
#setwd("/Users/rds/Google Drive/PhD/chapters/DME/RESULTS/")

dir.create(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density))
setwd(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density))

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


################################################################################################
####################### mapping the integer network node list to the igraph network node list
#The hash of netgraph 2 permuted gene names
nethash<-hash(keys=V(netgraph)$name, values=seq(1:length(V(netgraph))))
namehash<-hash(keys=networkGeneList, values= seq(1:length(networkGeneList)))

#netgraph node IDs for the networkGenelist genes in network genelist order
networkGeneListIndex<-sapply(seq(1:length(networkGeneList)), function(x) nethash[[networkGeneList[x]]])


########## nrow should be the number of permutations being run by the script
#the matrix to keep the module diversity for the permuted modules
permModules<- matrix(,nrow=NPerms, ncol= length(integerModuleTable))


######################################################################################
#######################################################################################

#load the permDegrees object

##make sure the file path is correct and the tart and end permutation numbers 1_1000 or 1001_2000 etc...
load(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/permutations/permDegrees_huppi2noSL_M0_SD5_",args[1],"_",args[2],"_.RData"))


#####the samplediversity permutation loop
### NPerms is the numer of permutations in the permDegrees object

for( m in seq(1:NPerms)){
  
  permutedNetworkGeneList<-names(permDegrees[[m]][networkGeneListIndex])
  #reassign networkgenelist names to be that of the new graph based on the PermDegrees names
  
  
  ##########the permutation test #################################      
###############################################################################
####### FINDS THE MEAN MUTSIGCV P-VALUE FOR EACH SUBNETWORK ###################
###############################################################################

  
  subnetworks<-lapply(lapply(integerModuleTable, as.integer),
                      function(i) permutedNetworkGeneList[i + 1L])# using the new permutedNetworkGenelist mapping
  

  #the mean mutsig score for the genes in the pseudoclique
  meanMutsig<-sapply(subnetworks, function(x) mean(scores[match(x, scores[,1]),2]))
  
    #the median mutsig score for the genes in the pseudoclique
#  medianMutsig<-sapply(subnetworks, function(x) median(scores[match(x, scores[,1]),2]))
 

 permModules[m,]<-meanMutsig
  ####Put the permuted module Diversity objects in to a matrix of subnetworks by N permutations
  ## get the value of the 50th (for 1000 perms) and compare the original value to that. If it is higher, then the prginal network has more mapped samples than we expect. 
  
}

#write out to a peruted sample diverity file which can be concatenated with the other permutation files
# make sure the 1_10 11_20 permutation numbers are correct.
write.table(permModules, paste0("permutedMutSigCV_",network.name,"_M0_SD5_",DME.density,"_",args[1],"_",args[2],".txt"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
##I need to sign the names back to the vertex names in the network.
##Then I need to redefine the gene sets so I can find the number of samples covered by a random network.











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




#########decide upon the k number of overlaps  to define the k-pseudoclique #########

setwd(paste0("/home/rsutherlandbrc/networks/huppi2_noSL/mutsig/",network.name,DME.density))
ks<- 3:7
##k  to define the k-communities
for(t in ks){
  k=t 
  #print(k)

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

  
  #the place I will keep the permuted sample diversity for each of the k communities
  kpermModules<- matrix(,nrow=NPerms, ncol=length(kcommunities))
  #I need to put my permutation loop here seq 1:nperms.
  
  
  #loop over the numer of permutations
  for (m in seq(1:NPerms)){
  permutedNetworkGeneList<-names(permDegrees[[m]][networkGeneListIndex])
  #reassign networkgenelist names to be that of the new graph based on the PermDegrees names
  
  #redefine subnetworks using the permuted Network gene list. This is vital forthe permutations
  subnetworks<-lapply(lapply(integerModuleTable, as.integer),
                      function(i) permutedNetworkGeneList[i + 1L])# using the new permutedNetworkGenelist mapping
  
  #The genes in each of the k communities
  kpcCommunities<-lapply(kcommunities, function(x) unique(unlist(subnetworks[c(which(kpcgraphClusters$membership==x))])))
  
  
  
  #################I can put this in to a function##############################
  
###############################################################################
####### FINDS THE MEAN MUTSIGCV P-VALUE FOR EACH SUBNETWORK ###################
##############################################################################

  #the mean mutsig score for the genes in the pseudoclique
  medianMutsig<-sapply(kpcCommunities, function(x) median(scores[match(x, scores[,1]),2]))

#  rTest<-unlist(lapply(kpcCommunities, function(x) paste(match(x, ranksGenesTest), collapse = "|")), use.names= FALSE)
  
  #add the above mutsigcV pvalues permutations test to the permutation
  kpermModules[m,]<- medianMutsig
  }
  
  #write out to a peruted sample diverity file which can be concatenated with the other permutation files
  # make sure the 1_10 11_20 permutation numbers are correct.
  write.table(kpermModules, paste0("permutedMutSigCV_huppi2noSL_M0_SD5_",DME.density,"_pseudo_",k,"_",args[1],"_",args[2],".txt"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  ##I need to sign the names back to the vertex names in the network.
  ##Then I need to redefine the gene sets so I can find the number of samples covered by a random network.

}


} #the outer densities loop


