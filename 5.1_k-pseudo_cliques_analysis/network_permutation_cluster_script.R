#########################################################################################
#########################################################################################
library(igraph)
args<- commandArgs(TRUE)
#args<-c(3001,4001)
args<- as.numeric(args)
print(args)
# load the networks
network.path  <- "/home/rsutherlandbrc/networks/huppi2_noSL/"
network.name  <- "huppi2_noSL"
network.file  <- paste0(network.name,".simple",collapse=NULL)

#set the results path
setwd("/home/rsutherlandbrc/networks/huppi2_noSL/permutations/")

network<- scan(file = paste0(network.path, network.file, collapse=NULL), skip = 1L, what = list(nameA = character(), nameB = character()), sep = "\t")


##plot the subnetwork which an adjusted pvalue below 0.1.
edgeList<-do.call(cbind,network)
netgraph<-graph.edgelist(edgeList, directed=FALSE)

netgraph2<-netgraph#this is the graph in which I reassign the node names
V(netgraph2)$name<-rep("",length(V(netgraph)$name))
#######
## get the unique degrees

uniqueDegrees<-sort(unique(degree(netgraph, V(netgraph))))
##The genes which are of each degree, I can swap the names


######This is the degree contrained label shuffling function...
#This is the standard deviation of which I should shuffle node names because some nodes are of unique degree
degreeSD<- sd(degree(netgraph,V(netgraph)))

permDegrees<-list()
permGeneNames<-list()

NPerms<-1000


for(j in seq(1:NPerms)){
  #set a new random seed each time I permute the network
  usedNames<-vector()# the vector for keeping track ofthe genes which have already been selected 
 see<-j-1+args[1] 
print(see) 
 set.seed(j-1+args[1])# modify this so that for each permutation of the network the seed setting value increases by one

  #in order to get proper random swapping of the node labels I start from a different degree for each label swapping procedure
  shuffleUniqueDegrees<-sample(uniqueDegrees)
  for(i in seq_along(shuffleUniqueDegrees)){
    degreeN<-shuffleUniqueDegrees[i]
    
    genesofDegreeN<-which(degree(netgraph,V(netgraph))==degreeN)
    #namesGenesOfDegreeN<-names(genesofDegreeN)# the names of the genes to be shuffled
    for (k in seq_along(genesofDegreeN)){
      node2swap<- genesofDegreeN[k]#the node label I want to swap  
      
      #genesforselection to reassign node names
      
      genes4selection<-setdiff(V(netgraph)$name, usedNames)
      genes4selection<-setdiff(genes4selection, names(node2swap))
      if(length(genes4selection)==0){
        V(netgraph2)$name[node2swap]<-names(node2swap)
      }else{
        
        
        ###get the list of gene names which have not yet been assigned
        #unassignedNames<<-setdiff(V(netgraph)$name, usedNames)
        #unassignedNames<<-union(unassignedNames, names(genesofDegreeN)[k])
        #### the absolute d of unassigned genes to the current degree N
        
        
        
        absdiff<-sort(abs(degreeN-degree(netgraph,V(netgraph)[genes4selection])))
        names(absdiff)<-unlist(sapply(unique(absdiff), function(x)  sample(names(which(absdiff==x)))))
        
       # print(length(absdiff))
        if(length(absdiff)>1){
          absRankOrder<- seq(1:length(absdiff))
          names(absRankOrder)<-names(absdiff)
          absRankOrderP<-sapply(absRankOrder, function(x) pnorm(x, mean=0, sd=5, lower.tail=FALSE))
          
          newName<-names(sample(absRankOrder, 1,prob = absRankOrderP))
          
        }else{
          newName<-sample(names(absdiff), 1)
        }
        
        V(netgraph2)$name[node2swap]<- newName#This should be returned by the function
      }
      #add the newly assigned gene name to the list of those which can no longer be assigned.
      usedNames<<-c(usedNames,newName)
      
      #if the node to be swapped is not in the set of genes already asigned to nodes in the permutation, re-add node2swap to the genes available for selection
      #if(node2swap%in%usedNames!TRUE){
      #  genes4selection1<-c(names(node2swap), genes4selection)
      #}
      
    }
    #This checks if the neigborhoods for all of the genes of degreeN are identical in the original graph and for their name reassigned version in netgraph2
    #sapply(seq(1:length(usedNames)), function(x) identical(unlist(neighborhood(netgraph,1,names(genesofDegreeN)[x])),unlist(neighborhood(netgraph2,1,usedNames[x]))))
  }
  
  permGeneNames[[j]]<-V(netgraph2)$name
  
  write(permGeneNames[[j]],paste0("huppi2_noSL_Dconstd_M0_SD5_",args[1],"_",args[2],args[3],".txt"),sep="\t",ncolumns=3666, append=TRUE)
  
  #for creating the diference in degree histogram
  permDegrees[[j]]<-degree(netgraph2, V(netgraph2))
  
  #permDegrees[[j]]<-permDegrees[[j]]-degree(netgraph, V(netgraph))
  
}


save(permDegrees, file=paste("permDegrees_huppi2noSL_M0_SD5",args[1],args[2],args[3],".RData", sep = "_"))
