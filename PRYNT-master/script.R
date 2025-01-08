########
#Script prioritization methods
#######
#This script describe prioritization methods applied on a PPI network using a urinary proteomics profiles from renal diseases 
#published  : 

#I illustrate my workflow here with the example of urinary profile from ADPKD
#published  : Bakun M, Niemczyk M, Domanski D, et al. Urine proteome of autosomal dominant polycystic kidney disease patients. Clin Proteomics. 2012;9(1):13. Published 2012 Dec 11. doi:10.1186/1559-0275-9-13


####Packages####
library(igraph)
library(zoo)
library(CTDquerier)
library(RandomWalkRestartMH)

####Functions####

shortest_path_ranking_function<-function(graph,seed,score=NULL,directed=T,distance_table=NULL){
  #calculate the closeness of each protein from seeds using the shortest path
  #graph : igraph object, containing the network to analyse
  #seed : character vector of DE proteins
  #score : NULL if we don't take weight in to account, or numeric vector containing weight of each relations (more reliable relations have a minimal score )
  #directed : boolean
  #distance_table : to have the possibility to calculate the distance table outside of the function
  
  if(is.null(distance_table)){
    if(directed){sens="out"}
    else{sens="all"}
    distance_table<-distances(graph,v = V(graph),to = V(graph), mode = sens,weights = score)
  }
  
  distance_table_to_seed<-distance_table[,colnames(distance_table)%in%seed]
  #0 means that the protein is DE. I put the distance of a DE to its self to 1.
  distance_table_to_seed[which(distance_table_to_seed==0,arr.ind = T)]<-Inf
  
  results_SP<-apply(X = distance_table_to_seed,MARGIN = 1,function(x){sum(1/x)})
  #resultat_closeness vector of adapted closeness centrality
  results_SP<-data.frame("node"=names(results_SP),"SP"=results_SP)
  results_SP<-results_SP[order(results_SP$node),]
  return(results_SP)
}

random_walk_ranking_function<-function(graph,seed,adjacency_matrix=NULL,multiplex_object=NULL){
  #calculate the closeness of each protein from seeds using the random walk 
  # using the RandomWalkRestartMH package
  #graph : igraph object, containing the network to analyse
  #seed : character vector of DE proteins
  
  #adjacency_matrix : to have the possibility to calculate the adjency matrix outside of the function
  #multiplex_object : to have the possibility to calculate the multiplex object outside of the function
  
  cat("calss of graph  : ", class(graph))
  if(is.null(adjacency_matrix)){
    multiplex_object <- create.multiplex(graph,Layers_Name=c("PPI"))
  }
  if(is.null(adjacency_matrix)){
    adjacency_matrix <- compute.adjacency.matrix(multiplex_object)
  }
  adjacency_matrix_norm <- normalize.multiplex.adjacency(adjacency_matrix)
  seed<-seed[seed%in%names(V(graph))]
  # We launch the algorithm with the default parameters (See details on manual)
  RWR_PPI_Results <- Random.Walk.Restart.Multiplex(adjacency_matrix_norm,
                                                   multiplex_object,seed)
  # We display the results
  results_random_walk<-RWR_PPI_Results$RWRM_Results
  colnames(results_random_walk)<-c("node","RW")
  node<-data.frame("node"=names(V(graph)))
  results_random_walk<-merge(x = node,y =results_random_walk,by="node",all.x = T )
  results_random_walk$RW[is.na(results_random_walk$RW)]<-0
  results_random_walk<-results_random_walk[order(results_random_walk$node),]
  
  return(results_random_walk)
}

direct_ranking_function<-function(graph,seed){
  #graph : igraph object, containing the network to analyse
  #seed : character vector of DE proteins
  edgelist<- as.data.frame(get.edgelist(graph),stringsAsFactors = F)
  
  upstream_neighbour<-edgelist[edgelist[,2]%in%seed,1]
  results_direct<-data.frame(table(upstream_neighbour),stringsAsFactors =F )
  results_direct<-results_direct[order(results_direct$Freq,decreasing = T),]
  colnames(results_direct)<-c("node","D")
  node<-data.frame("node"=names(V(graph)))
  results_direct<-merge(x = node,y =results_direct,by="node",all.x = T )
  results_direct$D[is.na(results_direct$D)]<-0
  results_direct<-results_direct[order(results_direct$node),]
  
  return(results_direct)
}

precision_plot_single<-function(precision,sequence,col){
  y<-precision
  x<-seq(0,1,length.out = length(y))
  AUC <- sum(diff(x)*rollmean(y,2))
  x<-sequence
  plot(x,y*100,pch=19,xlim = c(sequence[1],sequence[length(sequence)]),ylim=c(0,100),col=col,
       xlab = "Nombre de protéines",ylab="% de proteines adpkd")
  cord.x <- c(sequence[1],x,sequence[length(sequence)])
  cord.y <- c(0,y,0)
  rgb_color<-col2rgb(col)
  color<-rgb(rgb_color[1,1]/255,rgb_color[2,1]/255,rgb_color[3,1]/255,alpha=0.3,maxColorValue = 1) 
  polygon(cord.x,cord.y*100,col=color,density=NA)
  lines(x,y*100,type="l",pch=19,col=col)
  text(x =sequence[6],y=0.1*100,paste("AUC=",round(AUC,3)),col=col,cex = 1.8)
}


precision_function<-function(binary_ranking,step){
  sequence<-seq(from = step,to = length(binary_ranking),by = step)
  precision<-NULL
  for(i in sequence){
    precision<-c(precision,sum(binary_ranking[1:i])/length(binary_ranking[1:i]))
  }
  return(precision)  
}

precision_plot<-function(list_binary_results,step=5,size_top=100,legend=T,
                         color=rainbow(length(list_binary_results)),last=NULL){
  
  AUC_vector<-rep(0,length(list_binary_results))
  names(AUC_vector)<-names(list_binary_results)
  sequence<-1:length(list_binary_results)
  
  if(!is.null(last)) sequence<-c(sequence[-last],last)
  for(j in sequence){
    
    y<-precision_function(binary_ranking=list_binary_results[[j]][1:size_top],step)
    x<-seq(0,1,length.out = length(y))
    AUC_vector[j] <- sum(diff(x)*rollmean(y,2))
    if (length(list_binary_results)==1 & legend){
      precision_plot_single(precision =y,sequence = seq(step,size_top,step) ,col=color)
    }
    else{
      if(j==sequence[1]){
        plot(seq(step,size_top,step),y*100,pch=19,xlim = c(step,size_top),ylim=c(0,100),col=color[j],
             xlab = "Nombre de protéines",ylab="% de proteines adpkd",type="l")
      }
      else  {
        lines(seq(step,size_top,step),y*100,type="l",pch=19,col=color[j])
        
      }
      
      if(length(list_binary_results)<=5){
        lines(seq(step,size_top,step),y*100,type="p",pch=19,col=color[j])
        
      }
      if(legend){text(x=80, y = (100-10*j), labels = paste(names(list_binary_results)[[j]],"(auc=", round(AUC_vector[j],3),")"),
                      col=color[j],cex = 0.8)}
      
    }
  }
  res<-AUC_vector[order(AUC_vector,decreasing = T)]
  return( AUC_vector[order(AUC_vector,decreasing = T)])
  
  
}


####Network####
#Construction of the protein-protein interactions network from string database
#adress : https://string-db.org/cgi/download.pl?sessionId=0xLJ19QRr4LW
setwd("D:/olhumpus/string/")
network_actions <-read.table(file ="9606.protein.actions.v10.5.txt",header = T,sep = "\t",stringsAsFactors = F ) 

#annotation#
#details of each proteins. datas obtain with string database available on github
#annotation_string permit to obtain  to convert Ensembl name to protein names.
annotation<-read.table(file ="annotation_string2.txt",header = T,sep = "\t",stringsAsFactors = F,quote="" )
annotation<-annotation[,1:2]

#replace Ensembl names by regular protein names
network_actions<-merge(x = network_actions, y = annotation,by.x="item_id_a",by.y="protein_external_id")
network_actions<-network_actions[,-which(colnames(network_actions)=="item_id_a")]
colnames(network_actions)[which(colnames(network_actions)=="preferred_name")]<-"item_id_a"

network_actions<-merge(x = network_actions, y = annotation,by.x="item_id_b",by.y="protein_external_id")
network_actions<-network_actions[,-which(colnames(network_actions)=="item_id_b")]
colnames(network_actions)[which(colnames(network_actions)=="preferred_name")]<-"item_id_b"


#Select network with choosen parameters
network_actions_directed<-network_actions[network_actions$is_directional=="t",
                                          -which(colnames(network_actions)=="is_directional")]
network_actions_directed_acting<-network_actions_directed[network_actions_directed$a_is_acting=="t",
                                                          -which(colnames(network_actions_directed)=="a_is_acting")]

network_actions_directed_acting_900<-network_actions_directed_acting[which(network_actions_directed_acting$score>=900),
                                                                     -which(colnames(network_actions_directed_acting)=="score")]

#I don't use the action and mode columns
network_actions_directed_acting_900<-network_actions_directed_acting_900[,c("item_id_a","item_id_b")]


#remove redundant relations and self loops
network_actions_directed_acting_900<-unique(network_actions_directed_acting_900)
sum(network_actions_directed_acting_900$item_id_a==network_actions_directed_acting_900$item_id_b)

proteins_actions_directed_acting_900<-unique(c(network_actions_directed_acting_900$item_id_a,network_actions_directed_acting_900$item_id_b))

head(network_actions_directed_acting_900)
network<-network_actions_directed_acting_900

####experimental datas####

#list_DE_bakun2012 character vector of Differently Expressed proteins (DE) from Bankun et al. (2012) concerning ADPKD  available on github
#adress : https://github.com/Boizard/Prioritization_urinary_proteome/
#this list can be replace by your own dataset
DE_proteins<-read.table(file = "list_DE_bakun2012.csv",
                                     stringsAsFactors = F,sep = '\t',header = F)
DE_proteins<-DE_proteins[,1]




####Contexctualization of the network####
#add relations with DE_proteins under a score of 900 to the network
ind1<-network_actions_directed_acting$item_id_a%in%DE_proteins & 
  network_actions_directed_acting$item_id_b%in%proteins_actions_directed_acting_900

ind2<-network_actions_directed_acting$item_id_b%in%DE_proteins & 
  network_actions_directed_acting$item_id_a%in%proteins_actions_directed_acting_900

sup_relation<-network_actions_directed_acting[ind1|ind2,]
sup_relation<-sup_relation[,c("item_id_a","item_id_b")]
sup_relation<-unique(sup_relation)
network<-rbind(network_actions_directed_acting_900,sup_relation)
network<-unique(network)
proteins_network<-unique(c(network$item_id_a,network$item_id_b))


#The network is transform as an igraph object
graph<-graph_from_data_frame(network)
proteins_network<-names(V(graph))
number_proteins_network<-length(proteins_network)
DE_proteins_network<-DE_proteins[DE_proteins%in%proteins_network]

####Cliques calculation#######
#I use max_clique from igraph package
clique<-max_cliques(graph = graph,min=3)

#Transform the list to a matrix
clique_matrix<-matrix(data = "",nrow = length(clique),ncol =2 )
clique_matrix[,1]<-as.character(1:length(clique))
for( i in 1:length(clique)){
  clique_matrix[i,2]<-paste(names(clique[[i]]),sep= "",collapse = " / ")
  
}
#clique_matrix is a matrix, each line is a clique, and the second colum is the list of proteins is the clique

#I selection only the bigger clique for each protein
clique_max<-rep(0,length(V(graph)))
names(clique_max)<-names(V(graph))
clique_bigger<-rep(0,length(V(graph)))
names(clique_bigger)<-V(graph)

for (i in 1:length(clique)){
  longueuri<-clique_max[clique[[i]]]
  for( j in 1:length(longueuri)){
    if(longueuri[j]<length(clique[[i]])){
      clique_max[clique[[i]][j]]<-length(clique[[i]])
      clique_bigger[clique[[i]][j]]<-i}
    else if(longueuri[j]==length(clique[[i]])){
      clique_bigger[clique[[i]][j]]<-paste(clique_bigger[clique[[i]][j]],i,sep=" / ")
    }
  }
}

clique_bigger<-data.frame("node"=names(V(graph)),clique_bigger,"size_clique"=clique_max,stringsAsFactors = F)
#clique_bigger is a matrix, is line is a protein of the network
#2nd lines : bigger clique which the proteins is
#3rd lines : size of the clique


####Build a new network with clustered proteins####

clique_data<-clique_bigger$clique_bigger


#I put one condition clique have to be bigger than 2
#I accept that the "clique" is a group of several clique
clique_size<-table(clique_data)
clique_bigger_than_2<-(clique_data%in%names(clique_size[clique_size>2]))

#clustered_protein select only proteins which could be clustered in their bigger cliques
clustered_protein<-clique_bigger[clique_bigger_than_2&clique_bigger$clique_bigger!=0 ,]

#replace proteins names in the network by their clique numbers (ex: clique12)

clustered_clique<-unique(as.character(clustered_protein$clique_bigger))
network_clique<-network
for( i in 1:length(clustered_clique)){
  node_clique<-clustered_protein[clustered_protein$clique_bigger%in%clustered_clique[i],'node']
  network_clique[network_clique$item_id_a%in%node_clique,1]<-paste('clique',clustered_clique[i],sep = "")
  network_clique[network_clique$item_id_b%in%node_clique,2]<-paste('clique',clustered_clique[i],sep = "")
  
}

#A lot of relations are now redundant
#remove redundant relations
network_clique<-unique(network_clique)
#remove self-loops
network_clique<-network_clique[network_clique$item_id_a!=network_clique$item_id_b,]


#network  with clique
graph_clique<-graph_from_data_frame(network_clique)

protein_network_clique<-names(V(graph_clique))
number_protein_network_clique<-length(protein_network_clique)

#Selection of the DE proteins in network_clique
#I consider a clique as DE if at least one protein is DE in the clique
DE_proteins_network_clique<-DE_proteins

for( i in 1:length(DE_proteins_network_clique)){
  if(DE_proteins_network_clique[i]%in%clustered_protein$node){
    DE_proteins_network_clique[i]<-paste('clique',
                                         clustered_protein[clustered_protein$node==DE_proteins_network_clique[i],2],sep = "")
    
  }}

DE_proteins_network_clique<-unique(DE_proteins_network_clique)

####Calculation of ranking in the network####
results_SP<-shortest_path_ranking_function(graph = graph,seed = DE_proteins_network)
results_SP<-data.frame(results_SP,"rank_SP"=rank(-results_SP$SP))

results_RW<-random_walk_ranking_function(graph = graph,seed = DE_proteins_network)
results_RW<-data.frame(results_RW,"rank_RW"=rank(-results_RW$RW))

results_D<-direct_ranking_function(graph = graph,seed = DE_proteins_network)
results_D<-data.frame(results_D,"rank_D"=rank(-results_D$D))

head(results_SP)
head(results_RW)
head(results_D)



#Complete results on the PPI network
results_PPI<-merge(results_SP,results_RW,by = "node")
results_PPI<-merge(results_PPI,results_D,by = "node")
results_PPI<-merge(clustered_protein,results_PPI,by = "node",all.y =T )

#Combined scores
D_SP<-results_PPI$rank_D*results_PPI$rank_SP
D_RW<-results_PPI$rank_D*results_PPI$rank_RW
SP_RW<-results_PPI$rank_SP*results_PPI$rank_RW
results_PPI<-data.frame(results_PPI,D_SP,D_RW,SP_RW)


####Calculation of priorisation measures in the network with clique####
results_SP_clique<-shortest_path_ranking_function(graph = graph_clique,seed = DE_proteins_network_clique)
results_SP_clique<-data.frame(results_SP_clique,"rank_SP"=rank(-results_SP_clique$SP))

results_RW_clique<-random_walk_ranking_function(graph = graph_clique,seed = DE_proteins_network_clique)
results_RW_clique<-data.frame(results_RW_clique,"rank_RW"=rank(-results_RW_clique$RW))

results_D_clique<-direct_ranking_function(graph = graph_clique,seed = DE_proteins_network_clique)
results_D_clique<-data.frame(results_D_clique,"rank_D"=rank(-results_D_clique$D))

results_D_clique<-results_D_clique[order(results_D_clique$node),]
results_SP_clique<-results_SP_clique[order(results_SP_clique$node),]
results_RW_clique<-results_RW_clique[order(results_RW_clique$node),]

head(results_SP_clique)
head(results_RW_clique)
head(results_D_clique)

#Combined scores
D_SP<-results_D_clique$rank_D*results_SP_clique$rank_SP
D_RW<-results_D_clique$rank_D*results_RW_clique$rank_RW
SP_RW<-results_SP_clique$rank_SP*results_RW_clique$rank_RW

results_D_SP_clique<-data.frame("node"=results_D_clique$node,D_SP)
results_D_RW_clique<-data.frame("node"=results_D_clique$node,D_RW)
results_SP_RW_clique<-data.frame("node"=results_D_clique$node,SP_RW)

####selection of  best nodes for clique for each measure###

###D
select_node<-rep('',nrow(results_D_clique))
for( i in 1:nrow(results_D_clique)){
  node<-as.character(results_D_clique$node[i])
  if(substr(x =node,1,6)=="clique"){
    
    proteins_clique<-clustered_protein[clustered_protein$clique_bigger==substr(x = node,7,1000),"node"]
    results_proteins_clique<-results_D[results_D$node%in%proteins_clique,]
    
    select_node[i]<-as.character(results_proteins_clique$node[which.max(results_proteins_clique$D)[1]])
    
  }
  else{select_node[i]<-as.character(results_D_clique$node[i])}
}
results_D_clique<-cbind(select_node,results_D_clique)
results_D_clique<-results_D_clique[,-which(colnames(results_D_clique)== "node")]


###SP
select_node<-rep('',nrow(results_SP_clique))
for( i in 1:nrow(results_SP_clique)){
  node<-as.character(results_SP_clique$node[i])
  if(substr(x =node,1,6)=="clique"){
    
    proteins_clique<-clustered_protein[clustered_protein$clique_bigger==substr(x = node,7,1000),"node"]
    results_proteins_clique<-results_SP[results_SP$node%in%proteins_clique,]
    
    select_node[i]<-as.character(results_proteins_clique$node[which.max(results_proteins_clique$SP)[1]])
    
  }
  else{select_node[i]<-as.character(results_SP_clique$node[i])}
}
results_SP_clique<-cbind(select_node,results_SP_clique)
results_SP_clique<-results_SP_clique[,-which(colnames(results_SP_clique)== "node")]

###RW
select_node<-rep('',nrow(results_RW_clique))
for( i in 1:nrow(results_RW_clique)){
  node<-as.character(results_RW_clique$node[i])
  if(substr(x =node,1,6)=="clique"){
    
    proteins_clique<-clustered_protein[clustered_protein$clique_bigger==substr(x = node,7,1000),"node"]
    results_proteins_clique<-results_RW[results_RW$node%in%proteins_clique,]
    
    select_node[i]<-as.character(results_proteins_clique$node[which.max(results_proteins_clique$RW)[1]])
    
  }
  else{select_node[i]<-as.character(results_RW_clique$node[i])}
}
results_RW_clique<-cbind(select_node,results_RW_clique)
results_RW_clique<-results_RW_clique[,-which(colnames(results_RW_clique)== "node")]

###D_SP
select_node<-rep('',nrow(results_D_SP_clique))
for( i in 1:nrow(results_D_SP_clique)){
  node<-as.character(results_D_SP_clique$node[i])
  if(substr(x =node,1,6)=="clique"){
    
    proteins_clique<-clustered_protein[clustered_protein$clique_bigger==substr(x = node,7,1000),"node"]
    results_proteins_clique<-results_PPI[results_PPI$node%in%proteins_clique,]
    
    select_node[i]<-as.character(results_proteins_clique$node[which.min(results_proteins_clique$D_SP)[1]])
    
  }
  else{select_node[i]<-as.character(results_D_SP_clique$node[i])}
}
results_D_SP_clique<-cbind(select_node,results_D_SP_clique)
results_D_SP_clique<-results_D_SP_clique[,-which(colnames(results_D_SP_clique)== "node")]

###D_RW
select_node<-rep('',nrow(results_D_RW_clique))
for( i in 1:nrow(results_D_RW_clique)){
  node<-as.character(results_D_RW_clique$node[i])
  if(substr(x =node,1,6)=="clique"){
    
    proteins_clique<-clustered_protein[clustered_protein$clique_bigger==substr(x = node,7,1000),"node"]
    results_proteins_clique<-results_PPI[results_PPI$node%in%proteins_clique,]
    
    select_node[i]<-as.character(results_proteins_clique$node[which.min(results_proteins_clique$D_RW)[1]])
    
  }
  else{select_node[i]<-as.character(results_D_RW_clique$node[i])}
}
results_D_RW_clique<-cbind(select_node,results_D_RW_clique)
results_D_RW_clique<-results_D_RW_clique[,-which(colnames(results_D_RW_clique)== "node")]

###RP_RW
select_node<-rep('',nrow(results_SP_RW_clique))
for( i in 1:nrow(results_SP_RW_clique)){
  node<-as.character(results_SP_RW_clique$node[i])
  if(substr(x =node,1,6)=="clique"){
    
    proteins_clique<-clustered_protein[clustered_protein$clique_bigger==substr(x = node,7,1000),"node"]
    results_proteins_clique<-results_PPI[results_PPI$node%in%proteins_clique,]
    
    select_node[i]<-as.character(results_proteins_clique$node[which.min(results_proteins_clique$SP_RW)[1]])
    
  }
  else{select_node[i]<-as.character(results_SP_RW_clique$node[i])}
}
results_SP_RW_clique<-cbind(select_node,results_SP_RW_clique)
results_SP_RW_clique<-results_SP_RW_clique[,-which(colnames(results_SP_RW_clique)== "node")]

#The main result are the different ranking on the network clique
results<-merge(results_D_clique, results_SP_clique,by="select_node",all=T)
results<-merge(results, results_RW_clique,by="select_node",all=T)
results<-merge(results, results_D_SP_clique,by="select_node",all=T)
results<-merge(results, results_D_RW_clique,by="select_node",all=T)
results<-merge(results, results_SP_RW_clique,by="select_node",all=T)

#The best results has been observed using the Shortest Path ranking combined with the Random walk ranking

####Validation of results####

#CTbase

#I validate this result with CTD base
#'Polycystic Kidney, Autosomal Dominant' is the disease of interest
# It can be replace by you own disease of interest
main_disease <- query_ctd_dise( terms = c( 'Polycystic Kidney, Autosomal Dominant' ) )
disease_genes<-main_disease@gene_interactions@listData$Gene.Symbol
disease_genes_score<-main_disease@gene_interactions@listData$Inference.Score

result_CTDbase<-data.frame(disease_genes,disease_genes_score,stringsAsFactors = F)

####Precision curves
results_D_clique<-results_D_clique[order(results_D_clique$rank_D),]
binary_results_D<-as.numeric(as.character(results_D_clique$select_node)%in%result_CTDbase$disease_genes)

results_SP_clique<-results_SP_clique[order(results_SP_clique$rank_SP),]
binary_results_SP<-as.numeric(as.character(results_SP_clique$select_node)%in%result_CTDbase$disease_genes)

results_RW_clique<-results_RW_clique[order(results_RW_clique$rank_RW),]
binary_results_RW<-as.numeric(as.character(results_RW_clique$select_node)%in%result_CTDbase$disease_genes)

results_D_SP_clique<-results_D_SP_clique[order(results_D_SP_clique$D_SP),]
binary_results_D_SP<-as.numeric(as.character(results_D_SP_clique$select_node)%in%result_CTDbase$disease_genes)

results_D_RW_clique<-results_D_RW_clique[order(results_D_RW_clique$D_RW),]
binary_results_D_RW<-as.numeric(as.character(results_D_RW_clique$select_node)%in%result_CTDbase$disease_genes)

results_SP_RW_clique<-results_SP_RW_clique[order(results_SP_RW_clique$SP_RW),]
binary_results_SP_RW<-as.numeric(as.character(results_SP_RW_clique$select_node)%in%result_CTDbase$disease_genes)

list_binary_results<-list("D"=binary_results_D,"SP"=binary_results_SP,
                          "RW"=binary_results_RW,"D+SP"=binary_results_D_SP,"D+RW"=binary_results_D_RW,
                          "SP+RW"=binary_results_SP_RW)
precision_plot(list_binary_results = list_binary_results,step = 5,size_top = 100,legend = T)

