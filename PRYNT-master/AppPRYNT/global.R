####Packages####
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

#usePackage("igraph")
usePackage("zoo")
usePackage("DT")

usePackage("shinycssloaders")

usePackage_bioconductor <- function(p) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!is.element(p, installed.packages()[,1])) {
    BiocManager::install(pkgs = p, update = TRUE, force = TRUE)
  }
  # if (!is.element(p, installed.packages()[,1]) & p %in% c("Rgraphviz", "STRINGdb")) {
  #   BiocManager::install(pkgs = p, update = TRUE, force = TRUE)
  # }
  require(p, character.only = TRUE)
}


# Installation des packages requis
packages_to_install <- function() {
  # Configuration des options de mémoire pour R
  options(future.globals.maxSize = 4000 * 1024^2)  # Set to 4GB
  
  # Packages de base
  usePackage("devtools")
  usePackage("shiny")
  
  # Packages Bioconductor
  usePackage_bioconductor("STRINGdb")
  usePackage_bioconductor("Rgraphviz")
  
  # Packages GitHub avec gestion des timeouts
  options(timeout = 300)  # Augmente le timeout à 5 minutes
  
  # Packages GitHub
  if (!require("dnet")) {
    devtools::install_github("hfang-bristol/dnet", upgrade = "never")
  }
  
  if (!require("RandomWalkRestartMH")) {
    devtools::install_github("alberto-valdeolivas/RandomWalkRestartMH", upgrade = "never")
  }
}

packages_to_install()
library(RandomWalkRestartMH)
# usePackage("RandomWalkRestartMH")
# usePackage_bioconductor("RandomWalkRestartMH")
# usePackage_bioconductor("STRINGdb") 
#######Functions#######

#memory.limit(size =8078 )
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

random_walk_ranking_function <- function(graph, seed, adjacency_matrix=NULL, multiplex_object=NULL) {
  # La correction principale est ici
  if(is.null(multiplex_object)) {
    # Créer une liste contenant le graphe
    graph_list <- list(PPI = graph)
    # Créer le multiplex avec la liste
    multiplex_object <- create.multiplex(graph_list, Layers_Name=c("PPI"))
  }
  
  if(is.null(adjacency_matrix)) {
    adjacency_matrix <- compute.adjacency.matrix(multiplex_object)
  }
  
  adjacency_matrix_norm <- normalize.multiplex.adjacency(adjacency_matrix)
  seed <- seed[seed %in% names(V(graph))]
  
  RWR_PPI_Results <- Random.Walk.Restart.Multiplex(adjacency_matrix_norm,
                                                   multiplex_object, 
                                                   seed)
  
  results_random_walk <- RWR_PPI_Results$RWRM_Results
  colnames(results_random_walk) <- c("node", "RW")
  node <- data.frame("node" = names(V(graph)))
  results_random_walk <- merge(x = node, 
                               y = results_random_walk, 
                               by = "node", 
                               all.x = TRUE)
  results_random_walk$RW[is.na(results_random_walk$RW)] <- 0
  results_random_walk <- results_random_walk[order(results_random_walk$node), ]
  
  return(results_random_walk)
}

# random_walk_ranking_function<-function(graph,seed,adjacency_matrix=NULL,multiplex_object=NULL){
#   #calculate the closeness of each protein from seeds using the random walk 
#   # using the RandomWalkRestartMH package
#   #graph : igraph object, containing the network to analyse
#   #seed : character vector of DE proteins
#   
#   #adjacency_matrix : to have the possibility to calculate the adjency matrix outside of the function
#   #multiplex_object : to have the possibility to calculate the multiplex object outside of the function
#   
#   cat("calss of graph  : ", class(graph))
#   if(is.null(adjacency_matrix)){
#     multiplex_object <- create.multiplex(graph,Layers_Name=c("PPI"))
#   }
#   if(is.null(adjacency_matrix)){
#     adjacency_matrix <- compute.adjacency.matrix(multiplex_object)
#   }
#   adjacency_matrix_norm <- normalize.multiplex.adjacency(adjacency_matrix)
#   seed<-seed[seed%in%names(V(graph))]
#   ## We launch the algorithm with the default parameters (See details on manual)
#   RWR_PPI_Results <- Random.Walk.Restart.Multiplex(adjacency_matrix_norm,
#                                                    multiplex_object,seed)
#   # We display the results
#   results_random_walk<-RWR_PPI_Results$RWRM_Results
#   colnames(results_random_walk)<-c("node","RW")
#   node<-data.frame("node"=names(V(graph)))
#   results_random_walk<-merge(x = node,y =results_random_walk,by="node",all.x = T )
#   results_random_walk$RW[is.na(results_random_walk$RW)]<-0
#   results_random_walk<-results_random_walk[order(results_random_walk$node),]
#   
#   return(results_random_walk)
# }

downloaddataset<-function(x,file,cnames=T,rnames=T){
  ext<-tools::file_ext(file) #strsplit(x = file,split = "[.]")[[1]][2]
  if(ext=="csv"){
    if(sum(cnames,rnames)==2){write.csv(x,file)}
    else{write.table(x,file,col.names = cnames,row.names = F,sep=";",dec=".")}
  }
  if(ext=="xlsx"){
    write.xlsx(x,file,col.names = cnames,row.names =rnames )
    #writexl::write_xlsx(x,file, col_names = cnames)
  }
  
}



# Supprimer ces lignes à la fin du fichier global.R car elles causent l'erreur
# string_network<<-string_db$get_interactions(string_ids =annotation[,1] )
# string_network_900<<-string_network[string_network$combined_score>=900,]
# proteins_900<<-unique(c(string_network_900$from,string_network_900$to))

########download the PPI network#####
###from string package
print("loading String Database")
initializer <- function() {
  string_db <- STRINGdb$new(version="11", species=9606,
                            score_threshold=0, input_directory="")
  annotation <- string_db$get_proteins()
  string_network <- string_db$get_interactions(string_ids = annotation[,1])
  string_network_900 <- string_network[string_network$combined_score>=900,]
  proteins_900 <- unique(c(string_network_900$from, string_network_900$to))
  
  return(list(
    string_db = string_db,
    annotation = annotation,
    string_network = string_network,
    string_network_900 = string_network_900,
    proteins_900 = proteins_900
  ))
}

# Initialisation des données au démarrage
data <- initializer()

# Rendre les données disponibles globalement si nécessaire
string_db <- data$string_db
annotation <- data$annotation
string_network <- data$string_network
string_network_900 <- data$string_network_900
proteins_900 <- data$proteins_900