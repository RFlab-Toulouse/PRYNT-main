####Packages####
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

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
  require(p, character.only = TRUE)
}

# Configuration des options de mémoire et de calcul
options(future.globals.maxSize = 4000 * 1024^2)  # 4GB
options(timeout = 300)  # 5 minutes timeout

# Installation des packages nécessaires
packages_to_install <- function() {
  # Packages base et Bioconductor
  usePackage("devtools")
  usePackage("shiny")
  usePackage_bioconductor("STRINGdb")
  usePackage_bioconductor("Rgraphviz")
  
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

#######Functions#######
shortest_path_ranking_function <- function(graph, seed, score=NULL, directed=TRUE, distance_table=NULL) {
  if(is.null(distance_table)) {
    sens <- if(directed) "out" else "all"
    distance_table <- distances(graph, v = V(graph), to = V(graph), mode = sens, weights = score)
  }
  
  distance_table_to_seed <- distance_table[, colnames(distance_table) %in% seed]
  distance_table_to_seed[which(distance_table_to_seed == 0, arr.ind = TRUE)] <- Inf
  
  results_SP <- apply(X = distance_table_to_seed, MARGIN = 1, function(x) {sum(1/x)})
  results_SP <- data.frame("node" = names(results_SP), "SP" = results_SP)
  results_SP <- results_SP[order(results_SP$node), ]
  return(results_SP)
}

random_walk_ranking_function <- function(graph, seed, adjacency_matrix=NULL, multiplex_object=NULL) {
  if(is.null(multiplex_object)) {
    graph_list <- list(PPI = graph)
    multiplex_object <- create.multiplex(graph_list, Layers_Name=c("PPI"))
  }
  
  if(is.null(adjacency_matrix)) {
    adjacency_matrix <- compute.adjacency.matrix(multiplex_object)
  }
  
  # Utilisation de Matrix explicitement pour la normalisation
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

downloaddataset <- function(x, file, cnames=TRUE, rnames=TRUE) {
  ext <- tools::file_ext(file)
  if(ext == "csv") {
    if(sum(cnames, rnames) == 2) {
      write.csv(x, file)
    } else {
      write.table(x, file, col.names = cnames, row.names = FALSE, sep = ";", dec = ".")
    }
  }
  if(ext == "xlsx") {
    write.xlsx(x, file, col.names = cnames, row.names = rnames)
  }
}

# Initialisation de la base de données STRING
print("loading String Database")
data <- local({
  string_db <- STRINGdb$new(version="11", species=9606,
                            score_threshold=0, input_directory="")
  annotation <- string_db$get_proteins()
  string_network <- string_db$get_interactions(string_ids = annotation[,1])
  string_network_900 <- string_network[string_network$combined_score >= 900,]
  proteins_900 <- unique(c(string_network_900$from, string_network_900$to))
  
  list(
    string_db = string_db,
    annotation = annotation,
    string_network = string_network,
    string_network_900 = string_network_900,
    proteins_900 = proteins_900
  )
})