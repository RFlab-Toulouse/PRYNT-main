source("global.R")

server <- function(input, output) {
  
  
  
  
  
  
  # output$testtabparameters<- reactive({
  INPUT<- reactive({     DP_input<<- input$caption
  print("input")
  DP_input<-strsplit(DP_input, "\n", fixed=FALSE)[[1]]
  DP_not_string<-DP_input[!DP_input%in%annotation$preferred_name]
  DP_not_string
  })
  
    output$value<-renderText({INPUT()})
  
  RESULTS<-eventReactive(input$confirmdatabutton,{
    DP_input<<- input$caption
    DP_input<-strsplit(DP_input, "\n", fixed=FALSE)[[1]]

    ########download the PPI network#####
    ###from string package
    print("network selection")
    string_db <- STRINGdb$new(version="11", species=9606,
                              score_threshold=0, input_directory="" )
    
    annotation<-string_db$get_proteins()
    
    
    string_network<-string_db$get_interactions(string_ids =annotation[,1] )
    string_network_900<-string_network[string_network$combined_score>=900,]
    proteins_900<-unique(c(string_network_900$from,string_network_900$to))
    
    #DP repertoriated in string database
    DP_string<-annotation$protein_external_id[annotation$preferred_name%in%DP_input]
    
    
    ####Contextualization of the network####
    #add relations with DP_string under a score of 900 to the network
    print("contextualisation")
    
    ind1<-string_network$from%in%DP_string &string_network$to%in%proteins_900
    
    ind2<-string_network$from%in%DP_string & string_network$from%in%proteins_900
    
    sup_relation<-string_network[ind1|ind2,]
    # sup_relation<-sup_relation[,c("from","to","score")]
    context_network<-rbind(string_network_900,sup_relation)
    context_network<-unique(context_network)
    proteins_network<-unique(c(context_network$from,context_network$to))
    
    ####Cliques calculation#######
    #I use max_clique from igraph package
    print("calculation of cliques")
  
    graph<-graph_from_data_frame(context_network) # network is transform as an igraph object
    clique<-max_cliques(graph = graph,min=3)
    
    #Transform the list to a matrix
    clique_matrix<-matrix(data = "",nrow = length(clique),ncol =2 )
    clique_matrix[,1]<-as.character(1:length(clique))
    for( i in 1:length(clique)){
      clique_matrix[i,2]<-paste(names(clique[[i]]),sep= "",collapse = " / ")
      
    }
    #clique_matrix is a matrix, each line is a clique, and the second column is the list of proteins is the clique
    
    #I select only the bigger clique for each protein
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
    network_clique<-context_network
    for( i in 1:length(clustered_clique)){
      node_clique<-clustered_protein[clustered_protein$clique_bigger%in%clustered_clique[i],'node']
      network_clique[network_clique$from%in%node_clique,1]<-paste('clique',clustered_clique[i],sep = "")
      network_clique[network_clique$to%in%node_clique,2]<-paste('clique',clustered_clique[i],sep = "")
      
    }
    
    #A lot of relations are now redundant
    #remove redundant relations
    network_clique<-unique(network_clique)
    #remove self-loops
    network_clique<-network_clique[network_clique$from!=network_clique$to,]
    
    
    #network  with clique
    graph_clique<-graph_from_data_frame(network_clique)
    
    protein_network_clique<-names(V(graph_clique))
    number_protein_network_clique<-length(protein_network_clique)
    
    #Selection of the DP proteins in network_clique
    #I consider a clique as DP if at least one protein is DP in the clique
    DP_string_network_clique<-DP_string
    
    for( i in 1:length(DP_string_network_clique)){
      if(DP_string_network_clique[i]%in%clustered_protein$node){
        DP_string_network_clique[i]<-paste('clique',
                                           clustered_protein[clustered_protein$node==DP_string_network_clique[i],2],sep = "")
        
      }}
    
    DP_string_network_clique<-unique(DP_string_network_clique)
    
    #########
    ####Calculation of prioritization measures in the network with clique####
    print("prioritization")
    results_SP_clique<-shortest_path_ranking_function(graph = graph_clique,seed = DP_string_network_clique)
    results_SP_clique<-data.frame(results_SP_clique,"rank_SP"=rank(-results_SP_clique$SP))
    
    results_RW_clique<-random_walk_ranking_function(graph = graph_clique,seed = DP_string_network_clique)
    results_RW_clique<-data.frame(results_RW_clique,"rank_RW"=rank(-results_RW_clique$RW))
    
    results_SP_clique<-results_SP_clique[order(results_SP_clique$node),]
    results_RW_clique<-results_RW_clique[order(results_RW_clique$node),]
    
    head(results_SP_clique)
    head(results_RW_clique)
    
    #Combined scores
    SP_RW<-results_SP_clique$rank_SP*results_RW_clique$rank_RW
    
    results_SP_RW_clique<-data.frame("node"=results_SP_clique$node,SP_RW)
    
    ####selection of  best nodes for clique for each measure###
    
    ####Calculation of ranking in the network####
    results_SP<-shortest_path_ranking_function(graph = graph,seed = DP_string)
    results_SP<-data.frame(results_SP,"rank_SP"=rank(-results_SP$SP))
    
    results_RW<-random_walk_ranking_function(graph = graph,seed = DP_string)
    results_RW<-data.frame(results_RW,"rank_RW"=rank(-results_RW$RW))
    
    head(results_SP)
    head(results_RW)
    
    #Combined scores
    SP_RW<-results_SP$rank_SP*results_RW$rank_RW
    results_SP_RW<-data.frame("node"=results_SP$node,SP_RW)
    
    ####RP_RW####
    select_node<-rep('',nrow(results_SP_RW_clique))
    for( i in 1:nrow(results_SP_RW_clique)){
      node<-as.character(results_SP_RW_clique$node[i])
      if(substr(x =node,1,6)=="clique"){
        
        proteins_clique<-clustered_protein[clustered_protein$clique_bigger==substr(x = node,7,1000),"node"]
        results_proteins_clique<-results_SP_RW[results_SP_RW$node%in%proteins_clique,]
        
        select_node[i]<-as.character(results_proteins_clique$node[which.min(results_proteins_clique$SP_RW)[1]])
        
      }
      else{select_node[i]<-as.character(results_SP_RW_clique$node[i])}
    }
    results_SP_RW_clique<-cbind(select_node,results_SP_RW_clique)
    results_SP_RW_clique<-results_SP_RW_clique[,-which(colnames(results_SP_RW_clique)== "node")]
    
    ####The main result are the different ranking on the network clique#####
    print("done")
    

    final<-merge(x = annotation[,1:2],y= results_SP_RW_clique,by.x="protein_external_id",by.y="select_node")
    final<-final[order(final$SP_RW),2:3]
    rownames(final)<-1:nrow(final)
    final<-data.frame("rank"=1:nrow(final),final)
    colnames(final)<-c("rank","Gene symbol","score")
    print(final[1:10,])
    
    list("PRYNT"=final)
   
    
  })
  
  output$table =renderDataTable( {RESULTS()$PRYNT }, options = list(    "orderClasses" = F,
                                                                 "responsive" = F,
                                                                 "pageLength" = 10,
                                                                 "rownames"=TRUE,
                                                                 "colnames"=TRUE))
  
  output$downloadtable <- downloadHandler(
    filename = function() { paste('dataset', '.','csv', sep='') },
    content = function(file) {
      downloaddataset(   RESULTS()$PRYNT, file) })
  
  
}






