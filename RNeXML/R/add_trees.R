#' add_trees
#' 
#' add_trees
#' @param phy a phylo object, multiPhylo object, or list of 
#'  mulitPhylo to be added to the nexml
#' @param nexml a nexml object to which we should append this phylo.
#'  By default, a new nexml object will be created.  
#' @param append_to_existing_otus logical, indicating if we should 
#'  make a new OTU block (default) or append to the existing one. 
#' @return a nexml object containing the phy in nexml format. 
#' @export 
#' @examples 
#' library("geiger")
#' data(geospiza)
#' geiger_nex <- add_trees(geospiza$phy)
add_trees <- function(phy, 
                      nexml=new("nexml"), 
                      append_to_existing_otus=FALSE){
  nexml <- as(nexml, "nexml")

  phy <- standardize_phylo_list(phy)
  ## handle multiPhlyo cases
  new_taxa <- unlist(sapply(phy, function(y)
                      sapply(y, function(z) 
                        z$tip.label)))

  nexml <- add_otu(nexml, new_taxa, append=append_to_existing_otus)
  otus_id <- nexml@otus[[length(nexml@otus)]]@id
  nexml <- add_trees_block(nexml, phy, otus_id)
  nexml
}


##################### phylo -> nexml ###############

setAs("phylo", "nexml", function(from){
      add_trees(from)
})

setAs("multiPhylo", "nexml", function(from){
      add_trees(from)
})


standardize_phylo_list <- function(phy){
  if(is(phy, "list") && 
     (is(phy[[1]], "list") || is(phy[[1]], "multiPhylo")) && 
     is(phy[[1]][[1]], "phylo")){
    phy 
  } else if(is(phy, "multiPhylo") || (is(phy, "list") && is(phy[[1]], "phylo"))) { 
    list(phy)
  } else if(is(phy, "phylo")) {
    phy <- list(phy)
    class(phy) <- "multiPhylo"
    list(phy) 
  } else {
    # desperate
    phy <- list(as(phy, "phylo"))
    class(phy) <- "multiPhylo"
    list(phy) 
  }
}

add_trees_block <- function(nexml, phy, otus_id){
  phy <- standardize_phylo_list(phy)
  ## all trees will use the same  
  otu_map <- reverse_map(get_otu_maps(nexml))[[otus_id]]

  trees <- lapply(phy, function(trs){
         tree_id <- nexml_id("ts")
         new("trees", 
             id = tree_id,
             about = paste0("#", tree_id),
             otus = otus_id,
             tree = new("ListOftree", 
                        lapply(trs, function(tr)
                           fromPhylo(tr, otu_map)))
            )
  })

  ## Append to any existing trees nodes 
  nexml@trees <- new("ListOftrees", c(nexml@trees, trees))
  nexml
}


# Main routine to generate NeXML from ape:phylo
fromPhylo <- function(phy, otu_map){
  
  node_ids <- sapply(unique(as.numeric(phy$edge)), 
                     function(i) nexml_id("n")) 
  names(node_ids) <- as.character(unique(as.numeric(phy$edge)))

  ## Generate the "ListOfedge" made of "edge" objects
  edges <- 
    lapply(1:dim(phy$edge)[1], 
           function(i){
            edge_id <- nexml_id("e")
            source <- node_ids[as.character(phy$edge[i,1])]
            target <- node_ids[as.character(phy$edge[i,2])]
            e <- new("edge", 
                     source = source, 
                     target = target, 
                     id = edge_id,
                     about = paste0("#", edge_id))
           if(!is.null(phy$edge.length))
             e@length <- as.numeric(phy$edge.length[i])
           e
           }
  )
  edges <- new("ListOfedge", edges)
  ## Generate the ListOfnode made of "node" objects
  ## In doing so, generate otu_id numbers for tip nodes
  nodes <- lapply(unique(as.numeric(phy$edge)), function(i){
    node_id <- node_ids[as.character(i)] 
    if(is.na(phy$tip.label[i]))
      new("node", id = node_id, about = paste0("#", node_id))
    else if(is.character(phy$tip.label[i])){
      otu_id <- otu_map[phy$tip.label[i]]
      new("node", 
          id = node_id, 
          about = paste0("#", node_id), 
          otu = otu_id)  
    }
  })
  ## FIXME how about naming non-tip labels?  
  nodes <- new("ListOfnode", nodes)

 

  ## Create the "tree" S4 object
  tree_id <- nexml_id("tree") 
  tree <- new("tree", 
      node = nodes, 
      edge = edges,
      'xsi:type' = 'FloatTree',
      id = tree_id,
      about = paste0("#", tree_id))  
}



