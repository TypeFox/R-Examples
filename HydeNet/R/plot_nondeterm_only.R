#* Utility functions for plotting Hyde Networks without the deterministic nodes

#* plot_nondeterm_only is the main function in this process.  It systematically
#* searches each node for deterministic parents.  When a deterministic parent
#* is found, it scans that node's parents, and so on until there are no more
#* deterministic parents in the chain.
#* 
#* It then uses the information about parents to fit a new HydeNetwork, 
#* pulls the attributes relevant to plotting from the original, 
#* and returns the modified object for plotting.

plot_nondeterm_only <- function(network){
  non_determ_nodes <- network$nodes[!vapply(network$nodeType,
                                            function(x) x == "determ",
                                            logical(1))]
  non_determ_nodes <- non_determ_nodes[!is.na(non_determ_nodes)]
  
  form_parts <- 
    vapply(non_determ_nodes, 
           non_determ_parents_form,
           character(1),
           network)
  form <- as.formula(paste0(" ~ ", paste0(form_parts, collapse = " + ")))
  newNet <- HydeNetwork(form)
  newNet$nodeType <- network$nodeType[newNet$nodes]
  newNet$nodeUtility <- network$nodeUtility[newNet$nodes]
  newNet$nodeDecision <- network$nodeDecision[newNet$nodes]
  newNet
}

#* Creates a portion of the HydeNet formula pertaining to a specific node
#* The while loop searches up the parent-chain until none of the 
#* ancestry shows up as deterministic.

non_determ_parents_form <- function(node, network){
  parent_types <- parent_type(node, network)
  
  while(any(determ_parent(parent_types))){
    have_determ <- names(parent_types)[determ_parent(parent_types)]
    
    super_parents <- do.call("c", lapply(names(parent_types[have_determ]),
                                         parent_type,
                                         network))
    
    parent_types <- c(parent_types[!names(parent_types) %in% have_determ], super_parents)
  }
  
  form <- paste0(node, " | ", paste0(unique(names(parent_types)), collapse = " * "))
  gsub(" [|] $", "", form)
}

#* Utility function to extract the nodeType of the parent
parent_type <- function(node, network){
  network$nodeType[network$parents[[node]]]
}

#* Return a logical vector for the parent nodes indicating if they are deterministic
determ_parent <- function(parent_types){
  vapply(parent_types, 
             function(pt) "determ" %in% pt, 
             logical(1))
}

