sister <- function(phy, node, type = "terminal"){
  
  # checks and definitions
  # ----------------------
  if ( !inherits(phy, "phylo") ) 
    stop ("'phy' is not of class 'phylo'")
  tips <- 1:Ntip(phy)
  if ( is.character(node) )
    node <- which(phy$tip.label %in% node)
  
  if ( node == (Ntip(phy) + 1) )
    stop("node = ", node, " is root node")
  else {
    obj <- phy$edge[, 1][phy$edge[, 2] == node] # getmrca
    obj <- descendants(phy, obj, type = type) # get whole sister clade
    obj <- setdiff(obj, node) # eliminate node 
  }
  obj
}
