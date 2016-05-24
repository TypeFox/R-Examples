descendants <- function(phy, node, type = "t", ignore.tip = TRUE, labels = FALSE){
	
  # checks and definitions
  # ----------------------
  if ( !inherits(phy, "phylo") ) stop ("'phy' is not of class phylo")
  type <- match.arg(type, c("both", "internal", "terminal"))
	tips <- seq_along(phy$tip.label)
  
  # 'node' is a tip 
  # ---------------
  if ( node <= max(tips) ){
    if ( ignore.tip ) x <- node
    else stop("node ", node, " is not an internal node") 
  }
  # normal procedure when 'node' is internal
  # ----------------------------------------
  else {
    x <- phy$edge[phy$edge[,1] == node, 2]
    repeat{
      xx <- x
      x <- sort(unique(c(x, phy$edge[,2][phy$edge[,1] %in% x])))
      if (identical(x, xx)) break
    }
    
    # apply 'type' argument:
    # -----------------------------------------
    if (type == "terminal") {
      x <- intersect(x, tips)
      if (labels) {
        x <- phy$tip.label[x]
      }
    }
    if (type == "internal") x <- setdiff(x, tips)
  }
	return(x)
}