identify.cna <- function(x, labels=NULL, cna=NULL, ...){

  ## Be carefull with input argument order
  ##  - 'labels' can take any input and screw up priniting
  ##      e.g. if you pass cna as the second argument!
  ## Should this perhaps be able to take just a cna object as input
  ##   - Possible if cna object has layout defined
  ##   - Could take extra layout option for custom graphs
  ## x <- plot(net)
  ## 
  ## d <- identify.cna(x, cna=net)
  ## d <- identify.cna(x, labels=summary(net)$members)
 

  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
    stop("igraph package missing: Please install, see: ?install.packages")
  }

  if(dim(x)[2] != 2){
    stop("'x' object must be a Nx2 numeric matrix")
  }
  
  x.norm <- igraph::layout.norm(x, -1, 1, -1, 1)

  if( !is.null(labels) ) {
    ## Use input labels
    inds <- identify(x.norm[,1], x.norm[,2], labels, ...)
    return( labels[inds] ) 
  } else {
    if(is.null(cna)) {
      ## Use standard labels
      inds <- identify(x.norm[,1], x.norm[,2], ...)
      return(inds)
    } else {
      ## Take labels from cna object!!
      labels.all <- summary.cna(cna)
      labels.short <- labels.all$tbl$members
      labels.full <- labels.all$members
      inds <- identify(x.norm[,1], x.norm[,2], labels.short, ...)
      return( labels.full[inds] )
    }
  }
}
