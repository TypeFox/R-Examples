
##==============================================================================
##  Creates a matrix with the flows
##==============================================================================

Flowmatrix <- function(lim, web=NULL) {

  flowmatrix <- lim$Flowmatrix
  if (is.null(web))
    if (!is.null(lim$Cost) || !is.null(lim$Profit))
      web <- Linp(lim)$X
    else
      web <- Lsei.lim(lim,parsimonious=TRUE)$X

  X          <- as.vector(web)
  Xpos       <- pmax(0.,X)
  ii         <- which(flowmatrix >0,arr.ind=TRUE)
  flowmatrix[ii]<-Xpos[lim$Flowmatrix[ii]]

  Xneg       <- -1*pmin(0.,X)

  if( sum(Xneg) > 0)
    flowmatrix[ii[,c(2,1)]]<-flowmatrix[ii[,c(2,1)]]+Xneg[lim$Flowmatrix[ii]]

  return(flowmatrix)

}

