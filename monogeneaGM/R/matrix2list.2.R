#' Partitioning a matrix by row labels into objects of a list
#'
#' This function partitions a matrix according to row labels and assigns the partitioned
#' submatrices as objects of a list.
#' @param x a matrix with row labels, typically species names 
#' @return a list of objects; the number of objects is equal to the length of levels of the rownames of \code{x}.
#' Each object is a matrix with the same row names.
#' @details The output from this function is passed as input for \code{boxplotSort}.
#' @seealso \code{\link{boxplotSort}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' data(pwed_pd)
#' data(spcolmap)
#'
#' pwed_pd <- matrix2list.2(pwed_pd)
#'
#' cladeI <- spcolmap$species[spcolmap$host %in% "M.buchanani"]
#' #We just want to look at distance between LM1 and LM3 in for dorsal anchor
#' boxplotSort(lapply(pwed_pd, function(k) k[,which(colnames(k)=="D1_3")]), 
#' italic=TRUE, col=c("dodgerblue","violetred"), clade=cladeI,
#' ylab=expression(paste("Length ", "(", italic(mu),"m", ")")))
#'
#' #Separation of two lineage seems possible at 15 micrometers
#' abline(h=15)
#'

matrix2list.2 <- function(x) {

splabel <- levels(as.factor(rownames(x)))
templist <- vector("list",length(splabel))
for(k in 1:length(templist)){
templist[[k]] <- x[which(substr(rownames(x),1,20) == splabel[k]),]
}

names(templist) <- splabel
return(templist)

}



