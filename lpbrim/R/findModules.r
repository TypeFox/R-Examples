#' @title Find modules
#' @description This function takes a matrix, and performs a search for the best
#' partition in modules.
# Example with parallel
# M <- matrix(rbinom(10000, 1, 0.01), ncol=100)
# M <- M[rowSums(M)>0, colSums(M)>0]
# registerDoMC(3)
# mod <- findModules(M, iter=10, sparse=TRUE, .parallel=TRUE)
#' @export
#' @param M An adjacency matrix
#' @param iter Number of optimization runs to do
#' @param sparse Whether the matrix should be made sparse
#' @param ... Other arguments
#' @examples
#' M <- matrix(rbinom(100, 1, 0.3), ncol=10)
#' M <- M[rowSums(M)>0, colSums(M)>0]
#' mod <- findModules(M, iter=2, sparse=FALSE)
findModules = function(M,iter=50,sparse=TRUE, ...)
{
   if(is.null(rownames(M))) rownames(M) <- paste('r',c(1:NROW(M)),sep='')
   if(is.null(colnames(M))) colnames(M) <- paste('c',c(1:NCOL(M)),sep='')
   if(sparse) M <- Matrix::Matrix(M, sparse=TRUE)
   ModulOutput <- plyr::alply(c(1:iter),1, function(x) bBRIM(M), ...)
   Qs <- unlist(plyr::laply(ModulOutput, function(x)x$Q))
   maxQs <- which.max(Qs)
   return(ModulOutput[[maxQs]])
}
