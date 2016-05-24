gcd.coef<-function(mat, indices, pcindices = NULL)
{
#   
#       calcula o GCD entre um subconjunto ("indices") de variaveis
#       e um subconjunto ("pcindices") das CPs de todas as variaveis, 
#       cuja matriz de covariancias e "mat".

#  error checking

         if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")
         if  (!is.null(pcindices) & (sum(!(as.integer(pcindices) == pcindices)) > 0)) stop("\n The PC indices must be integers")

         if (!is.matrix(mat)) {
         stop("Data is missing or is not given in matrix form")}
     if (dim(mat)[1] != dim(mat)[2]) {
         mat<-cor(mat)
         warning("Data must be given as a covariance or correlation matrix. \n It has been assumed that you wanted the correlation matrix of the \n data matrix which was supplied.")
     }        

# body of function

# initializations

      if (!is.null(pcindices)) {
       if (!is.vector(pcindices)) stop("If Principal Components are user-specified, only one set of PCs is allowed for each function call")
      }
  
     dvsmat <- svd(mat)
      tr<-function(mat){sum(diag(mat))}
      gcd.1d<-function(mat,indices, pcindices){
             if (is.null(pcindices)) {pcindices <- 1:sum(!indices==0)}
             indices<-indices[!indices == 0]
             svdapprox <- function(mat, indices) {
             t(dvsmat$v[, indices] %*% (t(dvsmat$u[, indices]) * dvsmat$d[indices]))
                           }
             sum(diag(solve(mat[indices, indices]) %*% svdapprox(mat, 
             pcindices)[indices, indices]))/sqrt(length(indices) * 
             length(pcindices))
        }
      dimension<-length(dim(indices))

# output for each dimension of input array

      if (dimension > 1){
         gcd.2d<-function(mat,subsets,pcindices){
             apply(subsets,1,function(indices){gcd.1d(mat,indices,pcindices)})
            }  
           if (dimension > 2) {               
            gcd.3d<-function(mat,array3d,pcindices){
             apply(array3d,3,function(subsets){gcd.2d(mat,subsets,pcindices)})
            }
            output<-gcd.3d(mat,indices,pcindices)
           }
           if (dimension == 2) {output<-gcd.2d(mat,indices,pcindices)}
      }

      if (dimension < 2) {output<-gcd.1d(mat,indices,pcindices)}
      output
}
