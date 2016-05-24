rv.coef<-function(mat, indices)
{
#    Computes Escoufier's RV-coefficient for the configuration of
#    points defined by n observations of a set of p variables, and by
#    the regression of all variables on a subset of k variables
#    (given by \code{indices}). 

#  error checking

     if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")
     if (!is.matrix(mat)) {
         stop("Data is missing or is not given in matrix form")}
     if (dim(mat)[1] != dim(mat)[2]) {
         mat<-cor(mat)
         warning("Data must be given as a covariance or correlation matrix. \n It has been assumed that you wanted the correlation matrix of the \n data matrix which was supplied.")
       }
      tr<-function(mat){sum(diag(mat))}
      rv.1d<-function(mat,indices){
             mat2 <- (mat %*% mat)[indices, indices]
             invmatk <- solve(mat[indices, indices])
             sqrt(tr(mat2 %*% invmatk %*% mat2 %*% invmatk)/tr(mat %*% mat))
        }
      dimension<-length(dim(indices))
      if (dimension > 1){
         rv.2d<-function(mat,subsets){
             apply(subsets,1,function(indices){rv.1d(mat,indices)})
            }  
             if (dimension > 2) {
               rv.3d<-function(mat,array3d){
                   apply(array3d,3,function(subsets){rv.2d(mat,subsets)})
                 }
               output<-rv.3d(mat,indices)
              }
              if (dimension == 2) {output<-rv.2d(mat,indices)}
      }

      if (dimension < 2) {output<-rv.1d(mat,indices)}
      output
}
