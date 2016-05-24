rm.coef<-function(mat, indices)
{

#   Computes the matrix correlation between data matrices and their 
#   regression on a subset of their variables. Expected input is a
#   variance-covariance (or correlation) matrix. 

#  error checking


  if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")
  if (!is.matrix(mat)) {
         stop("Data is missing or is not given in matrix form")}
     if (dim(mat)[1] != dim(mat)[2]) {
         mat<-cor(mat)
         warning("Data must be given as a covariance or correlation matrix. \n It has been assumed that you wanted the correlation matrix of the \n data matrix which was supplied.")
       }
      tr<-function(mat){sum(diag(mat))}
      rm.1d<-function(mat,indices){
        sqrt(tr((mat %*% mat)[indices,indices] %*% solve(mat[indices,indices]))/tr(mat))
        }
      dimension<-length(dim(indices))
      if (dimension > 1){
         rm.2d<-function(mat,subsets){
             apply(subsets,1,function(indices){rm.1d(mat,indices)})
            }  
             if (dimension > 2) {
               rm.3d<-function(mat,array3d){
                   apply(array3d,3,function(subsets){rm.2d(mat,subsets)})
                 }
               output<-rm.3d(mat,indices)
              }
              if (dimension == 2) {output<-rm.2d(mat,indices)}
      }

      if (dimension < 2) {output<-rm.1d(mat,indices)}
      output
}
