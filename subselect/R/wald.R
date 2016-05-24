wald.coef<-function(mat,H,indices, tolval=10*.Machine$double.eps, tolsym=1000*.Machine$double.eps)
{
#   Function wald.coef
#   Computes the value of Wald statistic for testing the significance of the omited variables 
#   in a generalized linear model 

#  error checking

#  mat and indices

  if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")

  p <- dim(mat)[2]
  validmat(mat,p,tolval,tolsym)
    


# checks on r and H

  validnovcrit(mat,criterion="WALD",H,r=1,p,tolval,tolsym)


# Computing the criterion value

      tr <- function(mat){sum(diag(mat))}
      waldallvar <- tr(solve(mat,H))

      wald.1d <- function(mat,H,indices){ waldallvar - tr(solve(mat[indices,indices],H[indices,indices])) }
      dimension<-length(dim(indices))
      if (dimension > 1){
         wald.2d<-function(mat,H,subsets){
             apply(subsets,1,function(indices){
			
			wald.1d(mat,H,indices)})
            }  
             if (dimension > 2) {
               wald.3d<-function(mat,H,array3d){
                   apply(array3d,3,function(subsets){wald.2d(mat,H,subsets)})
                 }
               output<-wald.3d(mat,H=H,indices)
              }
              if (dimension == 2) {output<-wald.2d(mat=mat,H,indices)}
      }

      if (dimension < 2) {output<-wald.1d(mat,H,indices)}
      output

}

