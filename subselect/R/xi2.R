xi2.coef<-function(mat,H,r,indices, tolval=10*.Machine$double.eps, tolsym=1000*.Machine$double.eps)
{
#   Function xi_2.coef
#   Computes the xi^2 index of "effect magnitude".  
#   This criterion is equivalent to the maximization of Bartlett-Pillai trace statistic.
#   Expected input: a variance-covariance (or correlation) matrix, 
#   Effect descrption matrix (H) and its rank (r)

#  error checking

#  mat and indices

  if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")

  p <- dim(mat)[2]
  validmat(mat,p,tolval,tolsym,allowsingular=TRUE,algorithm="none")
    


# checks on r and H

  validnovcrit(mat,criterion="XI_2",H,r,p,tolval,tolsym)


# Computing the criterion value

      tr<-function(mat){sum(diag(mat))}

      xi2.1d<-function(mat,H,r,indices){
		l <- min(r,length(indices))
	tr(solve(mat[indices,indices],H[indices,indices]))/l	
       }
      dimension<-length(dim(indices))
      if (dimension > 1){
         xi2.2d<-function(mat,H,r,subsets){
             apply(subsets,1,function(indices){
			
			xi2.1d(mat,H,r,indices)})
            }  
             if (dimension > 2) {
               xi2.3d<-function(mat,H,r,array3d){
                   apply(array3d,3,function(subsets){xi2.2d(mat,H,r,subsets)})
                 }
               output<-xi2.3d(mat,H=H,r,indices)
              }
              if (dimension == 2) {output<-xi2.2d(mat=mat,H,r,indices)}
      }

      if (dimension < 2) {output<-xi2.1d(mat,H,r,indices)}
      output
}

