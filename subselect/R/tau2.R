tau2.coef<-function(mat,H,r,indices,tolval=10*.Machine$double.eps,tolsym=1000*.Machine$double.eps)
{
#   Function tau_2.coef
#   Computes the tau^2 index of "effect magnitude".  
#   This criterion to the minimization of Wilks lambda statistic.
#   Expected input: a variance-covariance (or correlation) matrix, 
#   Effect descrption matrix (H) and its rank (r)

#  error checking

#
#  mat and indices
#

  if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")

  p <- dim(mat)[2]
  validmat(mat,p,tolval,tolsym,allowsingular=TRUE,algorithm="none")

#
# checks on r and H
#

  validnovcrit(mat,criterion="TAU_2",H,r,p,tolval,tolsym)

#
# Computing the criterion value
#

	 tau2.1d<-function(mat,H,r,indices){
		l <- min(r,length(indices))
	    if (length(indices)==1) { detE <- mat[indices,indices]-H[indices,indices]
						detmat <- mat[indices,indices]  }
          else
              { detE <- det(mat[indices,indices]-H[indices,indices])
						detmat <- det(mat[indices,indices])  }

	1 - (detE/detmat)^(1/l)	
       }
      dimension<-length(dim(indices))
      if (dimension > 1){
         tau2.2d<-function(mat,H,r,subsets){
             apply(subsets,1,function(indices){
			
			tau2.1d(mat,H,r,indices)})
            }  
             if (dimension > 2) {
               tau2.3d<-function(mat,H,r,array3d){
                   apply(array3d,3,function(subsets){tau2.2d(mat,H,r,subsets)})
                 }
               output<-tau2.3d(mat,H=H,r,indices)
              }
              if (dimension == 2) {output<-tau2.2d(mat=mat,H,r,indices)}
      }

      if (dimension < 2) {output<-tau2.1d(mat,H,r,indices)}
      output
}

