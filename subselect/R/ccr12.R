ccr12.coef<-function(mat,H,r,indices,tolval=10*.Machine$double.eps,tolsym=1000*.Machine$double.eps)
{
#   Function ccr1_2.coef
#   Computes the first squared canonical correlation.  
#   This criterion is equivalent to the maximization of Roy first root.
#   Expected input: a variance-covariance (or correlation) matrix, 
#   Effect descrption matrix (H) and its rank (r)

#  error checking

#  mat and indices

  if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")

  p <- dim(mat)[2]
  validmat(mat,p,tolval,tolsym,allowsingular=TRUE,algorithm="none")


# checks on r and H

  validnovcrit(mat,criterion="CCR1_2",H,r,p,tolval,tolsym)
  
#  Computing the criterion value  

      ccr12.1d<-function(mat,H,r,indices){
		
	Re(eigen(solve(mat[indices,indices],H[indices,indices]))$values[1])
		
       }
      dimension<-length(dim(indices))
      if (dimension > 1){
         ccr12.2d<-function(mat,H,r,subsets){
             apply(subsets,1,function(indices){
			
			ccr12.1d(mat,H,r,indices)})
            }  
             if (dimension > 2) {
               ccr12.3d<-function(mat,H,r,array3d){
                   apply(array3d,3,function(subsets){ccr12.2d(mat,H,r,subsets)})
                 }
               output<-ccr12.3d(mat,H=H,r,indices)
              }
              if (dimension == 2) {output<-ccr12.2d(mat=mat,H,r,indices)}
      }

      if (dimension < 2) {output<-ccr12.1d(mat,H,r,indices)}
      output
}

