zeta2.coef<-function(mat,H,r,indices,tolval=10*.Machine$double.eps,tolsym=1000*.Machine$double.eps)
{
#   Function zeta_2.coef
#   Computes the zeta^2 index of "effect magnitude".  
#   This criterion is equivalent to the maximization of Hotteling-Lawley trace statistic.
#   Expected input: a variance-covariance (or correlation) matrix, 
#   Effect descrption matrix (H) and its rank (r)

#  error checking

#  mat and indices

  if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")


  p <- dim(mat)[2]
  validmat(mat,p,tolval,tolsym,allowsingular=TRUE,algorithm="none")


# checks on r and H

  validnovcrit(mat,criterion="ZETA_2",H,r,p,tolval,tolsym)

# Computing the criterion value


      tr<-function(mat){sum(diag(mat))}

      zeta2.1d<-function(mat,H,r,indices){
		l <- min(r,length(indices))
	V <-tr(solve(mat[indices,indices]-H[indices,indices],H[indices,indices]))
	V/(V+l)	
       }
      dimension<-length(dim(indices))
      if (dimension > 1){
         zeta2.2d<-function(mat,H,r,subsets){
             apply(subsets,1,function(indices){
			
			zeta2.1d(mat,H,r,indices)})
            }  
             if (dimension > 2) {
               zeta2.3d<-function(mat,H,r,array3d){
                   apply(array3d,3,function(subsets){zeta2.2d(mat,H,r,subsets)})
                 }
               output<-zeta2.3d(mat,H=H,r,indices)
              }
              if (dimension == 2) {output<-zeta2.2d(mat=mat,H,r,indices)}
      }

      if (dimension < 2) {output<-zeta2.1d(mat,H,r,indices)}
      output
}

