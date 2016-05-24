"getAniSuper" <- function (A, ani) 
{
	A_p <- A
	if(length(ani$super) > 0) {
	  for(i in 1:length(ani$super)) {
	      A_p <- cbind(A_p, rowSums( A[, ani$super[[i]]])) 
	  }
        }
	A_p
}
