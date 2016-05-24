#Perform Cholesky decomposition and do some adjustments if the decomposition fails.
#Described in Appendeix:NORTA RVG Methods step 3.2 according to the reference paper.
Perform_chol <- function(r_matrix, w_k_bar,
                         invcdfnames, paramslists) {
   ndim <- ncol(r_matrix)
   mk <- nrow(w_k_bar)
   tmp_paramslists <- list()
   #For better Cholesky decomposition, and  the while statement can stop
   #when ncol( r_adjust_matrix ) equals 2, that is , b equals (ndim-2)
   #replace lower.tri and upper.tri elements of r_matrix which equal 1 by 0.999
   r_matrix[lower.tri(r_matrix)] <- 0
   upper_elements <- r_matrix[upper.tri(r_matrix)]
   upper_elements[upper_elements==1] <- 0.999
   upper_elements[upper_elements==-1] <- -0.999
   r_matrix[upper.tri(r_matrix)] <- upper_elements
   r_matrix <- r_matrix + t(r_matrix)
   diag(r_matrix) <- 1
   r_adjust_matrix <- r_matrix
   b<- 0
   while (1) {
   chol_decompose_val <- try(chol(r_adjust_matrix), silent = TRUE)
      #Test whether the Cholesky decomposition fail
   if( class(chol_decompose_val) == "try-error" ) {
   b <- b + 1
   r_adjust_matrix <- r_matrix[1:(ndim - b),1:(ndim - b)]
   } else break
   }
   # compute simultaneously for (1:(ndim-b)) , rest shoud be computed separately
  mk_correlated_standard_norm <- w_k_bar[ ,1:(ndim - b)] %*% chol_decompose_val
  y_estimator <- matrix(rep(0,ndim * ndim ), nrow = ndim)
  diag(y_estimator) <- 1
  if (mk < 60)
    y_estimator[1:(ndim - b),1:(ndim - b)] <- y_crude_estimator(mk_correlated_standard_norm,
                                                           invcdfnames[1:(ndim - b)],
                                                           paramslists[1:(ndim - b)])
  else
    y_estimator[1:(ndim - b),1:(ndim - b)] <- y_control_estimator(mk_correlated_standard_norm,
                                                            r_matrix[1:(ndim - b),1:(ndim - b)],
                                                           invcdfnames[1:(ndim - b)],
                                                           paramslists[1:(ndim - b)])
  if (b != 0) {
     for (i in ndim:(ndim - b + 1) )
       for (j in 1:(i - 1)) {
         #The Cholesky decomposition will always succeed here for 1->0.99 transformation
         r_mat<- matrix(c(1,r_matrix[j,i],r_matrix[j,i],1), nrow = 2)
         mk_correlated_standard_norm<- w_k_bar[ ,c(j,i)] %*% chol(r_mat)
         if (mk < 60)
         tmp <- y_crude_estimator(mk_correlated_standard_norm,
                                  invcdfnames[c(i,j)], paramslists[c(i,j)])[1,2]
         else
         tmp <- y_control_estimator(mk_correlated_standard_norm, r_mat,
                                    invcdfnames[c(i,j)], paramslists[c(i,j)])[1,2]
         y_estimator[j,i] <- y_estimator[i,j] <- tmp
       }
  }
  y_estimator
}
