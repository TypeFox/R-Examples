para_trans_inv <-
function(A) {
	
  Q = dim(A)[1];
  B = exp(A);
  D = 1/(1+apply(B,1,sum));
  transmat = cbind(B*matrix(D,Q,Q-1),D)

}
