loglik_nh_inp <-
function(Y,covar,xi, nh_transition) {
  M = dim(xi)[1]
  T = dim(xi)[3]
  N.samples = dim(xi)[4]

  res = deplie2(Y,M);
  transmat = para_trans_inv(res$trans);
  f = 0
  f1 =0
  for (ex in (1:N.samples)) { 
  	transition=nh_transition(array(covar[ ,ex,],c(T,1,dim(covar)[3])),res$par,transmat);     
    w = which(transition<1e-15,arr.ind=T)
    if (length(w[,1])>0) {for (k in 1:length(w[,1])) {
    		transition[w[k,1],w[k,2],w[k,3]] = 1e-15}}
    ltr = log(transition)
    f = f+sum(ltr*xi[,,,ex])
    }
  f=-f;
}
