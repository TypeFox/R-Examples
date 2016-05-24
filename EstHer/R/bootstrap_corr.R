bootstrap_corr <- function(Y,Z,K,eta_hat,sigma2_hat,level,nb_cores)
{
N=ncol(Z)
n=length(Y)
res_eig=eigen(prod_cpp(Z)/N)
lambda=res_eig$values
O=res_eig$vectors
Y_tilde=t(O)%*%Y
gamma_vec_moins=(eta_hat*sigma2_hat*lambda+(1-eta_hat)*sigma2_hat)^(-0.5)
gamma_vec_plus=(eta_hat*sigma2_hat*lambda+(1-eta_hat)*sigma2_hat)^(0.5)
gamma_chap_moins=diag(gamma_vec_moins)
gamma_chap_plus=diag(gamma_vec_plus)
Y_bis=gamma_chap_moins%*%Y_tilde

eta_chap_boot=c()
resample_parallel <- function(i) 
{
  cat('+')
  Y1=sample(Y_bis,replace=T)
  Y1_tilde=gamma_chap_plus%*%Y1
  Y_boot=O%*%Y1_tilde
  return(estim_herit(Y_boot,Z)$heritability)
}
eta_chap_boot<- unlist(mclapply(1:K,resample_parallel, mc.cores=nb_cores))
cat('\n')

eta_boot_ord=eta_chap_boot[order(eta_chap_boot,decreasing=F)]
cut_eta=floor(level*K/2)
CI_low=eta_boot_ord[(cut_eta+1)]
CI_up=eta_boot_ord[(K-cut_eta)]
list(CI_up=min(CI_up,1),CI_low=max(CI_low,0))
}
