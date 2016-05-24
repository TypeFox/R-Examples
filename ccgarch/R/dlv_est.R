# the gradient of the GARCH part of the log-likelihood function
dlv.est <- function(par, dvar, model){
   ndim <- dim(dvar)[2]
   In <- diag(ndim)
   par <- c(par, In[lower.tri(In)])
   para <- p.mat(par, model, ndim)
   dl <- analytical.grad(para$a, para$A, para$B, In, dvar, model=model)
   dl.row <- dim(dl)[1] - ndim*(ndim-1)/2
   rowSums(dl[1:dl.row,])
}
