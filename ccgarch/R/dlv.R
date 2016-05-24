# the gradient of the GARCH part of the log-likelihood function
dlv <- function(u, a, A, B, model){
   ndim <- dim(u)[2]
   dl <- analytical.grad(a, A, B, diag(ndim), u, model=model)
   dl.row <- dim(dl)[1] - ndim*(ndim-1)/2
   dl[1:dl.row,]
}
