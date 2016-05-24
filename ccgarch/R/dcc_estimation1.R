#*****************************************************************************************************************
# The 1st stage DCC estimation.
dcc.estimation1 <- function(dvar, a, A, B, model, method="BFGS"){
   nobs <- dim(dvar)[1]
   ndim <- dim(dvar)[2]
   if(model=="diagonal"){
    init <- sqrt(c(a, diag(A), diag(B)))
   } else {
    init <- sqrt(c(a, as.vector(A), as.vector(B)))
   }
   step1 <- optim(par=init, fn=loglik.dcc1, method=method, control=list(maxit=10^5, reltol=1e-5), dvar=dvar, model=model)
   step1
}
