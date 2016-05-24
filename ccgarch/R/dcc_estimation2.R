#*****************************************************************************************************************
# The 2nd step DCC estimation.
dcc.estimation2 <- function(dvar, para, gradient=0){ # dvar must be standardised residuals
   resta <- rbind(c(-1, -1), diag(2))
   restb <- c(-1, 0, 0)
   if(gradient!=0){
      step2 <- constrOptim(theta=para, f=loglik.dcc2, grad=grad.dcc2, ui=resta, ci=restb, mu=1e-5, dvar=dvar)
   } else {
      step2 <- constrOptim(theta=para, f=loglik.dcc2, grad=NULL, ui=resta, ci=restb, mu=1e-5, dvar=dvar)
   }
   step2
}
