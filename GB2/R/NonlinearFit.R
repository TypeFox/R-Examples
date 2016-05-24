nlsfit.gb2 <- function(med, ei4, par0=c(1/ei4[4],med,1,1), cva=1, bound1 = par0[1]*max(0.2,1-2*cva), bound2=par0[1]*min(2,1+2*cva),ei4w=1/ei4){
  
  a0 <- par0[1]
  b0 <- par0[2]
  p0 <- par0[3]
  q0 <- par0[4]
  
# nls fit 1
  fit1 <- nls(ei4 ~ main2.gb2(0.6, a, 1, ap, aq)[3:6], weights = ei4w, start = list( a=a0, ap=a0*p0, aq=a0*q0), 
          trace=FALSE, algorithm = "port", lower = c(bound1, 1, 2), upper = c(bound2, 100, 100), 
          control = nls.control(maxiter = 1000, tol = 1e-06, minFactor = 1/1024, printEval = FALSE, warnOnly = TRUE))

  an <- coef(fit1)[[1]]
  pn <- coef(fit1)[[2]]/an
  qn <- coef(fit1)[[3]]/an

# nls fit 2
  fit2 <- nls(med ~ qgb2(0.5, an, b, pn, qn), start = list(b=b0), 
          trace=FALSE, algorithm = "port", lower = c(0.01), upper = c(2*b0), 
          control = nls.control(maxiter = 1000, tol = 1e-06, minFactor = 1/1024, printEval = FALSE, warnOnly = TRUE))

  bn <- coef(fit2)[[1]]
  pars <- c(an,bn,pn,qn)
  return(list(pars, fit1, fit2))
}
