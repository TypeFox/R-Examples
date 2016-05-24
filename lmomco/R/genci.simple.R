"genci" <- function(...) {
   warning("function name genci() is deprecated, please use ",
           "genci.simple() instead")
   genci.simple(...)
}

"genci.simple" <-
function(para,n, f=NULL, level=0.90, edist="gno",
         nsim=1000, expand=FALSE, verbose=FALSE, showpar=FALSE, quiet=FALSE) {
  if(is.null(f)) f <- nonexceeds()
  if(! check.fs(f)) {
     warning("The provided nonexceedance probabilities are invalid")
     return()
  }
  if(! are.par.valid(para)) {
     warning("The distribution parameters are invalid")
     return()
  }
  if(! check.fs(level) | level == 1) { # invalid confidence level
      warning("argument 'ci' is not in [0,1)")
      return()
  }

  num.Fs <- length(f)
  ci_low <- ci_tru <- ci_hi  <- vector(mode = 'numeric', length=num.Fs)
  ci_l1  <- ci_l2  <- ci_t3  <- ci_t4  <- ci_t5  <- ci_low
  ci_md  <- ci_mu  <- ci_var <- ci_skw <- ci_low

  if(! quiet) cat(c(num.Fs,"-"), sep="")
  for(i in seq(1,num.Fs)) {
    CI <- qua2ci.simple(f[i], para, n, level=level, edist=edist, nsim=nsim,
                        verbose=verbose, showpar=showpar, empdist=TRUE)
    if(CI$ifail > 0) {
       ci_low[i] <- ci_tru[i] <- ci_hi[i]  <- NA
       ci_l1[i]  <- ci_l2[i]  <- ci_t3[i]  <- ci_t4[i] <- ci_t5[i]  <- NA
       ci_md[i]  <- ci_mu[i]  <- ci_var[i] <- ci_skw[i] <- NA
       next
    }
    ci_low[i] <- CI$lwr
    ci_tru[i] <- CI$true
    ci_hi[i]  <- CI$upr
    ci_l1[i]  <- CI$elmoms$lambdas[1]
    ci_l2[i]  <- CI$elmoms$lambdas[2]
    ci_t3[i]  <- CI$elmoms$ratios[3]
    ci_t4[i]  <- CI$elmoms$ratios[4]
    ci_t5[i]  <- CI$elmoms$ratios[5]
    ci_md[i]  <- CI$empdist$median
    ci_mu[i]  <- CI$empdist$epmoms$moments[1]
    ci_var[i] <- CI$empdist$epmoms$moments[2]^2 # notice that a variance is computed
    ci_skw[i] <- CI$empdist$epmoms$ratios[3]
    if(! quiet) cat(c(num.Fs-i,"-"), sep="")
  }
  if(! quiet) cat("\n")

  cis <- data.frame(nonexceed=f,
                    lwr=ci_low, true=ci_tru, upr=ci_hi, qua_med=ci_md,
                    qua_mean=ci_mu, qua_var=ci_var, qua_lam2=ci_l2)

  lmr <- data.frame(lambda1=ci_l1, lambda2=ci_l2,
                    tau3=ci_t3, tau4=ci_t4, tau5=ci_t5)
  pmr <- data.frame(mu=ci_mu, var=ci_var, skw=ci_skw)
  if(expand) {
    return(list(limits=cis, parent=para,
                edist=edist, elmoms=lmr, epmoms=pmr, epara=CI$epara,
                ifail=CI$ifail, ifailtext=CI$ifailtext))
  }
  else {
    return(cis)
  }
}
