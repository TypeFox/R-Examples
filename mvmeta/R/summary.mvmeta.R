###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
summary.mvmeta <-
function(object, ci.level=0.95, ...) {
#
################################################################################
#
  if(ci.level<=0||ci.level>=1) stop("'ci.level' must be within 0 and 1")
#
  coef <- object$coefficients
  vcov <- object$vcov
  dim <- object$dim
  Psi <- object$Psi
  lab <- object$lab
#  
###########################################################################
# FIXED EFFECTS ESTIMATES
#
  # COMPUTE STATISTICS FOR FIXED EFFECTS
  coef <- as.numeric(coef)
  coef.se <- sqrt(diag(vcov))
  zval <- coef/coef.se
  zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
  pvalue <- 2*(1-pnorm(abs(zval)))
  ci.lb <- coef-zvalci*coef.se
  ci.ub <- coef+zvalci*coef.se
  cilab <- paste(signif(ci.level,2)*100,"%ci.",c("lb","ub"),sep="")
#
  # GENERATE TABLE AS MATRIX
  tabfixed <- cbind(coef,coef.se,zval,pvalue,ci.lb,ci.ub)
   dimnames(tabfixed) <- list(if(dim$k>1L) colnames(vcov) else lab$p,
     c("Estimate","Std. Error","z","Pr(>|z|)",cilab))
#
  # CORRELATION MATRIX OF FIXED EFFECTS
  corFixed <- vcov/outer(coef.se,coef.se)
#  
###########################################################################
# RANDOM EFFECTS ESTIMATES (FOR RANDOM-EFFECTS MODELS)
#  
  corRandom <- if(object$method!="fixed") {
    # SD OF EACH RANDOM EFFECT
    ran.sd <- sqrt(diag(Psi))
    # CORRELATION MATRIX OF RANDOM EFFECTS
    Psi/outer(ran.sd,ran.sd)
  } else NULL
#  
###########################################################################
# QTEST STATISTICS
#
  qstat <- unclass(qtest(object))
#
###########################################################################
# 
  # DEFINE THE LIST
  keep <- match(c("vcov","Psi","bscov","df.res","rank","logLik","converged",
    "niter","negeigen","method","dim","df","lab","na.action","call","terms"),
    names(object),0L)
  out <- c(list(coefficients=tabfixed),object[keep],list(AIC=AIC(object),
    BIC=BIC(object),corFixed=corFixed,corRandom=corRandom,qstat=qstat,
    ci.level=ci.level))
#
  class(out) <- "summary.mvmeta"
#
  out
}
