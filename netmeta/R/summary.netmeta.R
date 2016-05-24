summary.netmeta <- function(object,
                            level=object$level,
                            level.comb=object$level.comb,
                            comb.fixed=object$comb.fixed,
                            comb.random=object$comb.random,
                            reference.group=object$reference.group,
                            all.treatments=object$all.treatments,
                            warn=object$warn,
                            ...){
  
  
  if (!inherits(object, "netmeta"))
    stop("Argument 'object' must be an object of class \"netmeta\"")
  
  if (length(warn)==0){
    warn <- TRUE
  }
  
  k <- object$k
  m <- object$m
  n <- object$n
  Q <- object$Q
  
  ci.lab <- paste(round(100*level.comb, 1), "%-CI", sep="")
  ##
  ci.comp <- meta::ci(object$TE, object$seTE, level)
  ci.comp.nma.fixed <- meta::ci(object$TE.nma.fixed,
                                object$seTE.nma.fixed, level)
  ci.comp.nma.random <- meta::ci(object$TE.nma.random,
                                 object$seTE.nma.random, level)
  ci.f <- meta::ci(object$TE.fixed , object$seTE.fixed , level.comb)
  ci.r <- meta::ci(object$TE.random, object$seTE.random, level.comb)
  
  ci.comp$studlab <- object$studlab
  ci.comp$treat1 <- object$treat1
  ci.comp$treat2 <- object$treat2
  ##
  ci.comp.nma.fixed$studlab <- object$studlab
  ci.comp.nma.fixed$treat1 <- object$treat1
  ci.comp.nma.fixed$treat2 <- object$treat2
  ci.comp.nma.fixed$leverage <- object$leverage.fixed
  ##
  ci.comp.nma.random$studlab <- object$studlab
  ci.comp.nma.random$treat1 <- object$treat1
  ci.comp.nma.random$treat2 <- object$treat2
  
  
  res <- list(comparison=ci.comp,
              comparison.nma.fixed=ci.comp.nma.fixed,
              comparison.nma.random=ci.comp.nma.random,
              fixed=ci.f, random=ci.r,
              studies=object$studies,
              narms=object$narms,
              k=k, m=m, n=n, Q=Q, df=object$df,
              tau=object$tau, I2=object$I2,
              sm=object$sm,
              ci.lab=ci.lab,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              level=level,
              level.comb=level.comb,
              seq=object$seq
              )
  
  if (reference.group!="" & missing(all.treatments))
    all.treatments <- FALSE
  
  if (reference.group !="")
    reference.group <- setref(reference.group, rownames(object$A.matrix))
  
  res$reference.group <- reference.group
  res$all.treatments <- all.treatments
  ##
  res$title   <- object$title
  
  res$call <- match.call()
  res$version <- packageDescription("netmeta")$Version
  ##
  class(res) <- "summary.netmeta"
  
  res
}
