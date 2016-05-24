limitmeta <- function(x,
                      method.adjust="beta0",
                      level=x$level, level.comb=x$level.comb,
                      backtransf=x$backtransf,
                      title=x$title, complab=x$complab, outclab=x$outclab){
  
  meta:::chkclass(x, "meta")
  
  
  TE <- x$TE
  seTE <- x$seTE
  tau <- x$tau
  w.random <- x$w.random
  k <- x$k
  Q <- x$Q
  sm <- x$sm
  ##
  seTE.tau  <- sqrt(1/w.random)
  ##
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  lower.random <- x$lower.random
  upper.random <- x$upper.random
  zval.random <- x$zval.random
  pval.random <- x$pval.random
  
  
  imeth <- charmatch(tolower(method.adjust),
                     c("beta0", "betalim", "mulim"), nomatch = NA)
  ##
  if(is.na(imeth))
    stop("Argument 'method.adjust' should be \"beta0\", \"betalim\", or \"mulim\"")
  ##
  method.adjust <- c("beta0", "betalim", "mulim")[imeth]
  
  
  ##
  ## Radial plot, slope best fit (beta-F)
  ##
  reg.f <- radialregression(TE, seTE, k)
  
  
  ##
  ## Generalized radial plot, slope best fit (beta-R)
  ##
  reg.r <- radialregression(TE, seTE.tau, k)
  ##
  alpha.r <- reg.r$intercept
  beta.r  <- reg.r$slope
  
  
  ##
  ## Conduct limit meta-analysis
  ##
  TE.limit <- beta.r + sqrt(tau^2/seTE.tau^2)*(TE - beta.r)
  seTE.limit <- seTE / 1 # 1 == "Infinity"
  ##
  m.lim <- metagen(TE.limit, seTE.limit, sm=sm)
  ##
  reg.l <- radialregression(m.lim$TE, m.lim$seTE, k)
  
  
  ##
  ##
  ## Conduct adjustment methods
  ##
  ##
  if (method.adjust=="beta0"){
    ##
    ## Expectation (beta-0)
    ##
    TE.adjust   <- as.vector(beta.r + tau*alpha.r)
    seTE.adjust <- as.vector(1/sd(sqrt(1/seTE^2))/sqrt(k-1))
  }
  else{
    if (method.adjust=="mulim"){
      ##
      ## Limit radial plot, slope through origin (mu-lim)
      ##
      TE.adjust   <- m.lim$TE.fixed
      seTE.adjust <- m.lim$seTE.fixed
    }
    else if (method.adjust=="betalim"){
      ##
      ## Limit radial plot, slope best fit (beta-lim)
      ##
      TE.adjust   <- as.vector(reg.l$slope)
      seTE.adjust <- as.vector(reg.l$se.slope)
    }
  }
  ##
  ci.adjust <- ci(TE.adjust, seTE.adjust, level=level.comb)
  ##
  lower.adjust <- ci.adjust$lower
  upper.adjust <- ci.adjust$upper
  zval.adjust <- ci.adjust$z
  pval.adjust <- ci.adjust$p
  ##
  if (inherits(x, c("metaprop"))){
    zval.adjust <- NA
    pval.adjust <- NA
  }
  
  
  ##
  ## Only recalculate RE confidence interval if argument 'level.comb'
  ## is not missing
  ##
  if (!missing(level.comb)){
    ci.r <- ci(TE.random, seTE.random, level=level.comb)
    ##
    lower.random <- ci.r$lower
    upper.random <- ci.r$upper
  }
  
  
  ## Ruecker et al. (2011), Biostatistics, pages 133-134
  ##
  Q.resid <- reg.f$sigma^2*(k-2)
  Q.small <- Q - Q.resid
  ##
  G.squared=1-reg.l$r.squared
  
  
  res <- list(TE=TE,
              seTE=seTE,
              ##
              TE.limit=TE.limit,
              seTE.limit=seTE.limit,
              ##
              studlab=x$studlab,
              ##
              TE.random=TE.random,
              seTE.random=seTE.random,
              lower.random=lower.random,
              upper.random=upper.random,
              zval.random=zval.random,
              pval.random=pval.random,
              w.random=w.random,
              tau=tau,
              ##
              TE.adjust=TE.adjust,
              seTE.adjust=seTE.adjust,
              lower.adjust=lower.adjust,
              upper.adjust=upper.adjust,
              zval.adjust=zval.adjust,
              pval.adjust=pval.adjust,
              ##
              alpha.r=alpha.r,
              beta.r=beta.r,
              ##
              Q=Q,
              Q.small=Q.small,
              Q.resid=Q.resid,
              G.squared=G.squared,
              ##
              level=level,
              level.comb=level.comb,
              ##
              k=k,
              sm=sm,
              method.adjust=method.adjust,
              ##
              title=title,
              complab=complab,
              outclab=outclab,
              ##
              call=match.call(),
              x=x)
  
  
  res$backtransf <- backtransf
  
  res$version <- utils::packageDescription("metasens")$Version
  
  class(res) <- c("limitmeta")
  
  res
}
