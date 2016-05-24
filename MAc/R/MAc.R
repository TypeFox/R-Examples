##==================             MAc             ================##      
##==================  Meta-Analysis Correlations ================##

# Package created by AC Del Re & William T. Hoyt
# This package contains all the relevant functions to conduct a
# correlational meta-analysis using standard procedures as
# described in Cooper, Hedges, & Valentine's Handbook of
# Research Synthesis and Meta-Analysis (2009).

# suggests('ggplot2')
# suggests('metafor')
# suggests('irr')
# suggests('R2wd')
# enhances('compute.es')

##=== New Functions ===##

# Overhauled entire package 06.29.10

mareg <- function(formula, var, data, ztor = FALSE, method = "REML", subset, ...) {
  UseMethod("mareg")
}

print.mareg <- function (x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    coef <- as.vector(x$b)
    names(coef) <- rownames(x$b)
  print.default(format(coef, digits = digits), print.gap = 2, 
            quote = FALSE)   
  #print.default(coef)
    invisible(x)
}
    
summary.mareg <- function(object, ...) {
  b <- object$b
  se <- object$se
  z <- object$zval
  ci.lower <- object$ci.lb
  ci.upper <- object$ci.ub
  p <- object$pval
  QE <- object$QE
  QEp <- object$QEp
  QE.df <- round(object$k - object$p, 0)
  QM <- object$QM
  QMp <- object$QMp
  QM.df <- round(object$m, 0)
  tau2_empty <- object$tau2_empty
  tau2 <- object$tau2
  R2 <- object$R2
  table <- cbind(b, se, z, ci.lower, ci.upper, p)
  colnames(table) <- c("estimate", "se", "z", "ci.u",
    "ci.l", "Pr(>|z|)")
  #rownames(table) <- object$predictors
  table2 <- cbind(QE, QE.df, QEp, QM, QM.df, QMp)
  colnames(table2) <- c("QE", "QE.df", "QEp", "QM", "QM.df", "QMp")
  model <- list(coef = table, fit = table2, tau2_empty=tau2_empty,
   tau2=tau2, R2=R2)
  rownames(model$fit) <-NULL
  class(model) <- "summary.mareg"
  return(model)
}


print.summary.mareg <- function(x, ...) {
  cat("\n Note: 
       If using r and the variance of r as input, be sure to 
       leave the default ztor = FALSE, otherwise the output will be
       inaccurate. If using z' (Fisher's z) and the variance of z'
       as the input, changing ztor = TRUE will convert z' back to r 
       (for interpretive purposes) after all analyses have been 
       conducted.", "", "\n","----","\n")
  cat("\n Model Results:", "", "\n", "\n")
  printCoefmat(x$coef, signif.stars = TRUE)
  cat("\n Heterogeneity & Fit:", "", "\n" ,"\n")
  print(round(x$fit,4))
  invisible(x)
}

r2 <- function(x) {
   UseMethod("r2")
}

r2.mareg <- function(x) {
  cat("\n Explained Variance:", "", "\n")
  cat("\n Tau^2      (total):", round(x$tau2_empty,4))
  cat("\n Tau^2      (model):", round(x$tau2,4))
  cat("\n R^2 (% expl. var.):", round(x$R2,2), "\n")
}

# formula based function for meta-regression!! 
mareg.default <- function(formula, var, data, ztor = FALSE, method = "REML",  
  subset, ...) {
    require('metafor', quietly=TRUE)  # requires metafor's rma function
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    args <- match(c("formula", "var", "data", "subset"),
      names(mf), 0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    terms <- attr(mf, "terms")
    yi <- model.response(mf)
    vi <- model.extract(mf, "var")
    mods <- model.matrix(terms, mf, contrasts)
    intercept <- which(colnames(mods) == "(Intercept)")
    if(length(intercept > 0)) mods <- mods[ , -intercept]
    model0 <- rma(yi=yi, vi=vi, method=method, ...)
    model <- rma(yi,vi=vi,mods=mods, method=method, ...)
    predictors <- colnames(mods)
    if(ztor == TRUE) {
      model0$var.z <- diag(model0$vb)
      model0$var <- (exp(2*model0$var.z)-1)/(exp(2*model0$var.z) + 1)
      model0$se <- sqrt(model0$var)
      model0$pval <- (exp(2*model0$pval)-1)/(exp(2*model0$pval) + 1)
      model0$b <- (exp(2*model0$b)-1)/(exp(2*model0$b) + 1)
      model0$ci.lb <- (exp(2*model0$ci.lb)-1)/(exp(2*model0$ci.lb) + 1)
      model0$ci.ub <- (exp(2*model0$ci.ub)-1)/(exp(2*model0$ci.ub) + 1)
      model0$zval <- (exp(2*model0$zval)-1)/(exp(2*model0$zval) + 1)
      model$var.z <- diag(model$vb)
      model$var <- (exp(2*model$var.z)-1)/(exp(2*model$var.z) + 1)
      model$se <- sqrt(model$var)
      model$pval <- (exp(2*model$pval)-1)/(exp(2*model$pval) + 1)
      model$b <- (exp(2*model$b)-1)/(exp(2*model$b) + 1)
      model$ci.lb <- (exp(2*model$ci.lb)-1)/(exp(2*model$ci.lb) + 1)
      model$ci.ub <- (exp(2*model$ci.ub)-1)/(exp(2*model$ci.ub) + 1)
      model$zval <- (exp(2*model$zval)-1)/(exp(2*model$zval) + 1)
    } 
    model$out <- with(model, list(estimate=b, se=se, z=zval, ci.l=ci.lb,
      ci.u=ci.ub,  "Pr(>|z|)"=pval, predictors=predictors))
    model$out2 <- with(model, data.frame(QE, QEp, QM, QMp))
    model$call <- call
    model$tau2_empty <- model0$tau2
    model$R2 <-  1-(model$tau2 / model$tau2_empty)
    class(model) <- c("mareg", "rma.uni", "rma")
    return(model)
} 

# function to print objects to Word
wd <- function(object, get = FALSE, new = FALSE, ...) {
  UseMethod("wd")
}

wd.default <- function(object, get = FALSE, new = FALSE, ...) {
  require('R2wd', quietly = TRUE)
  get <- get
  open <- wdGet(get)	# If no word file is open, it will start a new one
  if(new == TRUE){
    new <- wdNewDoc("c:\\Temp.doc")
  }
  else{ 
    new <- NULL
  }  	
  wd1 <- round(data.frame(estimate=object$b, se=object$se, z=object$zval,
   ci.lower=object$ci.lb, ci.upper=object$ci.ub,  "p"=object$pval), 4)
  title <- wdHeading(level = 2, " Model Results:")
  wd2 <-  round(data.frame(QE=object$QE, QEp=object$QEp, QM=object$QM, QMP=object$QMp), 4)
  obj1 <- wdTable(wd1, ...)
  title2 <- wdHeading(level = 2, "Heterogeneity & Fit:")
  obj2 <- wdTable(wd2, ...)
  out <- (list(open, new, title, obj1, title2, obj2))
  class(out) <- "wd"
  return(out)
}
wd.mareg <- function(object, get = FALSE, new = FALSE, ...) {
  require('R2wd', quietly = TRUE)
  get <- get
  open <- wdGet(get)	# If no word file is open, it will start a new one
  if(new == TRUE){
    new <- wdNewDoc("c:\\Temp.doc")
  }
  else{ 
    new <- NULL
  }  	
  QE.df <- round(object$k - object$p, 0)
  QM.df <- round(object$m, 0)
  wd1 <- round(data.frame(estimate=object$b, se=object$se, z=object$zval,
   ci.l=object$ci.lb, ci.u=object$ci.ub,  "p"=object$pval), 4)
  title <- wdHeading(level = 2, " Model Results:")
  wd2 <-  round(data.frame(QE=object$QE, QE.df = QE.df,
            QEp=object$QEp, QM=object$QM, QM.df=QM.df, QMp=object$QMp), 4)
  obj1 <- wdTable(wd1)
  title2 <- wdHeading(level = 2, "Heterogeneity & Fit:")
  obj2 <- wdTable(wd2, ...)
  out <- (list(open, new, title, obj1, title2, obj2))
  class(out) <- "wd.mareg"
  return(out)
}

wd.omni <- function(object, get = FALSE, new = FALSE, ...) {
  require('R2wd', quietly = TRUE)
  get <- get
  open <- wdGet(get)	# If no word file is open, it will start a new one
  if(new == TRUE){
    new <- wdNewDoc("c:\\Temp.doc")
  }
  else{ 
    new <- NULL
  }  	
    k <- object$k
    estimate <- object$estimate
    var <- object$var
    se <- object$se
    ci.l <- object$ci.l
    ci.u <- object$ci.u 
    z <- object$z
    p <- object$p
    Q <- object$Q
    Q.df <- object$df.Q
    Qp <- object$Qp
    I2 <- object$I2
    results1 <- round(data.frame(estimate, var, se, ci.l, ci.u, z, p),4)
    #results <- formatC(table, format="f", digits=digits)
    results1$k <- object$k
    results1 <- results1[c(8,1:7)]
    results2 <- round(data.frame(Q, Q.df, Qp),4)
    results2$I2 <- object$I2
    title <- wdHeading(level = 2, " Omnibus Model Results:")
    obj1 <- wdTable(results1)
    title2 <- wdHeading(level = 2, "Heterogeneity:")
    obj2 <- wdTable(results2)
    out <- (list(open, new, title, obj1, title2, obj2))
    class(out) <- "wd.omni"
    return(out)
}

wd.macat <- function(object, get = FALSE, new = FALSE, ...) {
  require('R2wd', quietly = TRUE)
  get <- get
  open <- wdGet(get)	# If no word file is open, it will start a new one
  if(new == TRUE){
    new <- wdNewDoc("c:\\Temp.doc")
  }
  else{ 
    new <- NULL
  }  	
    x1 <- object$Model
    x2 <- object$Heterogeneity
    mod <- x1$mod
    k <- x1$k
    estimate <- x1$estimate
    var <- x1$var
    se <- x1$se
    ci.l <- x1$ci.l
    ci.u <- x1$ci.u 
    z <- x1$z
    p <- x1$p
    Q <- x1$Q
    df <- x1$df
    p.h <- x1$p.h
    I2 <- x1$I2
    Qoverall <- x2$Q
    Qw <- x2$Qw
    Qw.df <- x2$df.w
    Qw.p <- x2$p.w
    Qb <- x2$Qb
    Qb.df <- x2$df.b
    Qb.p <- x2$p.b
    results1 <- round(data.frame(estimate, var, se, ci.l, ci.u, z, p,
      Q, df, p.h),4)
    #results <- formatC(table, format="f", digits=digits)
    results1$k <- x1$k
    results1$I2 <- x1$I2
    results1$mod <- x1$mod
    results1 <- results1[c(13,11, 1:10,12)]
    results2 <- round(data.frame(Q=Qoverall, Qw, Qw.df, Qw.p, Qb, Qb.df, Qb.p),4)
    results1[2:12] <-round(results1[2:12],3)
    title <- wdHeading(level = 2, " Categorical Moderator Results:")
    obj1 <- wdTable(results1)
    title2 <- wdHeading(level = 2, "Heterogeneity:")
    obj2 <- wdTable(results2)
    out <- (list(open, new, title, obj1, title2, obj2))
    class(out) <- "wd.macat"
    return(out)
}

##============ WITHIN STUDY AGGREGATION OF EFFECT SIZES =============##

# Automated within-study effect size aggregation function accounting for 
# dependencies among effect sizes g (unbiased estimate of d).
# Required inputs are n.1 (treatment sample size), n.2 (comparison sample size)
# id (study id), g (unbiased effect size).

aggrs <- function(r, cor = .50) {
  # Intermediate level function to compute weighted aggregation of effect sizes 
  k <- length(r)
   rbar <- cor
   r.xY <- sum(r)/(1 * sqrt(k  +  k*(k - 1) * rbar))
   return(r.xY)
}

# automated agg 
agg <- function(id, r, n, cor = .50, mod=NULL, data) {
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("id", "r", "n", "mod", "cor", "data"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  #mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.id <- mf[[match("id", names(mf))]]
  id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  mf.r <- mf[[match("r", names(mf))]]
  r <- eval(mf.r, data, enclos = sys.frame(sys.parent()))
  mf.n <- mf[[match("n", names(mf))]]
  n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  if(is.null(mod)){

    st <- unique(id)
    out <- data.frame(id=st)
    for(i in 1:length(st)) { 
    out$id[i] <- st[i]
    out$r[i] <- aggrs(r = r[id==st[i]], cor)
    out$n[i] <- round(mean(n[id==st[i]]),0)
    }
  }
  if(!is.null(mod)) {
  #additional agg functions for mods, etc
    st <- as.factor(id)
    st <- unique(st)
    um <- unique(mod)
    out <- data.frame(id = rep(st, rep(length(um), length(st))))
    out$mod <- rep(um, length(st))
    for (i in 1:length(st)) {
        for (j in 1:length(um)) {
            ro <- (i - 1) * length(um) + j
            m1 <- match(id, st[i], nomatch = 0)
            m2 <- match(mod, um[j], nomatch = 0)
            num <- sum(m1 * m2)
            out$r[ro] <- ifelse(num == 0, NA, aggrs(r = r[id == 
                st[i] & mod == um[j]], cor))
            out$n[ro] <- round(mean(n[id == st[i] & mod == um[j]]), 
                0)
        }
    }
    out <- out[is.na(out$r) == 0, ]
    } 
  return(out)
}  


##=== Add Fixed and Random Effects Weights ===##
# Required input is a data.frame with column names id (study id), 
# g (unbiased standardized mean diff ES),  and n.1 (group 1 sample size),
# n.2 (group 2 sample size).
 

# required inputs are g, var.g and will output weights, ci, etc
wgts <-  function(es, var.es, data) {  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("es", "var.es", "data"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  #mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.es <- mf[[match("es", names(mf))]]
  meta <- data
  es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var.es <- mf[[match("var.es", names(mf))]]
  var.es <- eval(mf.var.es, data, enclos = sys.frame(sys.parent()))
  meta$l.ci95 <- es-1.96*sqrt(var.es)     #create random ci for each study
  meta$u.ci95 <- es + 1.96*sqrt(var.es)
  meta$z.score <- es/sqrt(var.es)
  meta$wi <-  1/var.es  # computing weight for each study
  meta$wiTi <- meta$wi*es  # used to calculate omnibus
  meta$wiTi2 <- meta$wi*(es)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(es))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  meta$tau <- (Q-(k - 1))/comp  # Level 2 variance
  meta$var.tau <- meta$tau + var.es  # Random effects variance (within study var + between var) 
  meta$wi.tau <- 1/meta$var.tau  # Random effects weights
  meta$wiTi.tau <- meta$wi.tau*es
  meta$wiTi2.tau <- meta$wi.tau*(es)^2
  return(meta)
}


##================= FIXED AND RANDOM EFFECTS OMNIBUS ===============##
# Function to calculate fixed and random effects omnibus effect size for g,  
# outputing omnibus effect size,  variance,  standard error,  upper and lower 
# confidence intervals,  and heterogeneity test.
# Required input is a data.frame with column names id (study id), 
# g (unbiased standardized mean diff ES),  and n.1 (group 1 sample size),
# n.2 (group 2 sample size).

omni <-  function(es, var, data, type="weighted", method = "random", ztor = FALSE) {
  # Computes fixed and random effects omnibus effect size for correlations.  
  # Args:
  #   es:    vector of r or z' effect sizes
  #   var:  vector of variances of r or z'.
  #   data: data.frame where r (or z') and var for each study are located.
  #   type:  "weighted" or "unweighted". "weighted" is the default. Use the 
  #   unweighted variance method only if Q is rejected and is very large relative to k.   
  # Returns:
  #   Fixed or random (Restricted Maximal Likelihood) effects omnibus effect size, variance, standard error, 
  #   upper and lower confidence intervals, p-value, Q (heterogeneity test), I2
  #   (I-squared--proportion of total variation in tmt effects due to heterogeneity 
  #   rather than chance). 
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("es", "var", "data", "type", "method"),
    names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.es <- mf[[match("es", names(mf))]]
  es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var.es <- mf[[match("var", names(mf))]]
  var.es <- eval(mf.var.es, data, enclos = sys.frame(sys.parent()))
  l.ci95 <- es-1.96*sqrt(var.es)     #create random ci for each study
  u.ci95 <- es + 1.96*sqrt(var.es)
  z.score <- es/sqrt(var.es)
  wi <-  1/var.es  # computing weight for each study
  wiTi <- wi*es  # used to calculate omnibus
  wiTi2 <- wi*(es)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(es))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  tau <- (Q-(k - 1))/comp  # Level 2 variance
  var.tau <- tau + var.es  # Random effects variance (within study var + between var) 
  wi.tau <- 1/var.tau  # Random effects weights
  wiTi.tau <- wi.tau*es
  wiTi2.tau <- wi.tau*(es)^2
  k <- length(!is.na(es)) # number of studies
  df <- k-1 
  sum.wi <- sum(wi, na.rm=TRUE)  # used to calculate omnibus
  sum.wiTi <- sum(wiTi, na.rm=TRUE)  # used to calculate omnibus
  sum.wiTi2 <- sum(wiTi2, na.rm=TRUE)  # used to calculate omnibus
  T.agg <- sum.wiTi/sum.wi  # omnibus g 
  var.T.agg <- 1/sum.wi  # omnibus var.g
  se.T.agg <- sqrt(var.T.agg) 
  z.value  <-  T.agg/se.T.agg
  p.value <-  2*pnorm(abs(z.value), lower.tail=FALSE)
  lower.ci <- T.agg-1.96*se.T.agg
  upper.ci <- T.agg + 1.96*se.T.agg
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  # FE homogeneity test
  I2 <- (Q-(k-1))/Q  # I-squared 
  I2 <- ifelse(I2<0, 0, I2)        
  I2 <- paste(round(I2*100, 4),  "%",  sep="")                        
  p.homog <- pchisq(Q, df, lower=FALSE)  # <.05 = sig. heterogeneity  
  # random effects #
  sum.wi2 <- sum(wi^2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  sum.wi.tau <- sum(wi.tau, na.rm=TRUE)
  sum.wiTi.tau <- sum(wiTi.tau, na.rm=TRUE)
  sum.wiTi2.tau <- sum(wiTi2.tau, na.rm=TRUE)
  T.agg.tau <- sum.wiTi.tau/sum.wi.tau
  if(type == "weighted") {
    var.T.agg.tau <-  1/sum.wi.tau 
    se.T.agg.tau <- sqrt(var.T.agg.tau)
    z.valueR  <-  T.agg.tau/se.T.agg.tau
    p.valueR <-  2*pnorm(abs(z.valueR), lower.tail=FALSE)
    lower.ci.tau <- T.agg.tau-1.96*se.T.agg.tau
    upper.ci.tau <- T.agg.tau + 1.96*se.T.agg.tau
  }
  if(type == "unweighted") {  # unweighted variance method
    var.agg <- (sum(g^2)-sum(g)^2/k)/(k-1) #14.20
    #q.num <- (1/k)*sum(var.g)                                   
    unwgtvar.T.agg.tau <- var.agg-(1/k)*sum(var.g)  #14.22
    var.T.agg.tau <-  ifelse(unwgtvar.T.agg.tau <= 0, 0, unwgtvar.T.agg.tau)  #if var < 0,  its set to 0
    se.T.agg.tau <- sqrt(var.T.agg.tau)
    z.valueR  <-  T.agg.tau/se.T.agg.tau
    p.valueR <-  2*pnorm(abs(z.valueR), lower.tail=FALSE)
    lower.ci.tau <- T.agg.tau-1.96*se.T.agg.tau
    upper.ci.tau <- T.agg.tau + 1.96*se.T.agg.tau
  } 
  if(ztor == TRUE) {
      # fixed effects:
      var.T.agg <- (exp(2*var.T.agg)-1)/(exp(2*var.T.agg) + 1)
      se.T.agg <- sqrt(var.T.agg)
      p.value <- (exp(2*p.value)-1)/(exp(2*p.value) + 1)
      T.agg <- (exp(2*T.agg)-1)/(exp(2*T.agg) + 1)
      lower.ci <- (exp(2*lower.ci)-1)/(exp(2*lower.ci) + 1)
      upper.ci <- (exp(2*upper.ci)-1)/(exp(2*upper.ci) + 1)
      z.value <- (exp(2*z.value)-1)/(exp(2*z.value) + 1)
      # random effects:
      var.T.agg.tau <- (exp(2*var.T.agg.tau)-1)/(exp(2*var.T.agg.tau) + 1)
      se.T.agg.tau <- sqrt(var.T.agg.tau)
      p.valueR <- (exp(2*p.valueR)-1)/(exp(2*p.valueR) + 1)
      T.agg.tau <- (exp(2*T.agg.tau)-1)/(exp(2*T.agg.tau) + 1)
      lower.ci.tau <- (exp(2*lower.ci.tau)-1)/(exp(2*lower.ci.tau) + 1)
      upper.ci.tau <- (exp(2*upper.ci.tau)-1)/(exp(2*upper.ci.tau) + 1)
      z.valueR <- (exp(2*z.valueR)-1)/(exp(2*z.valueR) + 1)

    } 

  Fixed <- list(k=k, estimate=T.agg,  var=var.T.agg,  se=se.T.agg, 
                ci.l=lower.ci,  ci.u=upper.ci,  z=z.value,  p=p.value,
                Q=Q, df.Q=df, Qp=p.homog, I2=I2, call=call)
  Random <- list(k=k, estimate=T.agg.tau, var=var.T.agg.tau,  
                 se=se.T.agg.tau,  ci.l=lower.ci.tau, ci.u=upper.ci.tau, 
                 z=z.valueR, p=p.valueR, Q=Q, df.Q=df,  Qp=p.homog, 
                 I2=I2, call=call)
  if(method == "fixed") {
    omni.data <- Fixed
    row.names(omni.data) <- NULL
} 
  if(method == "random") {
    omni.data <- Random
    row.names(omni.data) <- NULL
}  
  class(omni.data) <-"omni"
  return(omni.data)
}

print.omni <- function (x, digits = 4, ...){
    k <- x$k
    estimate <- x$estimate
    var <- x$var
    se <- x$se
    ci.l <- x$ci.l
    ci.u <- x$ci.u 
    z <- x$z
    p <- x$p
    Q <- x$Q
    df.Q <- x$df.Q
    Qp <- x$Qp
    I2 <- x$I2
    results1 <- round(data.frame(estimate, var, se,ci.l, ci.u, z, p),4)
    #results <- formatC(table, format="f", digits=digits)
    results1$k <- x$k
    results1 <- results1[c(8,1:7)]
    results2 <- round(data.frame(Q, df.Q, Qp),4)
    results2$I2 <- x$I2
    cat("\n Note: 
       If using r and the variance of r as input, be sure to 
       leave the default ztor = FALSE, otherwise the output will be
       inaccurate. If using z' (Fisher's z) and the variance of z'
       as the input, changing ztor = TRUE will convert z' back to r 
       (for interpretive purposes) after all analyses have been 
       conducted.", "", "\n","----","\n")
    cat("\n Model Results:", "", "\n", "\n")
    print(results1)
    cat("\n Heterogeneity:", "", "\n","\n")
    print(results2)
    invisible(x)
}



##================= Categorical Moderator Analysis ================##

macat <- function (es, var, mod, data, method= "random", ztor = FALSE) {
    # Computes single predictor categorical moderator analysis. Computations
    # derived from chapter 15, Cooper et al. (2009). 
    # Args:
    #   meta: data.frame with g (standardized mean diff),var.g (variance of g)
    #   mod: Categorical moderator variable used for moderator analysis.    
    # Returns:
    #   Fixed or random effects moderator means per group, k per group, 
    #   95% confidence intervals, z-value, p-value, variances, 
    #   standard errors, Q, df(Q), and I-squared.
    #   Also outputs Q-statistic, Q-within & between, df(Qw & Qb), and 
    #   homogeneity p-value within & between levels.
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    args <- match(c("es", "var","mod", "data","method"),
    names(mf), 0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf.es <- mf[[match("es", names(mf))]]
    meta <- data
    meta$es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
    mf.var.es <- mf[[match("var", names(mf))]]
    meta$var.es <- eval(mf.var.es, data, enclos = sys.frame(sys.parent()))
    mf.mod <- mf[[match("mod", names(mf))]]
    meta$mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
    meta$mod <- as.character(meta$mod)
      meta$k <- 1
      meta$TW <- 1/meta$var.es
      meta$TWD <- meta$TW * meta$es
      meta$TWDS <- meta$TWD * meta$es
      temp <- aggregate(meta[c("k", "TW", "TWD", "TWDS")], by = meta["mod"], 
          FUN = sum)
      temp$mod <- as.character(temp$mod)
      lastrows <- dim(temp)[1] + 1
      temp[lastrows, 2:5] <- apply(temp[, -1], 2, FUN = sum)
      temp$mod[lastrows] <- "Overall"
      temp$es <- temp$TWD/temp$TW
      temp$var.es <- 1/temp$TW
      temp$se.es <- sqrt(temp$var.es)
      temp$Q <- temp$TWDS - temp$TWD^2/temp$TW
      temp$df <- temp$k - 1
      temp$z.value <- temp$es/temp$se.es
      temp$p.value <- 2 * pnorm(abs(temp$z.value), lower.tail = FALSE)
      temp$p_homog <- ifelse(temp$df == 0, 1, pchisq(temp$Q, temp$df, 
          lower = FALSE))
      temp$I2 <- (temp$Q - (temp$df))/temp$Q
      temp$I2 <- ifelse(temp$I2 < 0, 0, temp$I2)
      temp$I2 <- paste(round(temp$I2 * 100, 0), "%", sep = "")
      # fixed effect homogeneity (used in RE for Q, Qw, df.w, p.w)
      kf <- temp$k[temp$mod == "Overall"]
      levelsf <- length(temp$mod) - 1
      Qb.dff <- levelsf - 1
      Qw.dff <- kf - levelsf
      Qsf <- temp$Q[temp$mod == "Overall"]
      Qwf <- sum(temp$Q[!temp$mod == "Overall"])
      Qw_p.valuef <- 1 - pchisq(Qwf, Qw.dff)
    if(method == "fixed") {
      out <- aggregate(meta[c("k", "TW", "TWD", "TWDS")], by = meta["mod"], 
          FUN = sum)
      out$mod <- as.character(out$mod)
      lastrow <- dim(out)[1] + 1
      out[lastrow, 2:5] <- apply(out[, -1], 2, FUN = sum)
      out$mod[lastrow] <- "Overall"
      out$es <- out$TWD/out$TW
      out$var.es <- 1/out$TW
      out$se.es <- sqrt(out$var.es)
      out$Q <- out$TWDS - out$TWD^2/out$TW
      out$df <- out$k - 1
      out$z.value <- out$es/out$se.es
      out$p.value <- 2 * pnorm(abs(out$z.value), lower.tail = FALSE)
      out$L.95ci <- out$es - 1.96 * out$se.es
      out$U.95ci <- out$es + 1.96 * out$se.es
      out$p_homog <- ifelse(out$df == 0, 1, pchisq(out$Q, out$df, 
          lower = FALSE))
      out$I2 <- (out$Q - (out$df))/out$Q
      out$I2 <- ifelse(out$I2 < 0, 0, out$I2)
      out$I2 <- paste(round(out$I2 * 100, 0), "%", sep = "")
      out$TW <- NULL
      out$TWD <- NULL
      out$TWDS <- NULL      
      k <- out$k[out$mod == "Overall"]
      levels <- length(out$mod) - 1
      Qb.df <- levels - 1
      Qw.df <- k - levels
      Qs <- out$Q[out$mod == "Overall"]
      Qw <- sum(out$Q[!out$mod == "Overall"])
      Qw_p.value <- 1 - pchisq(Qw, Qw.df)
      Qb <- Qs - Qw
      Qb_p.value <- 1 - pchisq(Qb, Qb.df)
      mod.Qstat <- data.frame(Qs, Qw, Qw.df, Qw_p.value, Qb, Qb.df, 
        Qb_p.value)
      names(mod.Qstat) <- c("Q", "Qw", "df.w", "p.w", "Qb", "df.b", 
        "p.b")
    }
    if(method=="random") {
      meta$TWS <- meta$TW^2
      out <- aggregate(meta[c("k", "TW", "TWS", "TWD", "TWDS")], 
          by = meta["mod"], FUN = sum)
      out$Q <- out$TWDS - out$TWD^2/out$TW
      Q <- sum(out$Q)
      out$df <- out$k - 1
      df <- sum(out$df)
      out$c <- out$TW - (out$TWS/out$TW)
      c <- sum(out$c)
      TSw <- (Q - df)/c
      tau <- ifelse(TSw < 0, 0, TSw)
      meta$var.tau <- meta$var.es + tau
      meta$TW.tau <- 1/meta$var.tau
      meta$TWD.tau <- meta$TW.tau * meta$es
      meta$TWDS.tau <- meta$TWD.tau * meta$es
      out <- aggregate(meta[c("k", "TW.tau", "TWD.tau", "TWDS.tau")], 
          by = meta["mod"], FUN = sum)
      lastrow <- dim(out)[1] + 1
      out[lastrow, 2:5] <- apply(out[, -1], 2, FUN = sum)
      out$mod[lastrow] <- "Overall"
      out$es <- out$TWD.tau/out$TW.tau
      out$var.es <- 1/out$TW.tau
      out$es <- out$TWD.tau/out$TW.tau
      out$se.es <- sqrt(out$var.es)
      out$Q <-  out$TWDS.tau - out$TWD.tau^2/out$TW.tau
      out$df <- out$k - 1
      out$z.value <- out$es/out$se.es
      out$p.value <- 2 * pnorm(abs(out$z.value), lower.tail = FALSE)
      out$L.95ci <- out$es - 1.96 * out$se.es
      out$U.95ci <- out$es + 1.96 * out$se.es
      out$p_homog <-  ifelse(out$df == 0, 1, pchisq(out$Q, out$df, 
        lower = FALSE))
      out$I2 <- temp$I2
      out$TW.tau <- NULL
      out$TWD.tau <- NULL
      out$TWDS.tau <- NULL 
      k <- out$k[out$mod == "Overall"]
      levels <- length(out$mod) - 1
      Qb.df <- levels - 1
      Qw.df <- k - levels
      Q <- out$Q[out$mod == "Overall"]
      Qw <- sum(out$Q[!out$mod == "Overall"])
      Qw_p.value <- 1 - pchisq(Qw, Qw.df)
      Qb <- Q - Qw
      Qb_p.value <- 1 - pchisq(Qb, Qb.df)
      Qs <- temp$Q[temp$mod == "Overall"]
      Qw <- sum(temp$Q[!temp$mod == "Overall"])
      Qw_p.value <- Qw_p.valuef
      mod.Qstat <- data.frame(Qs, Qw, Qw.df, Qw_p.value, Qb, Qb.df, 
        Qb_p.value)
      names(mod.Qstat) <- c("Q", "Qw", "df.w", "p.w", "Qb", "df.b", 
        "p.b")
      out$Q <- temp$Q
      out$df <- temp$df
      out$p_homog <- temp$p_homog 
    }
    colnames(out) <- c("mod", "k", "estimate", "var", "se", "Q", 
          "df", "z", "p", "ci.l", "ci.u", "p.h", "I2")
    out1 <- out[, c(1, 2, 3, 5, 4, 10, 11, 8, 9, 6, 7, 12, 13)]
    if(ztor == TRUE) {
      out1$estimate <- (exp(2*out$estimate)-1)/(exp(2*out$estimate) + 1)
      out1$var <- (exp(2*out$var)-1)/(exp(2*out$var) + 1)
      out1$se <- sqrt(out$var)
      out1$z <- (exp(2*out$z)-1)/(exp(2*out$z) + 1)
      out1$p <- (exp(2*out$p)-1)/(exp(2*out$p) + 1)
      out1$ci.l <- (exp(2*out$ci.l)-1)/(exp(2*out$ci.l) + 1)
      out1$ci.u <- (exp(2*out$ci.u)-1)/(exp(2*out$ci.u) + 1)
    } 
    out2 <- mod.Qstat
    out <- list(Model=out1, Heterogeneity=out2)
    class(out) <- "macat"
    return(out)
}

print.macat <- function (x, digits = 4, ...){
    x1 <- x$Model
    x2 <- x$Heterogeneity
    mod <- x1$mod
    k <- x1$k
    estimate <- x1$estimate
    var <- x1$var
    se <- x1$se
    ci.l <- x1$ci.l
    ci.u <- x1$ci.u 
    z <- x1$z
    p <- x1$p
    Q <- x1$Q
    df <- x1$df
    p.h <- x1$p.h
    I2 <- x1$I2
    Qoverall <- x2$Q
    Qw <- x2$Qw
    Qw.df <- x2$df.w
    Qw.p <- x2$p.w
    Qb <- x2$Qb
    Qb.df <- x2$df.b
    Qb.p <- x2$p.b
    results1 <- round(data.frame(estimate, var, se, ci.l, ci.u, z, p,
      Q, df, p.h),4)
    #results <- formatC(table, format="f", digits=digits)
    results1$k <- x1$k
    results1$I2 <- x1$I2
    results1$mod <- x1$mod
    results1 <- results1[c(13,11, 1:10,12)]
    results2 <- round(data.frame(Q=Qoverall, Qw, Qw.df, Qw.p, Qb, Qb.df, Qb.p),4)
    cat("\n Note: 
       If using r and the variance of r as input, be sure to 
       leave the default ztor = FALSE, otherwise the output will be
       inaccurate. If using z' (Fisher's z) and the variance of z'
       as the input, changing ztor = TRUE will convert z' back to r 
       (for interpretive purposes) after all analyses have been 
       conducted.", "", "\n","----","\n")
    cat("\n Model Results:", "", "\n", "\n")
    print(results1)
    cat("\n Heterogeneity:", "", "\n","\n")
    print(results2)
    invisible(x)
}


      
# Function for planned comparisons between 2 levels of moderator (fixed effects)

macatC <- function(x1, x2, es, var, mod, data,  method= "random", type= "post.hoc",
          ztor = FALSE) {
  # Directly compares 2 levels of a categorical moderator using a fixed effects model. 
  # Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with g (standardized mean diff),var.g (variance of g),
  #         n.1 (grp 1 sample size), n.2 (grp 2 sample size).
  #   mod: Categorical moderator variable used for moderator analysis.    
  #   x1: One level of categorical moderator
  #   x2: Other level (comparison group) of same categorical moderator
  #   method: "post.hoc" assumes the comparision was not planned prior to conducting
  #           the meta-analysis. The other option "planned" assumes you have planned 
  #           a priori to compare these levels of the categorical moderator. 
  #           Default is "post.hoc1". 
  # Returns:
  #   Random effects moderator means per group, mean difference (d), variance of difference,
  #   p-value, and 95% confidence intervals. 
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("x1", "x.2", "es", "var", "mod", "data"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.es <- mf[[match("es", names(mf))]]
  es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var <- mf[[match("var", names(mf))]]
  var.es <- eval(mf.var, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  modsig <- macat(es=es, var=var.es, mod=mod, data=data, method=method)
  modsig$mod <- as.factor(modsig$Model$mod)
  com1 <- levels(modsig$mod)[x1]  # first level of moderator
  com2 <- levels(modsig$mod)[x2]  # second level of moderator
  x1.es <- modsig$Model[modsig$mod==com1, "estimate"]
  x2.es <- modsig$Model[modsig$mod==com2, "estimate"]
  x1.var <- modsig$Model[modsig$mod==com1, "var"]
  x2.var <- modsig$Model[modsig$mod==com2, "var"]
  es <- (-1)*x1.es + 1*x2.es  # pg 288 (Cooper et al.,  2009)
  var <- (-1)^2*x1.var + (1)^2*x2.var
  df <- 2-1
  chi.sqr <- es^2/var
  #m <- meta 
  #m$mod <- mod
  #compl<-!is.na(m$mod)
  #meta<-m[compl,]
  #if(type == "post.hoc1") {  # post-hoc comparison (Tukey HSD method)
  #  fit<-aov(es~mod,weights=meta$wi)
  #  fit<-TukeyHSD(fit)
  #}
  if(type == "post.hoc") {  # post-hoc comparison (Scheffe method)
    z <- es/sqrt(var)
    z2 <- z^2
    levels <- length(levels(as.factor(mod)))
    df.post <- levels-1 
    p.value <- 1-pchisq(z2, df.post)
    L.95ci <- es-1.96*sqrt(var)
    U.95ci <- es + 1.96*sqrt(var)
    fit <- data.frame(es, var, p.value) # , L.95ci, U.95ci)
    names(fit) <- c("diff", "var.diff",  
                    "p") #, "CI.lower", "CI.lower" )
  }
  if (type == "planned") {  # planned comparison (a priori)
   df <- 2-1
   chi.sqr <- es^2/var
   p.value <- 1-pchisq(chi.sqr, df)
   L.95ci <- es-1.96*sqrt(var)
   U.95ci <- es + 1.96*sqrt(var)
   fit <- data.frame(es, var, p.value) # , L.95ci, U.95ci)
   names(fit) <- c("diff", "var.diff",  
                   "p") #, "ci.l", "ci.u" )
  }
  if(ztor == TRUE) {
      fit$diff <- (exp(2*fit$diff)-1)/(exp(2*fit$diff) + 1)
      fit$var.diff <- (exp(2*fit$var.diff)-1)/(exp(2*fit$var.diff) + 1)
      fit$p <- (exp(2*fit$p)-1)/(exp(2*fit$p) + 1)
      #fit$ci.l <- (exp(2*fit$ci.l)-1)/(exp(2*fit$ci.l) + 1)
      #fit$ci.u <- (exp(2*fit$ci.u)-1)/(exp(2*fit$ci.u) + 1)
    } 
  return(fit) 
}

# multifactor cat mod analysis [IN PROGRESS]:
#add relevant columns and it will works fine!
#MFCatMod <- function(meta, mod1, mod2) {
#  m <- Wifun(meta)
#  m$mod1 <- mod1
#  m$mod2 <- mod2
#  fixed <- ddply(m, c("mod1", "mod2"), summarise, sum.wi = sum(wi),
#           sum.wiTi = sum(wiTi), sum.wiTi2 = sum(wiTi2))
#  fixed$ES <- fixed$sum.wiTi/fixed$sum.wi
#  random <- ddply(m, c("mod1", "mod2"), summarise, sum.wi.tau = sum(wi.tau),
#           sum.wiTi.tau = sum(wiTi.tau), sum.wiTi2.tau = sum(wiTi2.tau))
#  random$ES <- random$sum.wiTi.tau/random$sum.wi.tau
#  out <- list(Fixed = fixed, Random = random)
#  return(out)
#}


##============= GRAPHICS =============##

# requires ggplot2

##=== Meta-regression scatterplot with weighted regression line ===##

plotcon <- function(es, var, mod, data, method= "random", modname=NULL, 
  title=NULL, ylim=c(0, 1), ...) {
  # Outputs a scatterplot from a fixed or random effects meta-regression (continuous and/or
  # categorical). Computations derived from chapter 14 and 15, Cooper et al. (2009). 
  # Args:
  #   data: data.frame with es (r or z'),var (variance),
  #         n.1 (grp 1 sample size), n.2 (grp 2 sample size).
  #   mod: Moderator variable used for meta-regression.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   modname: Name of moderator to appear on x-axis of plot. Default is NULL.
  #   title: Plot title. Default is NULL.
  #   ylim: Limits of y-axis with the first argrument minimal value and second maximum value.
  #         Default is c(0,1).
  # Returns:
  #   Scatterplot with fixed or random effects regression line and where size of points are 
  #   based on study weights--more precise studies are larger. The ggplot2 package outputs the 
  #   rich graphics. 
  require('ggplot2')
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("es", "var","mod", "data","method"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.es <- mf[[match("es", names(mf))]]
  m <- data
  m$es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var.es <- mf[[match("var", names(mf))]]
  m$var.es <- eval(mf.var.es, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  m$mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  m$mod <- as.character(m$mod)
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  meta$wi <-  1/meta$var.es
  meta$wi2 <- (1/meta$var.es)^2
  meta$wiTi <- meta$wi*meta$es
  meta$wiTi2 <- meta$wi*(meta$es)^2
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$es))                                  
  df <- k-1      
  tau <- (Q-(k - 1))/comp 				
  meta$var.tau <- meta$var.es + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$es
  meta$wiTi2.tau <- meta$wi.tau*(meta$es)^2
  if(method=="fixed") {
    congraph <- ggplot(meta,  aes(mod, es, weight=wi), na.rm=TRUE, ...) + 
    geom_point(aes(size=wi), alpha=.6, na.rm=TRUE) + 
    geom_smooth(aes(group=1), method = lm,  se = FALSE) + 
    xlab(modname) + ylab("Effect Size") +  
    ylim(ylim) +  
    opts(title=title, legend.position = "none")
  }
  if(method=="random") {
    congraph <- ggplot(meta,  aes(mod, es, weight=wi.tau), na.rm=TRUE) + 
    geom_point(aes(size=wi.tau), alpha=.6, na.rm=TRUE) + 
    geom_smooth(aes(group=1, weight=wi.tau), method = lm, se = FALSE,  na.rm=TRUE) + 
    xlab(modname) + 
    ylab("Effect Size")  + 
    #ylim(min(meta$z), 1) +  
    opts(title=title, legend.position = "none")
  }
  return(congraph)
}

##=== Categorical Moderator Graph ===##

# Intermediate level function to add mean to boxplot

stat_sum_single1  <-  function(fun,  geom="point",  weight=wi, ...) {
                               stat_summary(fun.y=fun,  shape=" + ",
                               geom=geom,  size = 5,  ...)      
}
stat_sum_single2  <-  function(fun,  geom="point",  weight=wi.tau, ...) {
                               stat_summary(fun.y=fun,  shape=" + ",  
                               geom=geom,  size = 5,  ...)      
}

plotcat <- function(es, var, mod, data,  method="random",  modname=NULL,  title=NULL, ...) {
  # Outputs a boxplot from a fixed or random effects moderator analysis.
  # Computations derived from chapter 14 and 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with es (r or z'),var.es (variance of es),
  #         n.1 (grp 1 sample size), n.2 (grp 2 sample size).
  #   mod: Categorical moderator variable used for analysis.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   modname: Name of moderator to appear on x-axis of plot. Default is NULL.
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Boxplot graph with median, interquartile range, max, min, and 
  #   outliers from a fixed or random effects categorical moderator analysis. Places
  #   jitter points for each study and the size of points are based on study weights--more 
  #   precise studies are larger. The ggplot2 package outputs the 
  #   rich graphics. 
  require('ggplot2')
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("es", "var","mod", "data","method"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.es <- mf[[match("es", names(mf))]]
  m <- data
  m$es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var.es <- mf[[match("var", names(mf))]]
  m$var.es <- eval(mf.var.es, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  m$mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  m$mod <- as.character(m$mod)
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  meta$wi <-  1/meta$var.es
  meta$wi2 <- (1/meta$var.es)^2
  meta$wiTi <- meta$wi*meta$es
  meta$wiTi2 <- meta$wi*(meta$es)^2
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$es))                                  
  df <- k-1      
  tau <- (Q-(k - 1))/comp 				
  meta$var.tau <- meta$var.es + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$es
  meta$wiTi2.tau <- meta$wi.tau*(meta$es)^2
  if(method=="fixed") {
    catmod <- ggplot(meta,  aes(factor(mod), es,weight = wi), na.rm=TRUE) + 
                    geom_boxplot(aes(weight = wi), outlier.size=2,na.rm=TRUE) + 
                    geom_jitter(aes(shape=factor(mod), size=wi.tau), alpha=.25) + 
                    #theme_bw() + 
                    xlab(modname) + 
                    ylab("Effect Size")  + 
                    opts(title=title, legend.position="none", na.rm=TRUE) #+ 
                    #stat_sum_single2(mean)
  }  
  if(method=="random") {
    catmod <- ggplot(meta,  aes(factor(mod), es,  weight=wi.tau), na.rm=TRUE) + 
                    geom_boxplot(outlier.size=2, aes(weight=wi.tau),na.rm=TRUE) + 
                    geom_jitter(aes(shape=factor(mod), size=wi.tau), alpha=.25) + 
                    #theme_bw() + 
                    xlab(modname) + 
                    ylab("Effect Size")  + 
                    opts(title=title, legend.position="none", na.rm=TRUE) # + 
                    #stat_sum_single2(mean)
  }
  return(catmod)
}

##===================== PUBLICATION BIAS ===================##
# Three approaches to assess for publication bias: (1) Fail
# Safe N, (2) Trim & Fill, and (3) Selection Modeling.

# Fail Safe N provides an estimate of the number of missing studies that 
# would need to exist to overturn the current conclusions.

PubBias <- function(data) {  # requires a data.frame having been analyzed by 'weights' 
  meta <- data               # function
  k <- length(meta$z.score)
  sum.z <- sum(meta$z.score)
  Z <- sum.z/sqrt(k)
  k0 <- round(-k + sum.z^2/(1.96)^2, 0)
  k.per <- round(k0/k, 0)
  out <- list("Fail.safe"=cbind(k,Z,k0, k.per))
  return(out)
}

##================== INTERRATER RELIABILITY ================##

# Kappa coefficients for inter-rater reliability (categorical variables)
# Imputs required are rater1 (first rater on Xi categorical variable)
# and rater2 (second rater on same Xi categorical variable)

Kappa <- function(rater1, rater2)  {
  # Computes Kappa coefficients for inter-rater reliability (categorical variables).
  # Args:
  #   rater1: First rater of categorical variable to be analyzed.
  #   rater2: Second rater on same categorical variable to be analyzed.
  # Returns:
  #   Kappa coefficients for inter-rater reliability (categorical variables).
  freq <- table(rater1, rater2)  # frequency table
  marg <- margin.table(freq) # total observations 
  marg2 <- margin.table(freq, 1)  # A frequencies (summed over rater2) 
  marg1 <- margin.table(freq, 2)  # B frequencies (summed over rater1)
  cellper <- prop.table(freq)  # cell percentages
  rowper <- margin.table(freq, 2)/margin.table(freq)  # row percentages 
  colper <- margin.table(freq, 1)/margin.table(freq)  # column percentages
  expected  <-  as.array(rowper) %*% t(as.array(colper)) 
  p.e <- sum(diag(expected))
  p.a_p.e <- sum(diag(cellper))- p.e
  p.e1 <- 1-p.e
  kappa <- p.a_p.e/p.e1
  names(kappa) <- "Kappa"
  #cat("\n")
  #cat("strong agreement:    kappa > .75", "", "\n")
  #cat("moderate agreement:  .40 > kappa < .75", "", "\n")
  #cat("weak agreement:      kappa < .40", "", "\n","\n")
  return(kappa)
}

# icc (created by Matthias Gamer for the 'irr' package)
icc <- function (ratings, model = c("oneway", "twoway"), type = c("consistency", 
    "agreement"), unit = c("single", "average"), r0 = 0, conf.level = 0.95) 
{
    ratings <- as.matrix(na.omit(ratings))
    model <- match.arg(model)
    type <- match.arg(type)
    unit <- match.arg(unit)
    alpha <- 1 - conf.level
    ns <- nrow(ratings)
    nr <- ncol(ratings)
    SStotal <- var(as.numeric(ratings)) * (ns * nr - 1)
    MSr <- var(apply(ratings, 1, mean)) * nr
    MSw <- sum(apply(ratings, 1, var)/ns)
    MSc <- var(apply(ratings, 2, mean)) * ns
    MSe <- (SStotal - MSr * (ns - 1) - MSc * (nr - 1))/((ns - 
        1) * (nr - 1))
    if (unit == "single") {
        if (model == "oneway") {
            icc.name <- "ICC(1)"
            coeff <- (MSr - MSw)/(MSr + (nr - 1) * MSw)
            Fvalue <- MSr/MSw * ((1 - r0)/(1 + (nr - 1) * r0))
            df1 <- ns - 1
            df2 <- ns * (nr - 1)
            p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
            FL <- (MSr/MSw)/qf(1 - alpha/2, ns - 1, ns * (nr - 
                1))
            FU <- (MSr/MSw) * qf(1 - alpha/2, ns * (nr - 1), 
                ns - 1)
            lbound <- (FL - 1)/(FL + (nr - 1))
            ubound <- (FU - 1)/(FU + (nr - 1))
        }
        else if (model == "twoway") {
            if (type == "consistency") {
                icc.name <- "ICC(C,1)"
                coeff <- (MSr - MSe)/(MSr + (nr - 1) * MSe)
                Fvalue <- MSr/MSe * ((1 - r0)/(1 + (nr - 1) * 
                  r0))
                df1 <- ns - 1
                df2 <- (ns - 1) * (nr - 1)
                p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
                FL <- (MSr/MSe)/qf(1 - alpha/2, ns - 1, (ns - 
                  1) * (nr - 1))
                FU <- (MSr/MSe) * qf(1 - alpha/2, (ns - 1) * 
                  (nr - 1), ns - 1)
                lbound <- (FL - 1)/(FL + (nr - 1))
                ubound <- (FU - 1)/(FU + (nr - 1))
            }
            else if (type == "agreement") {
                icc.name <- "ICC(A,1)"
                coeff <- (MSr - MSe)/(MSr + (nr - 1) * MSe + 
                  (nr/ns) * (MSc - MSe))
                a <- (nr * r0)/(ns * (1 - r0))
                b <- 1 + (nr * r0 * (ns - 1))/(ns * (1 - r0))
                Fvalue <- MSr/(a * MSc + b * MSe)
                a <- (nr * coeff)/(ns * (1 - coeff))
                b <- 1 + (nr * coeff * (ns - 1))/(ns * (1 - coeff))
                v <- (a * MSc + b * MSe)^2/((a * MSc)^2/(nr - 
                  1) + (b * MSe)^2/((ns - 1) * (nr - 1)))
                df1 <- ns - 1
                df2 <- v
                p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
                FL <- qf(1 - alpha/2, ns - 1, v)
                FU <- qf(1 - alpha/2, v, ns - 1)
                lbound <- (ns * (MSr - FL * MSe))/(FL * (nr * 
                  MSc + (nr * ns - nr - ns) * MSe) + ns * MSr)
                ubound <- (ns * (FU * MSr - MSe))/(nr * MSc + 
                  (nr * ns - nr - ns) * MSe + ns * FU * MSr)
            }
        }
    }
    else if (unit == "average") {
        if (model == "oneway") {
            icc.name <- paste("ICC(", nr, ")", sep = "")
            coeff <- (MSr - MSw)/MSr
            Fvalue <- MSr/MSw * (1 - r0)
            df1 <- ns - 1
            df2 <- ns * (nr - 1)
            p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
            FL <- (MSr/MSw)/qf(1 - alpha/2, ns - 1, ns * (nr - 
                1))
            FU <- (MSr/MSw) * qf(1 - alpha/2, ns * (nr - 1), 
                ns - 1)
            lbound <- 1 - 1/FL
            ubound <- 1 - 1/FU
        }
        else if (model == "twoway") {
            if (type == "consistency") {
                icc.name <- paste("ICC(C,", nr, ")", sep = "")
                coeff <- (MSr - MSe)/MSr
                Fvalue <- MSr/MSe * (1 - r0)
                df1 <- ns - 1
                df2 <- (ns - 1) * (nr - 1)
                p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
                FL <- (MSr/MSe)/qf(1 - alpha/2, ns - 1, (ns - 
                  1) * (nr - 1))
                FU <- (MSr/MSe) * qf(1 - alpha/2, (ns - 1) * 
                  (nr - 1), ns - 1)
                lbound <- 1 - 1/FL
                ubound <- 1 - 1/FU
            }
            else if (type == "agreement") {
                icc.name <- paste("ICC(A,", nr, ")", sep = "")
                coeff <- (MSr - MSe)/(MSr + (MSc - MSe)/ns)
                a <- r0/(ns * (1 - r0))
                b <- 1 + (r0 * (ns - 1))/(ns * (1 - r0))
                Fvalue <- MSr/(a * MSc + b * MSe)
                a <- (nr * coeff)/(ns * (1 - coeff))
                b <- 1 + (nr * coeff * (ns - 1))/(ns * (1 - coeff))
                v <- (a * MSc + b * MSe)^2/((a * MSc)^2/(nr - 
                  1) + (b * MSe)^2/((ns - 1) * (nr - 1)))
                df1 <- ns - 1
                df2 <- v
                p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
                FL <- qf(1 - alpha/2, ns - 1, v)
                FU <- qf(1 - alpha/2, v, ns - 1)
                lbound <- (ns * (MSr - FL * MSe))/(FL * (MSc - 
                  MSe) + ns * MSr)
                ubound <- (ns * (FU * MSr - MSe))/(MSc - MSe + 
                  ns * FU * MSr)
            }
        }
    }
    rval <- structure(list(subjects = ns, raters = nr, model = model, 
        type = type, unit = unit, icc.name = icc.name, value = coeff, 
        r0 = r0, Fvalue = Fvalue, df1 = df1, df2 = df2, p.value = p.value, 
        conf.level = conf.level, lbound = lbound, ubound = ubound), 
        class = "icclist")
    return(rval)
}

print.icclist <- function (x, ...) 
{
    icc.title <- ifelse(x$unit == "single", "Single Score Intraclass Correlation", 
        "Average Score Intraclass Correlation")
    cat(paste(" ", icc.title, "\n\n", sep = ""))
    cat(paste("   Model:", x$model, "\n"))
    cat(paste("   Type :", x$type, "\n\n"))
    cat(paste("   Subjects =", x$subjects, "\n"))
    cat(paste("     Raters =", x$raters, "\n"))
    results <- paste(formatC(x$icc.name, width = 11, flag = "+"), 
        "=", format(x$value, digits = 3))
    cat(results)
    cat("\n\n F-Test, H0: r0 =", x$r0, "; H1: r0 >", x$r0, "\n")
    Ftest <- paste(formatC(paste("F(", x$df1, ",", format(x$df2, 
        digits = 3), ")", sep = ""), width = 11, flag = "+"), 
        "=", format(x$Fvalue, digits = 3), ", p =", format(x$p.value, 
            digits = 3), "\n\n")
    cat(Ftest)
    cat(" ", round(x$conf.level * 100, digits = 1), "%-Confidence Interval for ICC Population Values:\n", 
        sep = "")
    cat(paste("  ", round(x$lbound, digits = 3), " < ICC < ", 
        round(x$ubound, digits = 3), "\n", sep = ""))
}


# Correction for Attenuation

atten <- function (es, xx, yy, data) {
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    args <- match(c("es", "xx", "yy", "data"),
    names(mf), 0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    meta <- data 
    mf.es <- mf[[match("es", names(mf))]]
    es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
    mf.xx <- mf[[match("xx", names(mf))]]
    xx <- eval(mf.xx, data, enclos = sys.frame(sys.parent()))
    mf.yy <- mf[[match("yy", names(mf))]]
    yy <- eval(mf.yy, data, enclos = sys.frame(sys.parent()))
    meta$es.corrected <-  es/(sqrt(xx)*sqrt(yy))
    return(meta)
}








##---------------- old functions -------------------
# 3.18.10  internal for GUI: convert to factor

facts <- function(meta, mod) {
  meta[,mod] <- factor(meta[,mod])
  return(meta)
}


# 2.24.10: added (1) conversion to d functions, (2) CatComp
# update
# CatMod update
# new element of function for mods 2.20.10. Big thanks to 
# Jim Holtman for the function to randomly select mods
# integrated into agg_r2()
# New function: 2.16.10
# Multiple predictor regression model fit
# (required inputs are lm saved as objects)

MRfit <- function( ...) {
  models <- list(...)
  fit<- do.call(anova, models)
  fit$R2 <- 1 - (fit$RSS / fit$RSS[1])
  fit$R2.change <- c(NA, diff(fit$R2))
  fit$R2[1] <- NA
  fit$F <- NULL
  names(fit) <- c("df.Q", "Qe", "predictors", "Q", "p", 
                  "R^2", "R^2.change")
  fit$p <- NULL
  return(fit)
}



# a formula for attenuation already created:  correct.cor(x, y)


##=== Preliminary Steps ===##

# Import data into R:
# 1. Save main data file (excel or spss) to .csv [e.g.,  see save options in excel]
# 2. Import .csv file into R by setting the working directory to the location of 
#    your data file,  e.g.:
#    setwd("C:/Users/User/Documents/TA Meta-Analy/Horvath_2009/ANALYSIS/12-10-09") 
#    and then import data,  e.g.:
#    data <- read.csv("Alliance_1-30-10.csv", header=TRUE, na.strings="") 

##==== Data Manipulation ====##


# set numeric variables to numeric, e.g.:
# data$r <- as.numeric(as.character(data$r))

# set categorical variables to factors or character,  e.g.:
# data$id <- as.character(data$id)

# fix data with errors in factor names, requires car package,e.g.:
# library(car)
# data$outcome3 <- recode(data$outcome2, 'c("?",  "adherence", "compliance", "depression", 
#                         "depression ", "wellbeing", "work", "GAS")="Other"; 
#                         c("GSI", "SCL", "BSI")="SCL"; c("dropout")="Dropout"; 
#                         else= "Other"')    

##============ COMPUTATIONS TO CALCULATE EFFECT SIZES ================##

# Formulas for computing d, var(d), g (bias removed from d), and var(g)
# in designs with independent groups.
# Section 12.3.1 & Table 12.1 (Cooper et al., 2009; pp. 226-228)

# Computing d and g, independent groups
# (12.3.1) Study reported: 
# m.1 (post-test mean of treatment), m.2 (post-test mean of comparison),
# sd.1 (treatment standard deviation at post-test), sd.2 (comparison 
# standard deviation at post-test), n.1 (treatment), n.2 (comparison/control).

d_to_g <- function(d, var.d, n.1, n.2) {
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  out<-cbind(g, var.g)
  return(out)
}

mean_to_d <- function(m.1,m.2,sd.1,sd.2,n.1, n.2) {
  s.within<-sqrt(((n.1-1)*sd.1^2+(n.2-1)*sd.2^2)/(n.1+n.2-2))
  d<-(m.1-m.2)/s.within
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var.d)
  return(out)
}

# (1) Study reported: 
# m.1 (post-test mean of treatment), m.2 (post-test mean of comparison),
# s.pooled (pooled standard deviation), n.1 (treatment), 
# n.2 (comparison/control).

mean_to_d2 <- function(m.1,m.2,s.pooled,n.1, n.2) {
  d<-(m.1-m.2)/s.pooled
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (2) Study reported: 
# t (t-test value of treatment v comparison), n.1 (treatment),
# n.2 (comparison/control).

t_to_d <- function(t, n.1, n.2) {
  d<-t*sqrt((n.1+n.2)/(n.1*n.2))
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (3) Study reported: 
# f (F-test value of treatment v comparison), n.1 (treatment),
# n.2 (comparison/control).


f_to_d <- function(f,n.1, n.2) {
  d<-sqrt(f*(n.1+n.2)/(n.1*n.2))
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (4) Study reported: 
# p-value (for ONE-tailed test), n.1 (treatment), n.2 (comparison/control).

p_to_d1 <- function(p, n.1, n.2) {
  pxtwo<-p*2
  df<-(n.1+n.2)-2  
  TINV<-qt((1-pxtwo/2),df)
  d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (5) Study reported: 
# p-value (for TWO-tailed test), n.1 (treatment), n.2 (comparison/control).

p_to_d2 <- function(p, n.1, n.2) {
  df<-(n.1+n.2)-2  
  TINV<-qt((1-p/2),df)
  d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# rtod converts Pearson r to Cohen's d (also retains (N-1)/N).
r_to_d <- function(r, N) { 2*r*sqrt((N-1)/(N*(1-r^2)))*abs(r)/r    }

# Formulas for computing d and var(d) in designs with independent groups
# using ANCOVA. Section 12.3.3 & Table 12.3 (Cooper et al., 2009; pp. 228-230).

# Computing d and g from ANCOVA
# (12.3.3) Study reported: 
# m.1.adj (adjusted mean of treatment from ANCOVA),
# m.2.adj (adjusted mean of comparison/control from ANCOVA),
# sd.adj (adjusted standard deviation), n.1 (treatment),
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

ancova_to_d1 <- function(m.1.adj,m.2.adj,sd.adj,n.1, n.2, R, q) {
  s.within<-sd.adj/sqrt(1-R^2)
  d<-(m.1.adj-m.2.adj)/s.within
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}


# Table 12.3 (1) Study reported: 
# m.1.adj (adjusted mean of treatment from ANCOVA),
# m.2.adj (adjusted mean of comparison/control from ANCOVA),
# s.pooled (pooled standard deviation), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

ancova_to_d2 <- function(m.1.adj, m.2.adj, s.pooled, n.1, n.2, R, q) {
  d<-(m.1.adj-m.2.adj)/s.pooled
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (2) Study reported: 
# t (t-test value from ANCOVA), n.1 (treatment),
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

tt.ancova_to_d <- function(t, n.1, n.2, R, q) {
  d<-t*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (3) Study reported: 
# f (F-test value from ANCOVA), n.1 (treatment),
# n.2 (comparison/control),R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

f.ancova_to_d<-function(f,n.1, n.2, R, q) {
  d<-sqrt(f*(n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (4) Study reported: 
# p-value (for ONE-tailed test, from ANCOVA), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

p.ancova_to_d1 <- function(p, n.1, n.2, R, q) {
  pxtwo<-p*2
  df<-(n.1+n.2)-2  
  TINV<-qt((1-pxtwo/2),df)
  d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (5) Study reported: 
# p-value (for TWO-tailed test, from ANCOVA), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

p.ancova_to_d2 <- function(p, n.1, n.2, R, q) {
  df<-(n.1+n.2)-2  
  TINV<-qt((1-p/2),df)
  d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# computing d from odds ratio

or_to_d <- function(or) {
  lor <- log(or)
  d <- lor*sqrt(3)/pi
  out <- cbind(lor, d)
  return(out)
}  

# computing d from log odds ratio

lor_to_d <- function(lor, var.lor) {
  d <- lor*sqrt(3)/pi
  var.d <- 3*var.lor/pi^2
  out <- cbind(d, var.d)
  return(out)
}  

# compute or from proportions

#prop_to_or <- function(p1, p2, n.ab, n.cd) {
#  or <-(p1*(1-p2))/(p2*(1-p1))
#  lor <- log(or)
#  var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
#  out <- cbind(or,lor,var.lor)
#  return(out)
#}

#prop_to_d <-function(p1, p2, n.ab, n.cd) {
#  or <-(p1*(1-p2))/(p2*(1-p1))
#  lor <- log(or)
#  var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
#  d <- lor*sqrt(3)/pi
#  var.d <- 3*var.lor/pi^2
#  out <- cbind(or,lor,var.lor, d, var.d)
#  return(out)
#}

# Odds Ratio to d: if have info for 'failure' in both conditions 
# (B = # tmt failure; D = # non-tmt failure) and the sample size
# for each group (n.1 & n.2 respectively):

fail_to_d <- function(B, D, n.1, n.0) {
  A <- n.1 - B  # tmt success
  B <- B        # tmt failure
  C <- n.0 - D  # non-tmt success
  D <- D        # non-tmt failure
  p1 <- A/n.1   # proportion 1 
  p2 <- C/n.0   # proportion 2
  n.ab <-  A+B  # n of A+B
  n.cd <-  C+D  # n of C+D        
  or <- (p1 * (1 - p2))/(p2 * (1 - p1))  # odds ratio
  lor <- log(or)  # log odds ratio
  var.lor <-  1/A + 1/B + 1/C + 1/D  # variance of log odds ratio
  #var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  d <- lor * sqrt(3)/pi  # conversion to d
  var.d <- 3 * var.lor/pi^2  # variance of d
  out <- cbind(or, lor, var.lor, d, var.d)
  return(out)
}

##============ COMPUTATIONS TO CALCULATE CORRELATIONS ================##

# Formulas for computing r in designs with independent groups. 
# Section 12.4 & Table 12.4 (Cooper et al., 2009; pp. 231-234).

# (1) Computing variance of r,  z', and variance of z':

var_r <-  function(r, n) ((1 - r^2)^2)/(n - 1)  #calulate variance of r     
r_to_z  <-  function(r) 0.5*log((1 + r)/(1 - r))  #convert to z' 
var_z <- function(n) 1 / (n - 3)  #variance of z'

# Formulas for computing r in designs with independent groups. 
# Section 12.4 & Table 12.4 (Cooper et al., 2009; pp. 231-234).

# (1) Study reported: 
# t (t-test value of differences between 2 groups), n (total sample size)

r_from_t <- function(t, n) {
  r <- sqrt((t^2)/(t^2 + n-2))
  var_r <- ((1-r^2)^2)/(n-1)
  out <- cbind(r, var_r)
  return(out)
}

# Converting d (mean difference) to r where n.tmt = n.comparison 
# (Section 12.5.4; pp. 234)

r_from_d <- function(d,  var.d,  a=4) {
  r <- d/sqrt((d^2) + a)
  var_r <- (a^2*var.d)/(d^2 + a)^3
  out <- cbind(r, var_r)
  return(out)
}

# Converting d to r where n.tmt (not) = n.comparison (Section 12.5.4; pp. 234)

r_from_d1 <- function(d,  n.1, n.2,  var.d) {
  a <- ((n.1 + n.2)^2)/(n.1*n.2)
  r <- d/sqrt((d^2) + a)
  var_r <- (a^2*var.d)/(d^2 + a)^3
  out <- cbind(r, var_r)
  return(out)
  }

# Converting Chi-squared statistic with 1 df to r

r_from_chi <- function(chi.sq,  n) sqrt(chi.sq/n)

##============ WITHIN STUDY AGGREGATION OF EFFECT SIZES =============##

# Functions for aggregating within-study effect sizes (accounting for dependencies).
# Required inputs are a data.frame containing id (study id number) and 
# r (correlations).  This function fixes the correlation between predictor 
# variables at .50 and will compute the correct aggregated effect size for
# all studies. Function implements Hunter & Schmidt (2004)
# approach to aggregation of dependent r's (see chapter 10,  pp. 435-8).

aggrs <- function(r, cor = .50) {
  # Intermediate level function to compute weighted aggregation of effect sizes 
  k <- length(r)
   rbar <- cor
   r.xY <- sum(r)/(1 * sqrt(k  +  k*(k - 1) * rbar))
   return(r.xY)
}

agg_r <- function(meta, cor = .50) {
  id <- meta$id
  n <- meta$n
  r <- meta$r
  st <- unique(id)
  out <- data.frame(id=st)
    for(i in 1:length(st)) { 
    out$id[i] <- st[i]
    out$r[i] <- aggrs(r = r[id==st[i]], cor)
    out$n[i] <- round(mean(n[id==st[i]]),0)
  }
  return(out)
}  



 agg_r2<-function (meta, mod, cor = 0.5){
    id <- meta$id
    n <- meta$n
    r <- meta$r
    st <- as.factor(id)
    st <- unique(st)
    #st <- unique(id)
    um <- unique(mod)
    out <- data.frame(id = rep(st, rep(length(um), length(st))))
    out$mod <- rep(um, length(st))
    for (i in 1:length(st)) {
        for (j in 1:length(um)) {
            ro <- (i - 1) * length(um) + j
            m1 <- match(id, st[i], nomatch = 0)
            m2 <- match(mod, um[j], nomatch = 0)
            num <- sum(m1 * m2)
            out$r[ro] <- ifelse(num == 0, NA, aggrs(r = r[id == 
                st[i] & mod == um[j]], cor))
            out$n[ro] <- round(mean(n[id == st[i] & mod == um[j]]), 
                0)
        }
    }
    out2 <- out[is.na(out$r) == 0, ]
    return(out2)
}

##=== Add Fixed and Random Effects Weights ===##
# Required input is a data.frame with column names id (study id), 
# r (correlation coefficient),  and n (sample size).
 
MetaR <-  function(meta, cor = .50) {  
  # Computes within study aggregation and adds fixed and random effects weights 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for 
  #   each study.
  #   Aggregated effect size (one row per study) based on Hunter and Schmidt's (2004)
  #   approach to aggregation of dependent r's (see chapter 10,  pp. 435-8).
  #   Adds study weights based of a Fisher's z transformation of r's (see chapter 14, 
  #   Cooper et al., 2009)
  meta <- agg_r(meta, cor)
  meta$var.r <- var_r(meta$r, meta$n)
  meta$var.r <-  ifelse(is.na(meta$var.r), ((1-meta$r^2)^2)/(meta$n-1), meta$var.r)
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))  #computing r to z' for each study
  meta$var.z <- 1/(meta$n-3)  # computing var.z for each study
  meta$l.ci95z <- meta$z-1.96*sqrt(meta$var.z)     #create random ci for each study
  meta$u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
  meta$l.ci95 <- (exp(2*meta$l.ci95z)-1)/(exp(2*meta$l.ci95z) + 1)
  meta$u.ci95 <- (exp(2*meta$u.ci95z)-1)/(exp(2*meta$u.ci95z) + 1)  
  meta$z.score <- meta$z/sqrt(meta$var.z)
  meta$p.value <- 2*pnorm(abs(meta$z.score), lower.tail=FALSE)
  meta$wi <-  1/meta$var.z  # computing weight for each study
  meta$wiTi <- meta$wi*meta$z  # used to calculate omnibus
  meta$wiTi2 <- meta$wi*(meta$z)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(meta$r))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  meta$tau <- (Q-(k - 1))/comp  # Level 2 variance
  meta$tau <- ifelse(meta$tau <= 0, 0, meta$tau)
  meta$var.tau <- meta$tau + meta$var.z  # Random effects variance (within study var + between var) 
  meta$wi.tau <- 1/meta$var.tau  # Random effects weights
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  meta$l.ci95z <- NULL 
  meta$u.ci95z <- NULL
 return(meta)
}

##================= FIXED AND RANDOM EFFECTS OMNIBUS ===============##
# Function to calculate fixed and random effects omnibus effect size for correlations,  
# outputing omnibus effect size,  variance,  standard error,  upper and lower 
# confidence intervals,  and heterogeneity test.
# Required inputs are a data frame with id (study id) and r (correlation 
# coefficient)variable names.

OmnibusES<-  function(meta,  var="weighted" ) {
  # Computes fixed and random effects omnibus effect size for correlations.  
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   var:  "weighted" or "unweighted". "weighted" is the default. Use the 
  #   unweighted variance method only if Q is rejected and is very large relative to k.   
  # Returns:
  #   Fixed and random effects omnibus effect size, variance, standard error, 
  #   upper and lower confidence intervals, p-value, Q (heterogeneity test), I2
  #   (I-squared--proportion of total variation in tmt effects due to heterogeneity 
  #   rather than chance). 
  meta <- MetaR(meta)
  k <- length(!is.na(meta$r)) # number of studies
  df <- k-1 
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate omnibus
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate omnibus
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate omnibus
  Tz.agg <- sum.wiTi/sum.wi  # omnibus z' 
  var.Tz.agg <- 1/sum.wi  # omnibus var.z
  se.Tz.agg <- sqrt(var.Tz.agg) 
  z.value  <-  Tz.agg/se.Tz.agg
  p.value <-  2*pnorm(abs(z.value), lower.tail=FALSE)
  T.agg <- (exp(2*Tz.agg)-1)/(exp(2*Tz.agg) + 1)  # fixed effect (FE) omnibus effect size (r)
  var.T.agg <- (exp(2*var.Tz.agg)-1)/(exp(2*var.Tz.agg) + 1)  # FE omnibus var.r
  lower.ci.Tz <- Tz.agg-1.96*se.Tz.agg
  upper.ci.Tz <- Tz.agg + 1.96*se.Tz.agg
  lower.ci <- (exp(2*lower.ci.Tz)-1)/(exp(2*lower.ci.Tz) + 1)  # FE lower CI
  upper.ci <- (exp(2*upper.ci.Tz)-1)/(exp(2*upper.ci.Tz) + 1)  # FE upper CI
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  # FE homogeneity test
  I2 <- (Q-(k-1))/Q  # I-squared 
  I2 <- ifelse(I2<0, 0, I2)        
  I2 <- paste(round(I2*100, 4),  "%",  sep="")                        
  p.homog <- pchisq(Q, df, lower=FALSE)  # <.05 = sig. heterogeneity  
  # random effects #
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
  sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE)
  sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
  Tz.agg.tau <- sum.wiTi.tau/sum.wi.tau
  if(var == "weighted") {
    var.Tz.agg.tau <-  1/sum.wi.tau 
    se.Tz.agg.tau <- sqrt(var.Tz.agg.tau)
    z.valueR  <-  Tz.agg.tau/se.Tz.agg.tau
    p.valueR <-  2*pnorm(abs(z.valueR), lower.tail=FALSE)
    T.agg.tau <- (exp(2*Tz.agg.tau)-1)/(exp(2*Tz.agg.tau) + 1)
    var.T.agg.tau <- (exp(2*var.Tz.agg.tau)-1)/(exp(2*var.Tz.agg.tau) + 1)
    lower.ci.Tz.tau <- Tz.agg.tau-1.96*se.Tz.agg.tau
    upper.ci.Tz.tau <- Tz.agg.tau + 1.96*se.Tz.agg.tau
    lower.ci.tau <- (exp(2*lower.ci.Tz.tau)-1)/(exp(2*lower.ci.Tz.tau) + 1)
    upper.ci.tau <- (exp(2*upper.ci.Tz.tau)-1)/(exp(2*upper.ci.Tz.tau) + 1)
  }
  if(var == "unweighted") {  # unweighted variance method
    var.agg <- (sum(meta$z^2)-sum(meta$z)^2/k)/(k-1) #14.20
    q.num <- (1/k)*sum(meta$var.z)                                   
    unwgtvar.Tz.agg.tau <- var.agg-q.num  #14.22
    var.Tz.agg.tau <-  ifelse(unwgtvar.Tz.agg.tau <= 0, 0, unwgtvar.Tz.agg.tau)  #if var < 0,  its set to 0
    se.Tz.agg.tau <- sqrt(var.Tz.agg.tau)
    z.valueR  <-  Tz.agg.tau/se.Tz.agg.tau
    p.valueR <-  2*pnorm(abs(z.valueR), lower.tail=FALSE)
    T.agg.tau <- (exp(2*Tz.agg.tau)-1)/(exp(2*Tz.agg.tau) + 1)
    var.T.agg.tau <- (exp(2*var.Tz.agg.tau)-1)/(exp(2*var.Tz.agg.tau) + 1)
    lower.ci.Tz.tau <- Tz.agg.tau-1.96*se.Tz.agg.tau
    upper.ci.Tz.tau <- Tz.agg.tau + 1.96*se.Tz.agg.tau
    lower.ci.tau <- (exp(2*lower.ci.Tz.tau)-1)/(exp(2*lower.ci.Tz.tau) + 1)  
    upper.ci.tau <- (exp(2*upper.ci.Tz.tau)-1)/(exp(2*upper.ci.Tz.tau) + 1)
  } 

  Fixed <- list(FixedEffects=c(k=k, r=T.agg,  var.r=var.T.agg,  se=se.Tz.agg, 
                l.ci=lower.ci,  u.ci=upper.ci,  z.value=z.value,  p.value=p.value,
                Q=Q, df.Q=df, p_homog=p.homog, I2=I2))
  Random <- list(RandomEffects=c(k=k, r=T.agg.tau, var.r=var.T.agg.tau,  
                 se=se.Tz.agg.tau,  l.ci=lower.ci.tau, u.ci=upper.ci.tau, 
                 z.value=z.valueR, p.value=p.valueR, Q=Q, df.Q=df,  p_homog=p.homog, 
                 I2=I2))
  omni.data <- as.data.frame(c(Fixed, Random))      
  omni.data$Omnibus <- c("k", "ES", "var.ES", "SE", "CI.lower", 
                         "CI.upper", "Z", "p", "Q", "df", "p.h", 
                         "I2")
  omni.data <- omni.data[c(3, 1, 2)]
  omni.data <- as.data.frame(omni.data)
  row.names(omni.data) <- NULL
  return(omni.data)
}

# Now,  if there is significant heterogeneity (p_homog < .05),  look for moderators.

##================= Categorical Moderator Analysis ================##
# Single predictor categorical moderator function (fixed effects)

CatModf <-  function(meta,  mod) {
  # Computes single predictor categorical moderator analysis. Computations derived from 
  # chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  # Returns:
  #   Fixed effects moderator means per group, k per group, 95% confidence intervals,
  #   z-value, p-value, variances, standard errors, Q, df(Q), and I-squared.
  meta$mod <- as.character(mod)
  meta$z  <- 0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$k <- 1
  meta$TW <- 1/meta$var.z
  meta$TWD <- meta$TW*meta$z
  meta$TWDS <- meta$TWD*meta$z
  out <- aggregate(meta[c("k","TW","TWD", "TWDS")], by=meta["mod"], FUN=sum)
  out$mod <- as.character(out$mod)
  lastrow <- dim(out)[1] + 1
  out[lastrow, 2:5] <- apply(out[,-1], 2, FUN=sum)
  out$mod[lastrow] <- "Overall"
  out$z <- out$TWD/out$TW
  out$var.z <- 1/out$TW
  out$se.z <- sqrt(out$var.z)
  out$Q <- out$TWDS - out$TWD^2/out$TW
  out$df <- out$k - 1
  out$z.value <- out$z/out$se.z
  out$p.value <- 2*pnorm(abs(out$z.value), lower.tail=FALSE)
  out$L.95ci <- out$z-1.96*out$se.z
  out$U.95ci  <- out$z+1.96*out$se.z
  out$p_homog <- ifelse(out$df==0,  1,  pchisq(out$Q, out$df, lower=FALSE)) 
  out$I2 <- (out$Q-(out$df))/out$Q   #I-squared  
  out$I2 <- ifelse(out$I2<0, 0, out$I2)     
  out$I2 <- paste(round(out$I2*100, 4),  "%",  sep="")
  # convert back to r
  out$r <- (exp(2*out$z)-1)/(exp(2*out$z) + 1)
  out$var.r <- (exp(2*out$var.z)-1)/(exp(2*out$var.z) + 1)
  out$se.r <- sqrt(out$var.r)
  out$L.95ci <- (exp(2*out$L.95ci)-1)/(exp(2*out$L.95ci) + 1)
  out$U.95ci <- (exp(2*out$U.95ci)-1)/(exp(2*out$U.95ci) + 1)
  out$TW <- NULL
  out$TWD <- NULL
  out$TWDS <- NULL
  out$z <- NULL
  out$var.z <- NULL
  out$se.z <- NULL
  colnames(out) <- c("Mod", "k", "Q", "df", "Z", "p", "CI.lower",
                    "CI.upper", "p.h", "I2","ES", "Var", "SE")
  out <- out[, c(1, 2, 11, 12, 13, 7, 8,5, 6, 3, 4, 9, 10 )] 
  return(out)
}

# Fixed effect single predictor categorical moderator Q-statistic

CatModfQ <- function(meta,  mod) {
  # Computes fixed effect Q-statistic (homogeneity test) for single predictor categorical 
  # moderator analysis. Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  # Returns:
  #   Fixed effects moderator Q-statistic, Q-within & between, df(Qw & Qb), and 
  #   homogeneity p-value within & between levels.
  mod.sig <- CatModf(meta, mod)
  k <- mod.sig$k[mod.sig$Mod=="Overall"]                           #number of studies
  levels <- length(mod.sig$Mod)-1
  Qb.df <- levels-1 
  Qw.df <- k-levels
  Q <- mod.sig$Q[mod.sig$Mod=="Overall"]         #overall heterogeneity Q
  Qw <- sum(mod.sig$Q[!mod.sig$Mod=="Overall"])  #overall within-group heterogeneity statistic
  Qw_p.value <- 1-pchisq(Qw, Qw.df)      
  Qb <- Q-Qw                                     #overall between-group heterogeneity
  Qb_p.value <- 1-pchisq(Qb, Qb.df)
  mod.Qstat <- data.frame(Q, Qw, Qw.df, Qw_p.value, Qb, Qb.df, Qb_p.value)
  names(mod.Qstat) <- c("Q", "Qw", "df.w",  "p.w", "Qb", "df.b", "p.b")
  return(mod.Qstat)
}
         
# Function for planned comparisons between 2 levels of moderator (fixed effects)

CatCompf <- function(meta,  mod, x1, x2,  method= "post.hoc1") {
  # Directly compares 2 levels of a categorical moderator using a fixed effects model. 
  # Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  #   x1: One level of categorical moderator
  #   x2: Other level (comparison group) of same categorical moderator
  #   method: "post.hoc" assumes the comparision was not planned prior to conducting
  #           the meta-analysis. The other option "planned" assumes you have planned 
  #           a priori to compare these levels of the categorical moderator. 
  #           Default is "post.hoc1". 
  # Returns:
  #   Random effects moderator means per group, mean difference (d), variance of difference,
  #   p-value, and 95% confidence intervals. 
  modsig <- CatModf(meta, mod)
  modsig$Mod <- as.factor(modsig$Mod)
  com1 <- levels(modsig$Mod)[x1]  # first level of moderator
  com2 <- levels(modsig$Mod)[x2]  # second level of moderator
  x1.es <- modsig[modsig$Mod==com1, "ES"]
  x2.es <- modsig[modsig$Mod==com2, "ES"]
  x1.var <- modsig[modsig$Mod==com1, "Var"]
  x2.var <- modsig[modsig$Mod==com2, "Var"]
  g <- (-1)*x1.es + 1*x2.es  # pg 288 (Cooper et al.,  2009)
  var <- (-1)^2*x1.var + (1)^2*x2.var
  df <- 2-1
  chi.sqr <- g^2/var
  m <- meta 
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  meta$z <- 0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  if(method == "post.hoc1") {  # post-hoc comparison (Tukey HSD method)
    fit<-aov(meta$z~meta$mod,weights=meta$wi)
    fit<-TukeyHSD(fit)
  }
  if(method == "post.hoc2") {  # post-hoc comparison (Scheffe method)
    z <- g/sqrt(var)
    z2 <- z^2
    levels <- length(levels(as.factor(mod)))
    df.post <- levels-1 
    p.value <- 1-pchisq(z2, df.post)
    L.95ci <- g-1.96*sqrt(var)
    U.95ci <- g + 1.96*sqrt(var)
    fit <- data.frame(g, var, p.value, L.95ci, U.95ci)
    names(fit) <- c("diff", "var.diff",  
                    "p", "CI.lower", "CI.lower" )
  }
  if (method == "planned") {  # planned comparison (a priori)
   df <- 2-1
   chi.sqr <- g^2/var
   p.value <- 1-pchisq(chi.sqr, df)
   L.95ci <- g-1.96*sqrt(var)
   U.95ci <- g + 1.96*sqrt(var)
   fit <- data.frame(g, var, p.value, L.95ci, U.95ci)
   names(fit) <- c("diff", "var.diff",  
                   "p", "CI.lower", "CI.lower" )
  }
  return(fit) 
}

CatModr <-  function(meta,  mod) {
  # Computes single predictor categorical moderator analysis. Computations derived from 
  # chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with g (standardized mean diff),var.g (variance of g),
  #         n.1 (grp 1 sample size), n.2 (grp 2 sample size).
  #   mod: Categorical moderator variable used for moderator analysis.    
  # Returns:
  #   Fixed effects moderator means per group, k per group, 95% confidence intervals,
  #   z-value, p-value, variances, standard errors, Q, df(Q), and I-squared.
  meta$mod <- as.character(mod)
  meta$z  <- 0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$k <- 1
  meta$TW <- 1/meta$var.z
  meta$TWS <- meta$TW^2
  meta$TWD <- meta$TW*meta$z
  meta$TWDS <- meta$TWD*meta$z
  out <- aggregate(meta[c("k","TW","TWS", "TWD", "TWDS")], by=meta["mod"], FUN=sum)
  out$Q <- out$TWDS - out$TWD^2/out$TW
  Q <- sum(out$Q)
  out$df <- out$k - 1
  df <- sum(out$df)
  out$c <- out$TW - (out$TWS/out$TW) # 19.37 ch. 19 Borenstein (2009)
  c <- sum(out$c)
  TSw <- (Q - df)/c  # 19.38 ch. 19 Borenstein (2009)  
  tau <- ifelse(TSw < 0, 0, TSw)
  # add the tau variance to each indiv study & then compute TW again
  meta$var.tau <- meta$var.z + tau
  meta$TW.tau <- 1/meta$var.tau
  meta$TWD.tau <- meta$TW.tau*meta$z
  meta$TWDS.tau <- meta$TWD.tau*meta$z
  out <- aggregate(meta[c("k","TW.tau","TWD.tau", 
                   "TWDS.tau")], by=meta["mod"], FUN=sum)
  lastrow <- dim(out)[1] + 1
  out[lastrow, 2:5] <- apply(out[,-1], 2, FUN=sum)
  out$mod[lastrow] <- "Overall"
  out$z <- out$TWD.tau/out$TW.tau
  out$var.z <- 1/out$TW.tau
  out$z <- out$TWD.tau/out$TW.tau
  out$se.z <- sqrt(out$var.z)
  out$Q <- out$TWDS.tau - out$TWD.tau^2/out$TW.tau
  out$df <- out$k - 1
  out$z.value <- out$z/out$se.z
  out$p.value <- 2*pnorm(abs(out$z.value), lower.tail=FALSE)
  out$L.95ci <- out$z-1.96*out$se.z
  out$U.95ci  <- out$z+1.96*out$se.z
  out$p_homog <- ifelse(out$df==0,  1,  pchisq(out$Q, out$df, lower=FALSE)) 
  out$I2 <- (out$Q-(out$df))/out$Q   # these values are not to be used 
  out$I2 <- ifelse(out$I2<0, 0, out$I2)     
  out$I2 <- paste(round(out$I2*100, 4),  "%",  sep="")
  out$TW.tau <- NULL
  out$TWD.tau <- NULL
  out$TWDS.tau <- NULL
  # convert back to r
  out$r <- (exp(2*out$z)-1)/(exp(2*out$z) + 1)
  out$var.r <- (exp(2*out$var.z)-1)/(exp(2*out$var.z) + 1)
  out$se.r <- sqrt(out$var.r)
  out$L.95ci <- (exp(2*out$L.95ci)-1)/(exp(2*out$L.95ci) + 1)
  out$U.95ci <- (exp(2*out$U.95ci)-1)/(exp(2*out$U.95ci) + 1)
  out$z <- NULL
  out$var.z <- NULL
  out$se.z <- NULL
  colnames(out) <- c("Mod", "k", "Q", "df", "Z", "p", "CI.lower",
                    "CI.upper", "p.h", "I2","ES", "Var", "SE")
  out <- out[, c(1, 2, 11, 12, 13, 7, 8,5, 6, 3, 4, 9, 10 )] 
  return(out)
}


#Q-statistic function (random effects)

CatModrQ <-  function(meta,  mod) {
  # Computes random effect Q-statistic (homogeneity test) for single predictor categorical 
  # moderator analysis. Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  # Returns:
  #   Random effects moderator Q-statistic, Q-within & between, df(Qw & Qb), and 
  #   homogeneity p-value within & between levels.
  mod.sig <- CatModr(meta, mod)
  k <- mod.sig$k[mod.sig$Mod=="Overall"]                         #number of studies
  levels <- length(mod.sig$Mod)-1
  Qb.df <- levels-1 
  Qw.df <- k-levels
  Q <- mod.sig$Q[mod.sig$Mod=="Overall"]         #overall heterogeneity Q
  Qw <- sum(mod.sig$Q[!mod.sig$Mod=="Overall"])  #overall within-group heterogeneity statistic
  Qw_p.value <- 1-pchisq(Qw, Qw.df)      
  Qb <- Q-Qw                                     #overall between-group heterogeneity
  Qb_p.value <- 1-pchisq(Qb, Qb.df)
  mod.Qstat <- data.frame(Qb, Qb.df, Qb_p.value)
  names(mod.Qstat) <- c("Qb", "df.b", "p.b")
  return(mod.Qstat)
}

## new function 2.20.10
# Integrated function with fixed and random es, Q, etc
CatMod <- function(meta, mod) {
  fixed <- CatModf(meta,mod)
  random <-CatModr(meta,mod)
  random$Q <- fixed$Q
  random$df<- fixed$df
  random$p.h <- fixed$p.h
  random$I2 <- fixed$I2
  Qf <- CatModfQ(meta,mod) # fixed effect Q
  Qr <- CatModrQ(meta,mod)
  out<- list(Fixed=fixed, Q.fixed= Qf, Random=random, Q.random= Qr)
  return(out)
}


# Function for planned comparisons between 2 levels of moderator (random effects)

CatCompr <- function(meta,  mod, x1, x2,  method= "post.hoc1") {
  # Directly compares 2 levels of a categorical moderator using a random effects model. 
  # Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  #   x1: One level of categorical moderator
  #   x2: Other level (comparison group) of same categorical moderator
  #   method: "post.hoc" assumes the comparision was not planned prior to conducting
  #           the meta-analysis. The other option "planned" assumes you have planned 
  #           a priori to compare these levels of the categorical moderator. 
  #           Default is "post.hoc". 
  # Returns:
  #   Random effects moderator means per group, mean difference (d), variance of difference,
  #   p-value, and 95% confidence intervals. 
  modsig <- CatModr(meta, mod)
  modsig$Mod <- as.factor(modsig$Mod)
  com1 <- levels(modsig$Mod)[x1]  # first level of moderator
  com2 <- levels(modsig$Mod)[x2]  # second level of moderator
  x1.es <- modsig[modsig$Mod==com1, "ES"]
  x2.es <- modsig[modsig$Mod==com2, "ES"]
  x1.var <- modsig[modsig$Mod==com1, "Var"]
  x2.var <- modsig[modsig$Mod==com2, "Var"]
  g <- (-1)*x1.es + 1*x2.es  # pg 288 (Cooper et al.,  2009)
  var <- (-1)^2*x1.var + (1)^2*x2.var
  df <- 2-1
  chi.sqr <- g^2/var
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- meta$mod
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  if(method == "post.hoc1") {  # post-hoc comparison (Tukey HSD method)
    fit<-aov(meta$z~meta$mod,weights=meta$wi.tau)
    fit<-TukeyHSD(fit)
  }
  if(method == "post.hoc2") {  # post-hoc comparison (Scheffe method)
    z <- g/sqrt(var)
    z2 <- z^2
    levels <- length(levels(as.factor(mod)))
    df.post <- levels-1 
    p.value <- 1-pchisq(z2, df.post)
    L.95ci <- g-1.96*sqrt(var)
    U.95ci <- g + 1.96*sqrt(var)
    fit <- data.frame(g, var, p.value, L.95ci, U.95ci)
      names(fit) <- c("diff", "var.diff",  
      "p", "CI.lower", "CI.lower" )
  }
  if (method == "planned") {  # planned comparison (a priori)
   df <- 2-1
   chi.sqr <- g^2/var
   p.value <- 1-pchisq(chi.sqr, df)
   L.95ci <- g-1.96*sqrt(var)
   U.95ci <- g + 1.96*sqrt(var)
   fit <- data.frame(g, var, p.value, L.95ci, U.95ci)
   names(fit) <- c("diff", "var.diff",  
   "p", "CI.lower", "CI.lower" )
  }
  return(fit) 
}

#integrated catcomp function outputting both fixed and random

CatComp <- function(meta, mod, x1=NULL, x2=NULL,  method="post.hoc1") {
  fixed <- CatCompf(meta,mod, x1, x2, method)
  random <-CatCompr(meta,mod, x1, x2, method)
  out<- list(Fixed=fixed, Random=random)
  return(out)
}

# multifactor cat mod analysis [IN PROGRESS]:
#add relevant columns and convert back to r--it will works fine!
#MFCatMod <- function(meta, mod1, mod2) {
#  m <- Wifun(meta)
#  m$mod1 <- mod1
#  m$mod2 <- mod2
#  fixed <- ddply(m, c("mod1", "mod2"), summarise, sum.wi = sum(wi),
#           sum.wiTi = sum(wiTi), sum.wiTi2 = sum(wiTi2))
#  fixed$ES <- fixed$sum.wiTi/fixed$sum.wi
#  random <- ddply(m, c("mod1", "mod2"), summarise, sum.wi.tau = sum(wi.tau),
#           sum.wiTi.tau = sum(wiTi.tau), sum.wiTi2.tau = sum(wiTi2.tau))
#  random$ES <- random$sum.wiTi.tau/random$sum.wi.tau
#  out <- list(Fixed = fixed, Random = random)
#  return(out)
#}

##==== META-REGRESSION FUNCTIONS (for continuous & categorical moderators)====##

# Meta-regression functions that correct the standard errors in OLS regressions 
# (Cooper,  2009; pp. 289-290). 
# These functions flexibly allows for single or multivariate predictors in the meta-regression 
# with continuous,  categorical,  or both moderator types simultaneously.

MAreg1 <- function(meta, mod, method="random") {  # Single predictor meta-regression
  # Computes single predictor fixed or random effects meta-regression (continuous or categorical). 
  # Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Moderator variable used for meta-regression.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  # Returns:
  #   Fixed or random effects beta coefficients,  adjusted standard errors, adjusted t-value, 
  #   95% confidence intervals, and adjusted p-value.
  m <- meta 
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- ddply(m,  .(id),  summarize,  r = aggrs(r), n=mean(n), mod = pick_one(mod))
  #meta <- agg_r2(m,m$mod)
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #        function(.data) .data[sample(nrow(.data), 1),]))
  meta$z <- 0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- meta$mod
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-(k - 1))/comp  # random effects variance
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  if(method == "fixed") {
    reg0 <- lm(meta$z~1, weights=meta$wi)  # empty model
    reg <- lm(meta$z~mod, weights=meta$wi)  # model with mod  
    df <- anova(reg)["Residuals",  "Df"]
    ms.error <- anova(reg)["Residuals",  "Mean Sq"]
    t.crit <- qt(.975,  df)
    newSE.z <- summary(reg)$coef[, 2]/sqrt(ms.error)  # 15.20 
    Bs.z <- summary(reg)$coef[, 1]
    Bs <-  (exp(2*Bs.z)-1)/(exp(2*Bs.z) + 1)
    newSE <- (exp(2*newSE.z)-1)/(exp(2*newSE.z) + 1)
    t.adj <- Bs/newSE
    p.adj <- 2*pnorm(abs(t.adj), lower.tail=FALSE)
    lower.ci <- Bs-(t.crit*newSE)  # 95% CI
    upper.ci <- Bs + (t.crit*newSE)  # 95% CI
    modelfit <- MRfit(reg0,reg)   # internal function to assess fit
    #sig <- symnum(p.adj,  corr = FALSE,  na = FALSE,  
    #              cutpoints <- c(0,  0.001,  0.01,  0.05,  0.1,  1), 
    #              symbols <- c("***",  "**",  "*",  ".",  " ")) 
  }
  # random effects #
  if(method == "random") {
    reg0 <- lm(meta$z~1, weights=meta$wi.tau) 
    reg <- lm(meta$z~mod, weights=meta$wi.tau)
    df  <- anova(reg)["Residuals",  "Df"]
    ms.error <- anova(reg)["Residuals",  "Mean Sq"]
    t.crit <- qt(.975,  df)
    newSE.z  <- summary(reg)$coef[, 2]/sqrt(ms.error)  # 15.20 
    Bs.z  <- summary(reg)$coef[, 1]  # Need an option to keep in z' 
    Bs <- (exp(2*Bs.z)-1)/(exp(2*Bs.z) + 1)
    newSE <- (exp(2*newSE.z)-1)/(exp(2*newSE.z) + 1)
    t.adj <- Bs/newSE
    p.adj <- 2*pnorm(abs(t.adj), lower.tail=FALSE)
    lower.ci <- Bs-(t.crit*newSE)  # 95% CI
    upper.ci <- Bs + (t.crit*newSE)  # 95% CI
    modelfit <- MRfit(reg0,reg)
    #sig <- symnum(p.adj,  corr = FALSE,  na = FALSE,  
    #              cutpoints <- c(0,  0.001,  0.01,  0.05,  0.1,  1), 
    #              symbols <- c("***",  "**",  "*",  ".",  " ")) 
  }
  out <- data.frame(b=Bs, SE=newSE,  t=t.adj, CI.lower=lower.ci, CI.Upper=upper.ci, 
                      p=p.adj)
      
  return(list(out, modelfit))
}

##=== Multivariate Meta-Regression ===##

MAreg2  <-  function(reg) {  # Multivariate meta-regression
  # Computes multiple predictor fixed or random effects meta-regression (continuous and/or
  # categorical). Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   reg: Weighted linear regression saved as an object (e.g., 
  #        reg <- lm(data$z ~ data$mod8 + data$mod1, weights= data$wi.tau). The outcome 
  #        variable is Fisher's z and the predictor moderators can be either continuous or
  #        categorical. Weight the regression by either the fixed or random effect weight
  #        (e.g., fixed=data$wi and random=data$wi.tau)
  # Returns:
  #   Fixed or random effects multivariate beta coefficients, adjusted standard errors, 
  #   adjusted t-value, 95% confidence intervals, and adjusted p-value.
  df  <-  anova(reg)["Residuals",  "Df"]
  ms.error  <-  anova(reg)["Residuals",  "Mean Sq"]
  t.crit  <-  qt(.975,  df)
  newSE.z  <-  summary(reg)$coef[, 2]/sqrt(ms.error)  
  Bs.z  <-  summary(reg)$coef[, 1]
  Bs <-  (exp(2*Bs.z)-1)/(exp(2*Bs.z) + 1)
  newSE <-  (exp(2*newSE.z)-1)/(exp(2*newSE.z) + 1)
  t.adj  <-  Bs/newSE
  p.adj  <-  2*pnorm(abs(t.adj), lower.tail=FALSE)
  lower.ci <- Bs-(t.crit*newSE)     #95% CI
  upper.ci <- Bs + (t.crit*newSE)     #95% CI
  #sig <- symnum(p.adj,  corr = FALSE,  na = FALSE,  
  #           cutpoints <- c(0,  0.001,  0.01,  0.05,  0.1,  1), 
  #           symbols <- c("***",  "**",  "*",  ".",  " ")) 
  out  <-  data.frame(b=Bs, SE=newSE,  t=t.adj, CI.lower=lower.ci, CI.Upper=upper.ci, 
                      p=p.adj)  
  return(out)
}


##============= GRAPHICS =============##

# requires ggplot2

##=== Meta-regression scatterplot with weighted regression line ===##

MAregGraph <- function(meta, mod,  method="random",  modname=NULL,  title=NULL, ylim=c(0, 1)) {
  # Outputs a scatterplot from a fixed or random effects meta-regression (continuous and/or
  # categorical). Computations derived from chapter 14 and 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Moderator variable used for meta-regression.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   modname: Name of moderator to appear on x-axis of plot. Default is NULL.
  #   title: Plot title. Default is NULL.
  #   ylim: Limits of y-axis with the first argrument minimal value and second maximum value.
  #         Default is c(0,1).
  # Returns:
  #   Scatterplot with fixed or random effects regression line and where size of points are 
  #   based on study weights--more precise studies are larger. The ggplot2 package outputs the 
  #   rich graphics. 
  require('ggplot2')
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- agg_r2(m,m$mod)
  #meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #        function(.data) .data[sample(nrow(.data), 1),]))
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- meta$mod
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  if(method=="fixed") {
    congraph <- ggplot(meta,  aes(mod, z, weight=wi), na.rm=TRUE) + 
    geom_point(aes(size=wi), alpha=.6, na.rm=TRUE) + 
    geom_smooth(aes(group=1), method = lm,  se = FALSE) + 
    xlab(modname) + ylab("Effect Size") +  
    ylim(ylim) +  
    opts(title=title, legend.position = "none")
  }
  if(method=="random") {
    congraph <- ggplot(meta,  aes(mod, z), na.rm=TRUE) + 
    geom_point(aes(size=wi.tau), alpha=.6, na.rm=TRUE) + 
    geom_smooth(aes(group=1, weight=wi.tau), method = lm, se = FALSE,  na.rm=TRUE) + 
    xlab(modname) + 
    ylab("Effect Size")  + 
    #ylim(min(meta$z), 1) +  
    opts(title=title, legend.position = "none")
  }
  return(congraph)
}

##=== Categorical Moderator Graph ===##

# Intermediate level function to add mean to boxplot

stat_sum_single1  <-  function(fun,  geom="point",  weight=wi, ...) {
                               stat_summary(fun.y=fun,  shape=" + ",
                               geom=geom,  size = 7,  ...)      
}
stat_sum_single2  <-  function(fun,  geom="point",  weight=wi.tau, ...) {
                               stat_summary(fun.y=fun,  shape=" + ",  
                               geom=geom,  size = 7,  ...)      
}

CatModGraph <- function(meta, mod,  method="random",  modname=NULL,  title=NULL) {
  # Outputs a boxplot from a fixed or random effects moderator analysis.
  # Computations derived from chapter 14 and 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for analysis.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   modname: Name of moderator to appear on x-axis of plot. Default is NULL.
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Boxplot graph with median, interquartile range, max, min, and 
  #   outliers from a fixed or random effects categorical moderator analysis. Places
  #   jitter points for each study and the size of points are based on study weights--more 
  #   precise studies are larger. The ggplot2 package outputs the 
  #   rich graphics. 
  require('ggplot2')
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- agg_r2(m,m$mod)
  #meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #        function(.data) .data[sample(nrow(.data), 1),]))
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- length(meta$mod)
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  #meta <- meta
  #meta$mod <- mod
  #meta <- meta[!is.na(mod), ] 
  if(method=="fixed") {
    catmod <- ggplot(meta,  aes(factor(mod),  z, weight=wi), na.rm=TRUE) + 
                    geom_boxplot(outlier.size=2, na.rm=TRUE) + 
                    geom_jitter(aes(shape=factor(mod), size=wi.tau), alpha=.25) + 
                    #theme_bw() + 
                    xlab(modname) + 
                    ylab("Effect Size")  + 
                    opts(title=title, legend.position="none", na.rm=TRUE) + 
                    stat_sum_single2(mean)

  }  
  if(method=="random") {
    catmod <- ggplot(meta,  aes(factor(mod),  z, weight=wi.tau), na.rm=TRUE) + 
                    geom_boxplot(outlier.size=2, na.rm=TRUE) + 
                    geom_jitter(aes(shape=factor(mod), size=wi.tau), alpha=.25) + 
                    #theme_bw() + 
                    xlab(modname) + 
                    ylab("Effect Size")  + 
                    opts(title=title, legend.position="none", na.rm=TRUE) + 
                    stat_sum_single2(mean)
  }
  return(catmod)
}

##=== Forrest Plot ===##

ForestPlot <- function(meta, method="random", title=NULL) {
  # Outputs a forest plot from a fixed or random effects omnibus analysis.
  # Computations derived from chapter 14, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size) for each study.
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Forest plot with omnibus effect size (fixed or random), point for each study 
  #   where size of point is based on the study's precision (based primarily on 
  #   sample size) and 95% confidence intervals. The ggplot2 package outputs the rich graphics.  
  require('ggplot2') 
  meta <- meta 
  meta$id <- factor(meta$id)  # , levels=rev(id))
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  sum.wi <- sum(meta$wi, na.rm=TRUE)        
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE) 
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)   
   if(method=="fixed") {  
     Tz.agg <- sum.wiTi/sum.wi              
     var.Tz.agg <- 1/sum.wi                  
     se.Tz.agg <- sqrt(var.Tz.agg)          
     T.agg <- (exp(2*Tz.agg)-1)/(exp(2*Tz.agg) + 1)  
     var.T.agg <- (exp(2*var.Tz.agg)-1)/(exp(2*var.Tz.agg) + 1) 
     Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                        
     k <- sum(!is.na(meta$r))                                   
     df <- k-1  
     omnibus <- data.frame(id="Omnibus", r=T.agg)
     meta$l.ci95z <- meta$z-1.96*sqrt(meta$var.z)     #create fixed ci for each study
     meta$u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
     meta$l.ci95 <- (exp(2*meta$l.ci95z)-1)/(exp(2*meta$l.ci95z) + 1)
     meta$u.ci95 <- (exp(2*meta$u.ci95z)-1)/(exp(2*meta$u.ci95z) + 1)
     forest <- ggplot(meta,  aes(y = id, x=r))+    #aes(y = factor(id, levels=rev(levels(id))),  x = r))  +  
                    geom_vline(xintercept=0) + 
                    geom_point(data=omnibus, colour="red", size=8, shape=23) + 
                    #geom_point(aes(size=wi)) + 
                    opts(title=title,  legend.position="none") + 
                    geom_errorbarh(aes(xmin = l.ci95,  xmax=u.ci95), size=.3, alpha=.6) + 
                    geom_vline(colour="red", linetype=2,  xintercept=T.agg) + 
                    xlim(-1, 1) + 
                    xlab("Effect Size") + 
                    #scale_y_discrete(breaks = NA, labels=NA)+  # supress y-labels
                    ylab(NULL) 
  }
  if(method == "random") {
    comp <- sum.wi-sum.wi2/sum.wi
    Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
    k <- sum(!is.na(meta$r))                                  
    df <- k-1      
    tau <- (Q-k + 1)/comp 				#random effects variance
    meta$var.tau <- meta$var.z + tau
    meta$wi.tau <- 1/meta$var.tau
    meta$wiTi.tau <- meta$wi.tau*meta$z
    meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
    sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
    sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE) 
    sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
    Tz.agg.tau <- sum.wiTi.tau/sum.wi.tau
    var.Tz.agg.tau <-  1/sum.wi.tau                         #the following is inaccurate 14.23:  (Q - df)/comp
    se.Tz.agg.tau <- sqrt(var.Tz.agg.tau)
    T.agg.tau <- (exp(2*Tz.agg.tau)-1)/(exp(2*Tz.agg.tau) + 1)
    var.T.agg.tau <- (exp(2*var.Tz.agg.tau)-1)/(exp(2*var.Tz.agg.tau) + 1)
    omnibus.tau  <- data.frame(id="Omnibus", r=T.agg.tau)
    meta$l.ci95z <- meta$z-1.96*sqrt(meta$var.z)     #create random ci for each study
    meta$u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
    meta$l.ci95 <- (exp(2*meta$l.ci95z)-1)/(exp(2*meta$l.ci95z) + 1)
    meta$u.ci95 <- (exp(2*meta$u.ci95z)-1)/(exp(2*meta$u.ci95z) + 1)
    forest <- ggplot(meta,  aes(y = id,x = r))+    #aes(y = factor(id, levels=rev(levels(id))),  x = r))  +  
                  geom_vline(xintercept=0) + 
                  geom_point(data=omnibus.tau, colour="red",  size=8,  shape=23) + 
                  #geom_point(aes(size=wi.tau)) + 
                  opts(title=title,  legend.position="none") + 
                  geom_errorbarh(aes(xmin = l.ci95,  xmax=u.ci95), size=.3, alpha=.6) + 
                  geom_vline(colour="red", linetype=2,  xintercept=T.agg.tau) + 
                  xlim(-1, 1) + 
                  xlab("Effect Size") + 
                  #scale_y_discrete(breaks = NA, labels=NA)+  # supress y-labels
                  ylab(NULL) 
  }
  return(forest)
}

##=== Funnel Plot ===## 

FunnelPlot <- function(meta, method="random",  title=NULL) {
  # Outputs a funnel plot from a fixed or random effects omnibus analysis to assess for
  # publication bias in the meta-analysis.
  # Computations derived from chapter 14, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size) for each study.
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Funnel plot with omnibus effect size (fixed or random), point for each study 
  #   where size of point is based on the study's precision (based primarily on sample 
  #   size) and standard error lines to assess for publication bias. 
  #   The ggplot2 package outputs the rich graphics.     
  require('ggplot2')
  meta <- MetaR(meta)
  sum.wi <- sum(meta$wi, na.rm=TRUE)       
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)   
  if(method=="fixed") {  
    meta$se.z <- sqrt(meta$var.z)           
    Tz.agg <- sum.wiTi/sum.wi                
    var.Tz.agg <- 1/sum.wi                  
    se.Tz.agg <- sqrt(var.Tz.agg)          
    Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                   
    k <- sum(!is.na(meta$r))                                  
    df <- k-1  
    l.ci95z <- meta$z-1.96*sqrt(meta$var.z)    
    u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
    omnibus.z <- Tz.agg
    funnel <- ggplot(meta,  aes(y = se.z,  x = z))  +  
                    geom_vline(colour="black", linetype=1, 
                    xintercept=omnibus.z) + 
                    #geom_point(aes(size=wi)) + 
                    opts(title=title,  legend.position="none") + 
                    xlim(-1.7, 1.7) + 
                    ylim(.028, .5) + 
                    xlab("Fisher's z") + 
                    ylab("Standard Error") + 
                    stat_abline(intercept=omnibus.z/1.96, slope=(-1/1.96)) + 
                    stat_abline(intercept=(-omnibus.z/1.96), slope=1/1.96) + 
                    scale_y_continuous(trans="reverse")

  }
  if(method == "random") {
    meta$se.z.tau <- sqrt(meta$var.tau)
    sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)
    comp <- sum.wi-sum.wi2/sum.wi
    sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
    sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE)
    sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
    Tz.agg.tau <- sum.wiTi.tau/sum.wi.tau
    var.Tz.agg.tau <-  1/sum.wi.tau                        
    se.Tz.agg.tau <- sqrt(var.Tz.agg.tau)
    omnibus.z.tau  <- Tz.agg.tau
    l.ci95z <- meta$z-1.96*sqrt(meta$var.tau)    
    u.ci95z <- meta$z + 1.96*sqrt(meta$var.tau)
    l.ci95 <- (exp(2*l.ci95z)-1)/(exp(2*l.ci95z) + 1)
    u.ci95 <- (exp(2*u.ci95z)-1)/(exp(2*u.ci95z) + 1)
    meta$l.ci95 <- ifelse(l.ci95<=-1, -1, l.ci95)
    meta$u.ci95 <- ifelse(u.ci95>=1, 1, u.ci95)
    funnel <- ggplot(meta,  aes(y = se.z.tau,  x = z))  +  
                    geom_vline(colour="black", linetype=1, 
                    xintercept=omnibus.z.tau) + 
                    #geom_point(aes(size=wi.tau)) + 
                    opts(title=title,  legend.position="none") + 
                    xlim(-1.7, 1.7) + 
                    ylim(.028, .5) + 
                    xlab("Fisher's z") + 
                    ylab("Standard Error") + 
                    stat_abline(intercept=omnibus.z.tau/1.96, slope=(-1/1.96)) + 
                    stat_abline(intercept=(-omnibus.z.tau/1.96), slope=1/1.96) + 
                    scale_y_continuous(trans="reverse")
  }
  return(funnel)
}

##== Multivariate Moderator Graphs ==##

MultiModGraph <- function(meta, conmod,  catmod, method="random",  
                          conmod.name=NULL,  title=NULL) {
  # Outputs a scatterplot and boxplot faceted by the categorical moderator from a 
  # fixed or random effects moderator analysis. Computations derived from chapter 
  # 14 and 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   conmod: Continuous moderator variable used for analysis.
  #   catmod: Categorical moderator variable used for analysis.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   conmod.name: Name of continuous moderator variable to appear on x-axis of plot. 
  #   Default is NULL.
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Multivariate moderator scatterplot graph faceted by categorical moderator levels. Also
  #   places a weighted regression line based on either a fixed or random effects analysis. 
  #   The ggplot2 packages outputs the rich graphics.  
  require('ggplot2')
  m <- meta
  m$conmod <- conmod
  m$catmod <- catmod
  compl <- !is.na(m$conmod)& !is.na(m$catmod)
  meta <- m[compl, ]
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- length(meta$mod)
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  if(method=="fixed") {
    multimod <- ggplot(meta, aes(conmod, z, weight=wi), na.rm=TRUE) + 
                       opts(title=title, legend.position="none", na.rm=TRUE) + 
                       facet_wrap(~catmod)  + 
                       geom_point( aes(size=wi, shape=catmod)) + 
                       geom_smooth(aes(group=1, weight=wi),
                                   method= lm, se=FALSE, na.rm=TRUE) +
                       ylab("Effect Size") + 
                       xlab(conmod.name)
  }  
  if(method=="random") {
    multimod <- ggplot(meta, aes(conmod, z, weight=wi.tau), na.rm=TRUE) + 
                              opts(title=title, legend.position="none", na.rm=TRUE) + 
                              facet_wrap(~catmod)  + 
                              geom_point(aes(size=wi.tau, shape=catmod)) + 
                              geom_smooth(aes(group=1, weight=wi.tau), 
                                          method = lm, se = FALSE,  na.rm=TRUE) + 
                              ylab("Effect Size") + 
                              xlab(conmod.name)
  }
  return(multimod)
}


##===================== PUBLICATION BIAS ===================##
# Three approaches to assess for publication bias: (1) Fail
# Safe N, (2) Trim & Fill, and (3) Selection Modeling.

# Fail Safe N provides an estimate of the number of missing studies that 
# would need to exist to overturn the current conclusions.

PubBias <- function(data) {
  k <- length(data$z.score)
  sum.z <- sum(data$z.score)
  Z <- sum.z/sqrt(k)
  k0 <- round(-k + sum.z^2/(1.96)^2, 0)
  k.per <- round(k0/k, 0)
  out <- list("Fail.safe"=cbind(k,Z,k0, k.per))
  return(out)
}




##====== Additional Functions ========##

# Function to reduce data set with complete data for 1 predictor 

ComplData <- function(meta, mod, type= "independent", cor = .50) {   
  # Outputs an aggregated data.frame that will remove any missing data from the data 
  # set. This is particularly useful to output non-missing data based on a specified
  # number of variables (generally in conjunction with the multivariate moderator
  # functions above)
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod1: Moderator variable wanting to be kept for further analysis.    
  # Returns:
  #   Reduced data.frame (with complete data) for the moderator entered into the 
  #   function while aggregating based on recommended procedured (Hunter & Schmidt,
  #   2004).  
  m <- meta
  m$mod <- mod
  compl <- !is.na(m$mod)
  m <- m[compl, ]
  if(type == "independent") {
    meta <- agg_r2(m,  m$mod, cor)
    meta <- do.call(rbind, lapply(split(meta, meta$id), 
          function(.data) .data[sample(nrow(.data), 1),]))
  }
  if(type == "dependent") {
    meta <- agg_r2(m,  m$mod, cor)
  }
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- length(meta$mod)
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  return(meta)
}
         
    

# Correction for Attenuation

Rho_TU<- function(r,xx,yy) {
  r.corrected<-r/(sqrt(xx)*sqrt(yy))
  return(r.corrected)
   }


CorAtten <- function (meta, xx, yy) 
{
    m <- meta
    meta$xx <- xx
    meta$yy <- yy
    meta$r <- ifelse(is.na(meta$xx & meta$yy), meta$r, Rho_TU(meta$r, 
        meta$xx, meta$yy))
    meta$z <- 0.5 * log((1 + meta$r)/(1 - meta$r))
    meta$var.r <- ((1 - meta$r^2)^2)/(meta$n - 1)
    meta$var.r2 <- ifelse(is.na(m$xx & m$yy), meta$var.r, meta$var.r/(xx * 
        yy))
    meta$var.z <- 1/(meta$n - 3)
    meta$wi <- 1/meta$var.z
    meta$wiTi <- meta$wi * meta$z
    meta$wiTi2 <- meta$wi * (meta$z)^2
    meta$k <- length(meta$mod)
    mod <- meta$mod
    sum.wi <- sum(meta$wi, na.rm = TRUE)
    sum.wi2 <- sum(meta$wi2, na.rm = TRUE)
    sum.wiTi <- sum(meta$wiTi, na.rm = TRUE)
    sum.wiTi2 <- sum(meta$wiTi2, na.rm = TRUE)
    comp <- sum.wi - sum.wi2/sum.wi
    Q <- sum.wiTi2 - (sum.wiTi^2)/sum.wi
    k <- sum(!is.na(meta$r))
    df <- k - 1
    tau <- (Q - k + 1)/comp
    meta$var.tau <- meta$var.z + tau
    meta$wi.tau <- 1/meta$var.tau
    meta$wiTi.tau <- meta$wi.tau * meta$z
    meta$wiTi2.tau <- meta$wi.tau * (meta$z)^2
    return(meta)
}


Wifun <-  function(meta) {  
  meta$var.r <- var_r(meta$r, meta$n)
  meta$var.r <-  ifelse(is.na(meta$var.r), ((1-meta$r^2)^2)/(meta$n-1), meta$var.r)
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))  #computing r to z' for each study
  meta$var.z <- 1/(meta$n-3)  # computing var.z for each study
  meta$l.ci95z <- meta$z-1.96*sqrt(meta$var.z)     #create random ci for each study
  meta$u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
  meta$l.ci95 <- (exp(2*meta$l.ci95z)-1)/(exp(2*meta$l.ci95z) + 1)
  meta$u.ci95 <- (exp(2*meta$u.ci95z)-1)/(exp(2*meta$u.ci95z) + 1)  
  meta$z.score <- meta$z/sqrt(meta$var.z)
  meta$p.value <- 2*pnorm(abs(meta$z.score), lower.tail=FALSE)
  meta$wi <-  1/meta$var.z  # computing weight for each study
  meta$wiTi <- meta$wi*meta$z  # used to calculate omnibus
  meta$wiTi2 <- meta$wi*(meta$z)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(meta$r))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  meta$tau <- (Q-(k - 1))/comp  # Level 2 variance
  meta$var.tau <- meta$tau + meta$var.z  # Random effects variance (within study var + between var) 
  meta$wi.tau <- 1/meta$var.tau  # Random effects weights
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  meta$l.ci95z <- NULL 
  meta$u.ci95z <- NULL
 return(meta)
}

