##=================              MAd               ================##      
##================= Meta-Analysis Mean Differences ================##

# Package created by AC Del Re & William T. Hoyt
# This package contains all the relevant functions to conduct a
# mean differences (d and g) meta-analysis using standard procedures
# as described in Cooper, Hedges, & Valentine's Handbook of
# Research Synthesis and Meta-Analysis (2009).


## Robust SE based on Hedges et al., (2010) Eq. 6, Research Synthesis Methods 
# Coded by Mike Cheung, PhD (author of metaSEM package) with modifications
# by AC Del Re

# If the correlations of the dependent effect sizes are unknown, one
#  approach is to conduct the meta-analysis by assuming that the effect
#  sizes are independent. A robust standard error is then calculated to
#  adjust for the dependence. You may refer to Hedges et. al., (2010) for
#  more information. I have coded it here for reference.
# 
#  Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance
#  estimation in meta-regression with dependent effect size estimates.
#  Research Synthesis Methods, 1(1), 39-65. doi:10.1002/jrsm.5
# model: object fitted by metafor()
# cluster: indicator for clusters of studies


robustSE <- function(model, cluster=NULL, CI=.95, digits=3) {
  # Number of clusters assumed independent if not specified
  # model$not.na: complete cases
  if (is.null(cluster)) {
    m=length(model$X[model$not.na,1])
  } else {
    #m=nlevels(unique(as.factor(cluster[model$not.na])))
    m=length(unique(as.factor(cluster)))
  }
  res2 <- diag(residuals(model)^2)
  X <- model$X
  b <- model$b
  W <- diag(1/(model$vi+model$tau2))     
  meat <- t(X) %*% W %*% res2 %*% W %*% X    # W is symmetric
  bread <- solve( t(X) %*% W %*% X)
  V.R <- bread %*% meat %*% bread            # Robust sampling covariance matrix
  p <- length(b)                             # no. of predictors including intercept
  se <- sqrt( diag(V.R)*m/(m-p) )            # small sample adjustment (Eq.7)
  tval <- b/se
  pval <- 2*(1-pt(abs(tval),df=(m-p)))
  crit <- qt( (1-CI)/2, df=(m-p), lower.tail=FALSE )
  ci.lb <- b-crit*se
  ci.ub <- b+crit*se
  round(data.frame(estimate=b, se=se, t=tval,  ci.l=ci.lb, ci.u=ci.ub, p=pval), digits)
}


####==== META-REGRESSION (for continuous & categorical moderators)====####

mareg <- function(formula, var, data, method = "REML", subset, digits=3, ...) {
  UseMethod("mareg")
}

print.mareg <- function (x, digits = 3, ...){
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    coef <- as.vector(x$b)
    names(coef) <- rownames(x$b)
  print(round(coef, digits))
  #print.default(coef)
    invisible(x)
}
    
summary.mareg <- function(object, digits = 3,...) {
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
  table <- round(cbind(b, se, z, ci.lower, ci.upper, p), digits)
  colnames(table) <- c("estimate", "se", "z", "ci.l",
    "ci.u", "p")
  #rownames(table) <- object$predictors
  table2 <- round(cbind(QE, QE.df, QEp, QM, QM.df, QMp), digits)
  colnames(table2) <- c("QE", "QE.df", "QEp", "QM", "QM.df", "QMp")
  model <- list(coef = table, fit = table2, tau2_empty=tau2_empty,
   tau2=tau2, R2=R2)
  rownames(model$fit) <-NULL
  class(model) <- "summary.mareg"
  return(model)
  #return(round(model, digits))
}

print.summary.mareg <- function(x, digits=3,...) {
  cat("\n Model Results:", "", "\n", "\n")
  printCoefmat(x$coef, signif.stars = TRUE)
  cat("\n Heterogeneity & Fit:", "", "\n" ,"\n")
  print(round(x$fit,digits))
  invisible(x)
}


# formula based function for meta-regression 
mareg.default <- function(formula, var, data, method = "REML", 
  subset, ...) {
    require('metafor', quietly=TRUE)  # requires metafor's rma function
    suppressMessages(library(metafor))
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
    model0 <- metafor::rma(yi=yi, vi=vi, method=method, ...)
    model <- metafor::rma(yi,vi=vi,mods=mods, method=method, ...)
    predictors <- colnames(mods)
    model$out <- with(model, list(estimate=b, se=se, z=zval, ci.l=ci.lb,
      ci.u=ci.ub,  "Pr(>|z|)"=pval, predictors=predictors))
    model$out2 <- with(model, data.frame(QE, QEp, QM, QMp))
    model$call <- call
    model$tau2_empty <- model0$tau2
    model$R2 <-  1-(model$tau2 / model$tau2_empty)
    class(model) <- c("mareg", "rma.uni", "rma")
    return(model)
} 

print.summary.mareg <- function(x, ...) {
  cat("\n Model Results:", "", "\n", "\n")
  printCoefmat(x$coef, signif.stars = TRUE)
  cat("\n Heterogeneity & Fit:", "", "\n" ,"\n")
  printCoefmat(x$fit, signif.stars = TRUE)
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


# function to print objects to Word
wd <- function(object, get = FALSE, new = FALSE, ...) {
  UseMethod("wd")
}

wd.default <- function(object, get = FALSE, new = FALSE, ...) {
  require('R2wd', quietly = TRUE)
  get <- get
  open <- R2wd::wdGet(get)	# If no word file is open, it will start a new one
  if(new == TRUE){
    new <- R2wd::wdNewDoc("c:\\Temp.doc")
  }
  else{ 
    new <- NULL
  }  	
  wd1 <- round(data.frame(estimate=object$b, se=object$se, z=object$zval,
   ci.lower=object$ci.lb, ci.upper=object$ci.ub,  "p"=object$pval), 4)
  title <- R2wd::wdHeading(level = 2, " Model Results:")
  wd2 <-  round(data.frame(QE=object$QE, QEp=object$QEp, QM=object$QM, QMP=object$QMp), 4)
  obj1 <- R2wd::wdTable(wd1, ...)
  title2 <- R2wd::wdHeading(level = 2, "Heterogeneity & Fit:")
  obj2 <- R2wd::wdTable(wd2, ...)
  out <- (list(open, new, title, obj1, title2, obj2))
  class(out) <- "wd"
  return(out)
}
wd.mareg <- function(object, get = FALSE, new = FALSE, ...) {
  require('R2wd', quietly = TRUE)
  get <- get
  open <- R2wd::wdGet(get)	# If no word file is open, it will start a new one
  if(new == TRUE){
    new <- R2wd::wdNewDoc("c:\\Temp.doc")
  }
  else{ 
    new <- NULL
  }  	
  QE.df <- round(object$k - object$p, 0)
  QM.df <- round(object$m, 0)
  wd1 <- round(data.frame(estimate=object$b, se=object$se, z=object$zval,
   ci.l=object$ci.lb, ci.u=object$ci.ub,  "p"=object$pval), 4)
  title <- R2wd::wdHeading(level = 2, " Model Results:")
  wd2 <-  round(data.frame(QE=object$QE, QE.df = QE.df,
            QEp=object$QEp, QM=object$QM, QM.df=QM.df, QMp=object$QMp), 4)
  obj1 <- R2wd::wdTable(wd1)
  title2 <- R2wd::wdHeading(level = 2, "Heterogeneity & Fit:")
  obj2 <- R2wd::wdTable(wd2, ...)
  out <- (list(open, new, title, obj1, title2, obj2))
  class(out) <- "wd.mareg"
  return(out)
}

wd.omni <- function(object, get = FALSE, new = FALSE, ...) {
  require('R2wd', quietly = TRUE)
  get <- get
  open <- R2wd::wdGet(get)	# If no word file is open, it will start a new one
  if(new == TRUE){
    new <- R2wd::wdNewDoc("c:\\Temp.doc")
  }
  else{ 
    new <- NULL
  }  	
    k <- object$k
    estimate <- object$estimate
    #var <- object$var
    se <- object$se
    ci.l <- object$ci.l
    ci.u <- object$ci.u 
    z <- object$z
    p <- object$p
    Q <- object$Q
    Q.df <- object$df.Q
    Qp <- object$Qp
    I2 <- object$I2
    results1 <- data.frame(estimate, se, ci.l, ci.u, z, p)
    #results <- formatC(table, format="f", digits=digits)
    results1$k <- object$k
    results1 <- results1[c(8,1:7)]
    results2 <- data.frame(Q, Q.df, Qp)
    results2$I2 <- object$I2
    title <- R2wd::wdHeading(level = 2, " Omnibus Model Results:")
    obj1 <- R2wd::wdTable(results1)
    title2 <- R2wd::wdHeading(level = 2, "Heterogeneity:")
    obj2 <- R2wd::wdTable(results2)
    out <- (list(open, new, title, obj1, title2, obj2))
    class(out) <- "wd.omni"
    return(out)
}

wd.macat <- function(object, get = FALSE, new = FALSE, ...) {
  require('R2wd', quietly = TRUE)
  get <- get
  open <- R2wd::wdGet(get)	# If no word file is open, it will start a new one
  if(new == TRUE){
    new <- R2wd::wdNewDoc("c:\\Temp.doc")
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
    results1 <- data.frame(estimate, var, se, ci.l, ci.u, z, p,
      Q, df, p.h)
    #results <- formatC(table, format="f", digits=digits)
    results1$k <- x1$k
    results1$I2 <- x1$I2
    results1$mod <- x1$mod
    results1 <- results1[c(13,11, 1:10,12)]
    results2 <- round(data.frame(Q=Qoverall, Qw, Qw.df, Qw.p, Qb, Qb.df, Qb.p),4)
    results1[2:12] <-round(results1[2:12],3)
    title <- R2wd::wdHeading(level = 2, " Categorical Moderator Results:")
    obj1 <- R2wd::wdTable(results1)
    title2 <- R2wd::wdHeading(level = 2, "Heterogeneity:")
    obj2 <- R2wd::wdTable(results2)
    out <- (list(open, new, title, obj1, title2, obj2))
    class(out) <- "wd.macat"
    return(out)
}

# 3.18.10 internal for GUI: convert to factor

facts <- function(meta, mod) {
  meta[,mod] <- factor(meta[,mod])
  return(meta)
}


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

##============ COMPUTING EFFECT SIZE ESTIMATE ==============##

compute_ds <- function(n.1, m.1, sd.1, n.2, m.2, sd.2, data,
                       denom = "pooled.sd") {
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    args <- match(c("n.1", "m.1", "sd.1", "n.2", "m.2", "sd.2", "data"),
    names(mf), 0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    meta <- data
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.m.1 <- mf[[match("m.1", names(mf))]]
    m.1 <- eval(mf.m.1, data, enclos = sys.frame(sys.parent()))
    mf.sd.1 <- mf[[match("sd.1", names(mf))]]
    sd.1 <- eval(mf.sd.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.m.2 <- mf[[match("m.2", names(mf))]]
    m.2 <- eval(mf.m.2, data, enclos = sys.frame(sys.parent()))
    mf.sd.2 <- mf[[match("sd.2", names(mf))]]
    sd.2 <- eval(mf.sd.2, data, enclos = sys.frame(sys.parent()))
  if(denom == "pooled.sd"){
    meta$s.within <- sqrt(((n.1-1)*sd.1^2+
                      (n.2-1)*sd.2^2)/(n.1 + n.2-2))# pooled.sd
    meta$d <-(m.1-m.2)/meta$s.within
    meta$var.d <- ((n.1+n.2)/(n.1*n.2))+ 
                  ((meta$d^2)/(2*(n.1+n.2)))
    #meta$se.d <- sqrt(meta$var.d)
  }
  if(denom == "control.sd") {
    meta$d <-(m.1-m.2)/sd.2  # control.sd in denominator
    meta$var.d <- ((n.1+n.2)/(n.1*n.2))+ 
                  ((meta$d^2)/(2*(n.1+n.2)))
    #meta$se.d <- sqrt(meta$var.d)
  }
  return(meta)
}

# if already have d and var.d and want to derive g and var.g: 

compute_gs <- function(d, var.d, n.1, n.2, data) {
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("d", "var.d", "n.1", "n.2", "data"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  meta <- data
  mf.d <- mf[[match("d", names(mf))]]
  d <- eval(mf.d, data, enclos = sys.frame(sys.parent()))
  mf.var.d <- mf[[match("var.d", names(mf))]]
  var.d <- eval(mf.var.d, data, enclos = sys.frame(sys.parent()))
  mf.n.1 <- mf[[match("n.1", names(mf))]]
  n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
  mf.n.2 <- mf[[match("n.2", names(mf))]]
  n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
  df <- (n.1+n.2)-2
  j <- 1-(3/(4*df-1))
  meta$g <- j*d
  meta$var.g <- j^2*var.d
  meta$se.g <- sqrt(var.d)
  return(meta)
}

# will compute both d and g
compute_dgs <- function(n.1, m.1, sd.1, n.2, m.2, sd.2, data,
                       denom = "pooled.sd") {
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    args <- match(c("n.1", "m.1", "sd.1", "n.2", "m.2", "sd.2", "data"),
    names(mf), 0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    meta <- data
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.m.1 <- mf[[match("m.1", names(mf))]]
    m.1 <- eval(mf.m.1, data, enclos = sys.frame(sys.parent()))
    mf.sd.1 <- mf[[match("sd.1", names(mf))]]
    sd.1 <- eval(mf.sd.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.m.2 <- mf[[match("m.2", names(mf))]]
    m.2 <- eval(mf.m.2, data, enclos = sys.frame(sys.parent()))
    mf.sd.2 <- mf[[match("sd.2", names(mf))]]
    sd.2 <- eval(mf.sd.2, data, enclos = sys.frame(sys.parent()))
  if(denom == "pooled.sd"){
    meta$s.within <- sqrt(((n.1-1)*sd.1^2+
                      (n.2-1)*sd.2^2)/(n.1+n.2-2))# pooled.sd
    meta$d <-(m.1-m.2)/meta$s.within
    meta$var.d <- ((n.1+n.2)/(n.1*n.2))+ 
                  ((meta$d^2)/(2*(n.1+n.2)))
    df <- n.1+n.2-2
    j <- 1-(3/(4*df-1))
    meta$g <- j*meta$d
    meta$var.g <- j^2*meta$var.d
    meta$se.g <- sqrt(meta$var.g)
  }
  if(denom == "control.sd") {
    meta$d <-(m.1-m.2)/sd.2  # control.sd in denominator
    meta$var.d <- (n.1+n.2)/(n.1*n.2)+ 
                  (meta$d^2)/(2*(n.1+n.2))
    df <- n.1+n.2-2
    j <- 1-(3/(4*df-1))
    meta$g <- j*meta$d
    meta$var.g <- j^2*meta$var.d
    meta$se.g <- sqrt(meta$var.g)
  }
  return(meta)
}


##============ WITHIN STUDY AGGREGATION OF EFFECT SIZES =============##

# Automated within-study effect size aggregation function accounting for 
# dependencies among effect sizes g (unbiased estimate of d).
# Required inputs are n.1 (treatment sample size), n.2 (comparison sample size)
# id (study id), g (unbiased effect size).

aggs <- function(es,  n.1, n.2, cor = .50) {
  n.1 <- mean(n.1)   
  n.2 <- mean(n.2)    
  N_ES <- length(es)
  corr.mat <- matrix (rep(cor, N_ES^2), nrow=N_ES)
  diag(corr.mat) <- 1
  es1es2 <- cbind(es) %*% es
  PSI <- (8*corr.mat + es1es2*corr.mat^2)/(2*(n.1+n.2))
  PSI.inv <- solve(PSI)
  a <- rowSums(PSI.inv)/sum(PSI.inv)
  var <- 1/sum(PSI.inv)
  es <- sum(es*a)
  out<-cbind(es,var, n.1, n.2)
  return(out)
  }


# automated agg 
# newAggFn.R includes separate functions
# aggs1() and agg1() are for g1 (pooled sd in denom);
# aggs2() and agg2() are for g2 (sdC in denom).

# aggs1() is for g1 (pooled sd in denominator)
aggs1 <- function (es, n.1, n.2, cor = 0.5)   {
  n.1 <- mean(n.1)
  n.2 <- mean(n.2)
  N_ES <- length(es)
  corr.mat <- if (is.matrix(cor)==TRUE) cor  else matrix(rep(cor, N_ES^2), nrow = N_ES)
  diag(corr.mat) <- 1
  es1es2 <- cbind(es) %*% es
  PSI <- (1/n.1 + 1/n.2)*corr.mat + 0.5*es1es2*corr.mat^2/(n.1 + n.2)
  PSI.inv <- solve(PSI)
  a <- rowSums(PSI.inv)/sum(PSI.inv)
  var1 <- 1/sum(PSI.inv)
  es1 <- sum(es * a)
  out <- cbind(es1, var1, n.1, n.2)
  return(out)
}

# aggs2() is for g2 (sdC in denom)
aggs2 <- function(es, n.1, n.2, cor = 0.5)  {
  n.1 <- mean(n.1)
  n.2 <- mean(n.2)
  N_ES <- length(es)
  # corr.mat <- matrix(rep(cor, N_ES^2), nrow = N_ES)
  corr.mat <- if (is.matrix(cor)==TRUE) cor  else matrix(rep(N_ES, N_ES^2), nrow = N_ES)
  diag(corr.mat) <- 1
  es1es2 <- cbind(es) %*% es
  PSI <- (1/n.1 + 1/n.2)*corr.mat + 0.5*es1es2*corr.mat^2/(n.2)
  PSI.inv <- solve(PSI)
  a <- rowSums(PSI.inv)/sum(PSI.inv)
  var2 <- 1/sum(PSI.inv)
  es2 <- sum(es * a)
  out <- cbind(es2, var2, n.1, n.2)
  return(out)
}

# agg1() uses aggs1() to aggregate
agg1 <- function (id,  es, var, n.1, n.2, cor = 0.5, mod = NULL, data) 
{
  mf <- match.call()
  args <- match(c("id", "es", "var", "n.1", "n.2", "mod", "cor", 
                  "data"), names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.id <- mf[[match("id", names(mf))]]
  id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  mf.es <- mf[[match("es", names(mf))]]
  es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var <- mf[[match("var", names(mf))]]
  var <- eval(mf.var, data, enclos = sys.frame(sys.parent()))
  mf.n.1 <- mf[[match("n.1", names(mf))]]
  n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
  mf.n.2 <- mf[[match("n.2", names(mf))]]
  n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  data$var1 <- var
  if (is.null(mod)) {
    st <- unique(id)
    out <- data.frame(id = st)
    for (i in 1:length(st)) {
      out$id[i] <- st[i]
      out$es1[i] <- aggs1(es = es[id == st[i]], n.1 = n.1[id == st[i]], 
                         n.2 = n.2[id == st[i]], cor)[1]
      
      n_id <- sum(match(id, st[i], nomatch = 0))
      
      out$var1[i] <- ifelse(n_id == 1, data$var1[id == st[i]], 
                           aggs1(es = es[id == st[i]], n.1 = n.1[id == st[i]], 
                                 n.2 = n.2[id == st[i]], cor)[2])
      out$n.1[i] <- round(mean(n.1[id == st[i]]), 0)
      out$n.2[i] <- round(mean(n.2[id == st[i]]), 0)
    }
  }
  if (!is.null(mod)) {
    st <- unique(id)
    um <- unique(mod)
    out <- data.frame(id = rep(st, rep(length(um), length(st))),
                      mod <- rep(um, length(st)))
    for (i in 1:length(st)) {
      for (j in 1:length(um)) {
        ro <- (i - 1) * length(um) + j
        m1 <- match(id, st[i], nomatch = 0)
        m2 <- match(mod, um[j], nomatch = 0)
        num <- sum(m1 * m2)
        out$g1[ro] <- ifelse(num == 0, NA, aggs1(es = es[id == st[i] & mod == um[j]], 
                                                 n.1 = n.1[id == st[i] & mod == um[j]], 
                                                 n.2 = n.2[id == st[i] & mod == um[j]], cor)[1])
        out$var1[ro] <- ifelse(num == 0, NA, data$var1[id == st[i] & mod == um[j]])
        out$var1[ro] <- ifelse(num > 1, aggs1(es = es[id == st[i] & mod == um[j]], 
                                             n.1 = n.1[id == st[i] & mod == um[j]], 
                                             n.2 = n.2[id == st[i] & mod == um[j]], cor)[2], out$var1[ro])
        out$n.1[ro] <- round(mean(n.1[id == st[i] & mod == um[j]]), 0)
        out$n.2[ro] <- round(mean(n.2[id == st[i] & mod == um[j]]), 0)
      }
    }
    out <- out[is.na(out$es1) == 0, ]
  }
  return(out)
}

# agg2() uses aggs2() to aggregate
agg2 <- function (id,  es, var, n.1, n.2, cor = 0.5, mod = NULL, data) 
{
  mf <- match.call()
  args <- match(c("id", "es", "var", "n.1", "n.2", "mod", "cor", 
                  "data"), names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.id <- mf[[match("id", names(mf))]]
  id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  mf.es <- mf[[match("es", names(mf))]]
  es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var <- mf[[match("var", names(mf))]]
  var <- eval(mf.var, data, enclos = sys.frame(sys.parent()))
  mf.n.1 <- mf[[match("n.1", names(mf))]]
  n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
  mf.n.2 <- mf[[match("n.2", names(mf))]]
  n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  data$var1 <- var
  if (is.null(mod)) {
    st <- unique(id)
    out <- data.frame(id = st)
    for (i in 1:length(st)) {
      out$id[i] <- st[i]
      out$es2[i] <- aggs2(es = es[id == st[i]], n.1 = n.1[id == st[i]], 
                         n.2 = n.2[id == st[i]], cor)[1]
      
      n_id <- sum(match(id, st[i], nomatch = 0))
      
      out$var2[i] <- ifelse(n_id == 1, data$var1[id == st[i]], 
                           aggs2(es = es[id == st[i]], n.1 = n.1[id == st[i]], 
                                 n.2 = n.2[id == st[i]], cor)[2])
      out$n.1[i] <- round(mean(n.1[id == st[i]]), 0)
      out$n.2[i] <- round(mean(n.2[id == st[i]]), 0)
    }
  }
  if (!is.null(mod)) {
    st <- unique(id)
    um <- unique(mod)
    out <- data.frame(id = rep(st, rep(length(um), length(st))),
                      mod <- rep(um, length(st)))
    for (i in 1:length(st)) {
      for (j in 1:length(um)) {
        ro <- (i - 1) * length(um) + j
        m1 <- match(id, st[i], nomatch = 0)
        m2 <- match(mod, um[j], nomatch = 0)
        num <- sum(m1 * m2)
        out$es2[ro] <- ifelse(num == 0, NA, aggs2(es = es[id == st[i] & mod == um[j]], 
                                                 n.1 = n.1[id == st[i] & mod == um[j]], 
                                                 n.2 = n.2[id == st[i] & mod == um[j]], cor)[1])
        out$var2[ro] <- ifelse(num == 0, NA, data$var1[id == st[i] & mod == um[j]])
        out$var2[ro] <- ifelse(num > 1, aggs2(es = es[id == st[i] & mod == um[j]], 
                                             n.1 = n.1[id == st[i] & mod == um[j]], 
                                             n.2 = n.2[id == st[i] & mod == um[j]], cor)[2], out$var1[ro])
        out$n.1[ro] <- round(mean(n.1[id == st[i] & mod == um[j]]), 0)
        out$n.2[ro] <- round(mean(n.2[id == st[i] & mod == um[j]]), 0)
      }
    }
    out <- out[is.na(out$es2) == 0, ]
  }
  return(out)
}

# Two more functions for BHHR approach
aggsB <- function(es, var, cor=.5) {  # BHHR method to compute gComp, vgComp
  SE <- sqrt(var)
  SEmat <- SE %*% t(SE)  # cross-products
  rmat <-  if (is.matrix(cor)==TRUE) cor  else matrix(rep(cor, length(SEmat)), ncol=length(SE))
  diag(rmat) <- rep(1, length(SE))  # Option: user could input rmat instead of r
  VCV <- SEmat*rmat
  esComp <- mean(es)
  varComp <- sum(VCV)*(1/length(VCV))
  out <- cbind(esComp,varComp)
  return(out)
}


aggB <- function (id, es, var, cor = 0.5, data)  {
  mf <- match.call()
  args <- match(c("id", "es", "var", "cor", "data"), names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.id <- mf[[match("id", names(mf))]]
  id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  mf.es <- mf[[match("es", names(mf))]]
  es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var <- mf[[match("var", names(mf))]]
  var <- eval(mf.var, data, enclos = sys.frame(sys.parent()))
  
  st <- unique(id)
  out <- data.frame(id = st)
  
  for (i in 1:length(st)) {
    out$es[i] <- aggsB(es = es[id == st[i]], var = var[id == st[i]], cor)[1]
    out$var[i] <- aggsB(es = es[id == st[i]], var = var[id == st[i]], cor)[2] 
  }
  
  return(out)
}


##-----

# UPDATED AGG FUNCTION (12-2013)
agg <- function(id, es, var, n.1=NULL, n.2=NULL, method="BHHR", cor = .50, mod=NULL, data) {
  mf <- match.call()
  args <- match(c("id", "es", "var", "n.1", "n.2", "mod", "cor", 
                  "data"), names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.id <- mf[[match("id", names(mf))]]
  id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  mf.es <- mf[[match("es", names(mf))]]
  es <- eval(mf.es, data, enclos = sys.frame(sys.parent()))
  mf.var <- mf[[match("var", names(mf))]]
  var <- eval(mf.var, data, enclos = sys.frame(sys.parent()))
  mf.n.1 <- mf[[match("n.1", names(mf))]]
  n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
  mf.n.1 <- mf[[match("n.1", names(mf))]]
  n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
  mf.n.2 <- mf[[match("n.2", names(mf))]]
  n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  if(min(cor)<0 | max(cor)>1){
    stop("Assumed correlation between measures should be positive or zero (and not > 1)")
  }
  if(method=="GO1"){
    if(length(n.1)==0 | length(n.2)==0) stop("missing n.1 and/or n.2 arguments") else 
      out <- agg1(id=id, es=es, var=var, n.1=n.1, n.2=n.2, cor, mod, data)
  } else
    if(method=="GO2"){
      if(length(n.1)==0 | length(n.2)==0) stop("missing n.1 and/or n.2 arguments") else
        out <- agg2(id=id, es=es, var=var, n.1=n.1, n.2=n.2, cor, mod, data)
    } else
      if(method=="BHHR"){
        out <- aggB(id=id, es=es, var=var, cor, data)
      } else  stop("method must be either 'GO1', 'GO2', or 'BHHR'")
  return(out)
}

##=== Add Fixed and Random Effects Weights ===##
# Required input is a data.frame with column names id (study id), 
# g (unbiased standardized mean diff ES),  and n.1 (group 1 sample size),
# n.2 (group 2 sample size).
 

# required inputs are g, var.g and will output weights, ci, etc
wgts <-  function(g, var.g, data) {  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("g", "var.g", "data"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  #mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.g <- mf[[match("g", names(mf))]]
  meta <- data
  g <- eval(mf.g, data, enclos = sys.frame(sys.parent()))
  mf.var.g <- mf[[match("var.g", names(mf))]]
  var.g <- eval(mf.var.g, data, enclos = sys.frame(sys.parent()))
  meta$l.ci95 <- g-1.96*sqrt(var.g)     #create random ci for each study
  meta$u.ci95 <- g + 1.96*sqrt(var.g)
  meta$z.score <- g/sqrt(var.g)
  meta$wi <-  1/var.g  # computing weight for each study
  meta$wiTi <- meta$wi*g  # used to calculate omnibus
  meta$wiTi2 <- meta$wi*(g)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(g))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  meta$tau <- (Q-(k - 1))/comp  # Level 2 variance
  meta$var.tau <- meta$tau + var.g  # Random effects variance (within study var + between var) 
  meta$wi.tau <- 1/meta$var.tau  # Random effects weights
  meta$wiTi.tau <- meta$wi.tau*g
  meta$wiTi2.tau <- meta$wi.tau*(g)^2
  return(meta)
}


##================= FIXED AND RANDOM EFFECTS OMNIBUS ===============##
# Function to calculate fixed and random effects omnibus effect size for g,  
# outputing omnibus effect size,  variance,  standard error,  upper and lower 
# confidence intervals,  and heterogeneity test.
# Required input is a data.frame with column names id (study id), 
# g (unbiased standardized mean diff ES),  and n.1 (group 1 sample size),
# n.2 (group 2 sample size).

omni <-  function(g, var, data, type="weighted", method = "random") {
  # Computes fixed and random effects omnibus effect size for correlations.  
  # Args:
  #   g:    vector of unbiased estimate of standardized mean diff (d)
  #   var:  vector of variances of g.
  #   data: data.frame where g and var for each study are located.
  #   type:  "weighted" or "unweighted". "weighted" is the default. Use the 
  #   unweighted variance method only if Q is rejected and is very large relative to k.   
  # Returns:
  #   Fixed or random (Restricted Maximal Likelihood) effects omnibus effect size, variance, standard error, 
  #   upper and lower confidence intervals, p-value, Q (heterogeneity test), I2
  #   (I-squared--proportion of total variation in tmt effects due to heterogeneity 
  #   rather than chance). 
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("g", "var", "data", "type", "method"),
    names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.g <- mf[[match("g", names(mf))]]
  g <- eval(mf.g, data, enclos = sys.frame(sys.parent()))
  mf.var.g <- mf[[match("var", names(mf))]]
  var.g <- eval(mf.var.g, data, enclos = sys.frame(sys.parent()))
  #mf <- eval.parent(mf)
  #terms <- attr(mf, "terms")
  #g <- model.extract(mf, "g")
  #var.g <- model.extract(mf, "var")
  l.ci95 <- g-1.96*sqrt(var.g)     #create random ci for each study
  u.ci95 <- g + 1.96*sqrt(var.g)
  z.score <- g/sqrt(var.g)
  wi <-  1/var.g  # computing weight for each study
  wiTi <- wi*g  # used to calculate omnibus
  wiTi2 <- wi*(g)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(g))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  tau <- (Q-(k - 1))/comp  # Level 2 variance
  var.tau <- tau + var.g  # Random effects variance (within study var + between var) 
  wi.tau <- 1/var.tau  # Random effects weights
  wiTi.tau <- wi.tau*g
  wiTi2.tau <- wi.tau*(g)^2
  k <- length(!is.na(g)) # number of studies
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
  p.homog <- pchisq(Q, df, lower.tail=FALSE)  # <.05 = sig. heterogeneity  
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
  Fixed <- list(k=k, estimate=T.agg,  se=se.T.agg, 
                ci.l=lower.ci,  ci.u=upper.ci,  z=z.value,  p=p.value,
                Q=Q, df.Q=df, Qp=p.homog, I2=I2, call=call)
  Random <- list(k=k, estimate=T.agg.tau,   
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

print.omni <- function (x, digits = 3, ...){
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
    results1 <- round(data.frame(estimate, se,ci.l, ci.u, z, p),digits)
    #results <- formatC(table, format="f", digits=digits)
    results1$k <- x$k
    results1 <- results1[c(7,1:2,5,3:4,6)]
    results2 <- round(data.frame(Q, df.Q, Qp),digits)
    results2$I2 <- x$I2
    cat("\n Model Results:", "", "\n", "\n")
    print(results1)
    cat("\n Heterogeneity:", "", "\n","\n")
    print(results2)
    invisible(x)
}



##================= Categorical Moderator Analysis ================##

macat <- function (g, var, mod, data, method= "random") {
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
    args <- match(c("g", "var","mod", "data","method"),
    names(mf), 0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf.g <- mf[[match("g", names(mf))]]
    meta <- data
    meta$g <- eval(mf.g, data, enclos = sys.frame(sys.parent()))
    mf.var.g <- mf[[match("var", names(mf))]]
    meta$var.g <- eval(mf.var.g, data, enclos = sys.frame(sys.parent()))
    mf.mod <- mf[[match("mod", names(mf))]]
    meta$mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
    meta$mod <- as.character(meta$mod)
      meta$k <- 1
      meta$TW <- 1/meta$var.g
      meta$TWD <- meta$TW * meta$g
      meta$TWDS <- meta$TWD * meta$g
      temp <- aggregate(meta[c("k", "TW", "TWD", "TWDS")], by = meta["mod"], 
          FUN = sum)
      temp$mod <- as.character(temp$mod)
      lastrows <- dim(temp)[1] + 1
      temp[lastrows, 2:5] <- apply(temp[, -1], 2, FUN = sum)
      temp$mod[lastrows] <- "Overall"
      temp$g <- temp$TWD/temp$TW
      temp$var.g <- 1/temp$TW
      temp$se.g <- sqrt(temp$var.g)
      temp$Q <- temp$TWDS - temp$TWD^2/temp$TW
      temp$df <- temp$k - 1
      temp$z.value <- temp$g/temp$se.g
      temp$p.value <- 2 * pnorm(abs(temp$z.value), lower.tail = FALSE)
      temp$p_homog <- ifelse(temp$df == 0, 1, pchisq(temp$Q, temp$df, 
          lower.tail = FALSE))
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
      out$g <- out$TWD/out$TW
      out$var.g <- 1/out$TW
      out$se.g <- sqrt(out$var.g)
      out$Q <- out$TWDS - out$TWD^2/out$TW
      out$df <- out$k - 1
      out$z.value <- out$g/out$se.g
      out$p.value <- 2 * pnorm(abs(out$z.value), lower.tail = FALSE)
      out$L.95ci <- out$g - 1.96 * out$se.g
      out$U.95ci <- out$g + 1.96 * out$se.g
      out$p_homog <- ifelse(out$df == 0, 1, pchisq(out$Q, out$df, 
          lower.tail = FALSE))
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
      meta$var.tau <- meta$var.g + tau
      meta$TW.tau <- 1/meta$var.tau
      meta$TWD.tau <- meta$TW.tau * meta$g
      meta$TWDS.tau <- meta$TWD.tau * meta$g
      out <- aggregate(meta[c("k", "TW.tau", "TWD.tau", "TWDS.tau")], 
          by = meta["mod"], FUN = sum)
      lastrow <- dim(out)[1] + 1
      out[lastrow, 2:5] <- apply(out[, -1], 2, FUN = sum)
      out$mod[lastrow] <- "Overall"
      out$g <- out$TWD.tau/out$TW.tau
      out$var.g <- 1/out$TW.tau
      out$g <- out$TWD.tau/out$TW.tau
      out$se.g <- sqrt(out$var.g)
      out$Q <-  out$TWDS.tau - out$TWD.tau^2/out$TW.tau
      out$df <- out$k - 1
      out$z.value <- out$g/out$se.g
      out$p.value <- 2 * pnorm(abs(out$z.value), lower.tail = FALSE)
      out$L.95ci <- out$g - 1.96 * out$se.g
      out$U.95ci <- out$g + 1.96 * out$se.g
      out$p_homog <-  ifelse(out$df == 0, 1, pchisq(out$Q, out$df, 
        lower.tail = FALSE))
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
    out2 <- mod.Qstat
    out <- list(Model=out1, Heterogeneity=out2)
    class(out) <- "macat"
    return(out)
}

print.macat <- function (x, digits = 3, ...){
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
      Q, df, p.h), digits)
    #results <- formatC(table, format="f", digits=digits)
    results1$k <- x1$k
    results1$I2 <- x1$I2
    results1$mod <- x1$mod
    results1 <- results1[c(13,11, 1:10,12)]
    results2 <- round(data.frame(Q=Qoverall, Qw, Qw.df, Qw.p, Qb, Qb.df, Qb.p),digits)
    cat("\n Model Results:", "", "\n", "\n")
    print(results1)
    cat("\n Heterogeneity:", "", "\n","\n")
    print(results2)
    invisible(x)
}


      
# Function for planned comparisons between 2 levels of moderator (fixed effects)

macatC <- function(x1, x2, g, var, mod, data,  method= "random", type= "post.hoc") {
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
  args <- match(c("x1", "x.2", "g", "var", "mod", "data"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.g <- mf[[match("g", names(mf))]]
  g <- eval(mf.g, data, enclos = sys.frame(sys.parent()))
  mf.var <- mf[[match("var", names(mf))]]
  var.g <- eval(mf.var, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  modsig <- macat(g=g, var=var.g, mod=mod, data=data, method=method)
  modsig$mod <- as.factor(modsig$Model$mod)
  com1 <- levels(modsig$mod)[x1]  # first level of moderator
  com2 <- levels(modsig$mod)[x2]  # second level of moderator
  x1.es <- modsig$Model[modsig$mod==com1, "estimate"]
  x2.es <- modsig$Model[modsig$mod==com2, "estimate"]
  x1.var <- modsig$Model[modsig$mod==com1, "var"]
  x2.var <- modsig$Model[modsig$mod==com2, "var"]
  g <- (-1)*x1.es + 1*x2.es  # pg 288 (Cooper et al.,  2009)
  var <- (-1)^2*x1.var + (1)^2*x2.var
  df <- 2-1
  chi.sqr <- g^2/var
  #m <- meta 
  #m$mod <- mod
  #compl<-!is.na(m$mod)
  #meta<-m[compl,]
  #if(type == "post.hoc1") {  # post-hoc comparison (Tukey HSD method)
  #  fit<-aov(g~mod,weights=meta$wi)
  #  fit<-TukeyHSD(fit)
  #}
  if(type == "post.hoc") {  # post-hoc comparison (Scheffe method)
    z <- g/sqrt(var)
    z2 <- z^2
    levels <- length(levels(as.factor(mod)))
    df.post <- levels-1 
    p.value <- 1-pchisq(z2, df.post)
    L.95ci <- g-1.96*sqrt(var)
    U.95ci <- g + 1.96*sqrt(var)
    fit <- data.frame(g, var, p.value) # , L.95ci, U.95ci)
    names(fit) <- c("diff", "var.diff",  
                    "p") #, "CI.lower", "CI.lower" )
  }
  if (type == "planned") {  # planned comparison (a priori)
   df <- 2-1
   chi.sqr <- g^2/var
   p.value <- 1-pchisq(chi.sqr, df)
   L.95ci <- g-1.96*sqrt(var)
   U.95ci <- g + 1.96*sqrt(var)
   fit <- data.frame(g, var, p.value) # , L.95ci, U.95ci)
   names(fit) <- c("diff", "var.diff",  
                   "p") #, "ci.l", "ci.u" )
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

plotcon <- function(g, var, mod, data, method= "random", modname=NULL, 
  title=NULL, ...) {
  # Outputs a scatterplot from a fixed or random effects meta-regression (continuous and/or
  # categorical).  
  require('ggplot2')
  require('metafor')
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("g", "var","mod", "data","method"),
      names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.g <- mf[[match("g", names(mf))]]
  g <- eval(mf.g, data, enclos = sys.frame(sys.parent()))
  mf.var.g <- mf[[match("var", names(mf))]]
  var.g <- eval(mf.var.g, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  data$g <- g
  data$var.g <- var.g
  data$mod <- mod
  if(method=="fixed") {
    m1 <- mareg(g~as.numeric(mod), var=var.g, method="FE", data=data)
    #preds <- predict(dat)
    #preds <- as.data.frame(unclass(preds))
    congraph <- ggplot2::ggplot(data=data,ggplot2::aes(x=mod, y=g, size=1/var.g), na.rm=TRUE) +
                  ggplot2::geom_point(alpha=.7) +   
                  ggplot2::geom_abline(intercept=m1$b[1], slope=m1$b[2]) + 
                  ggplot2::expand_limits(x = 0) +
                  ggplot2::xlab(modname) + ggplot2::ylab("Effect Size") +  
                  ggplot2::ggtitle(title) +
                  ggplot2::theme_bw() +
                  ggplot2::theme(legend.position = "none")
  }
  if(method=="random") {
    m1 <- mareg(g~as.numeric(mod), var=var.g, method="REML", data=data)
    #preds <- predict(dat)
    #preds <- as.data.frame(unclass(preds))
    congraph <- ggplot2::ggplot(data=data,ggplot2::aes(x=mod, y=g, size=1/var.g), na.rm=TRUE) +
      ggplot2::geom_point(alpha=.7) +   
      ggplot2::geom_abline(intercept=m1$b[1], slope=m1$b[2]) + 
      ggplot2::expand_limits(x = 0) +
      ggplot2::xlab(modname) + ggplot2::ylab("Effect Size") +  
      ggplot2::ggtitle(title) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
  }
  return(congraph)
}

##=== Categorical Moderator Graph ===##

# Intermediate level function to add mean to boxplot

wmean <- function (x, weights = NULL, normwt = "ignored", na.rm = TRUE) 
{
  if (!length(weights)) 
    return(mean(x, na.rm = na.rm))
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  sum(weights * x)/sum(weights)
}

plotcat <- function(g, var, mod, data,  modname=NULL,  title=NULL, ...) {
  require('ggplot2')
  #require('quantreg')
  #require('SparseM')
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("g", "var","mod", "data"),
  names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf.g <- mf[[match("g", names(mf))]]
  m <- data
  m$g <- eval(mf.g, data, enclos = sys.frame(sys.parent()))
  mf.var.g <- mf[[match("var", names(mf))]]
  m$var.g <- eval(mf.var.g, data, enclos = sys.frame(sys.parent()))
  var.g <- NULL
  mf.mod <- mf[[match("mod", names(mf))]]
  m$mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  m$mod <- as.character(m$mod)
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  catmod <- ggplot2::ggplot(meta,  ggplot2::aes(factor(mod), g, weight = 1/var.g), na.rm=TRUE) + 
    ggplot2::geom_boxplot(ggplot2::aes(weight = 1/var.g), outlier.size=2,na.rm=TRUE) + 
    ggplot2::geom_jitter(ggplot2::aes(shape=factor(mod), size=1/var.g), alpha=.25) + 
    ggplot2::geom_point(ggplot2::aes(x=factor(mod),y=wmean(g, weights= 1/var.g, na.rm=TRUE)),shape = 23, 
               size = 3)+
    ggplot2::xlab(modname) + 
    ggplot2::ylab("Effect Size")  +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  return(catmod)
}

##=== Forest Plot ===##
## REVISE

#ForestPlot <- function(meta, method="random", title=NULL) {
  # Outputs a forest plot from a fixed or random effects omnibus analysis.
  # Computations derived from chapter 14, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with id, g (standardized mean diff),var.g (variance of g).
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Forest plot with omnibus effect size (fixed or random), point for each study 
  #   where size of point is based on the study's precision (based primarily on 
  #   sample size) and 95% confidence intervals. The ggplot2 package outputs the rich graphics.  
#  require('ggplot2')
#  meta <- meta 
#  meta$id <- factor(meta$id) # , levels=rev(id))                                  
#  meta$wi <-  1/meta$var.g
#  meta$wi2 <- (1/meta$var.g)^2
#  meta$wiTi <- meta$wi*meta$g
#  meta$wiTi2 <- meta$wi*(meta$g)^2
#  meta$k <- meta$mod
#  sum.wi <- sum(meta$wi, na.rm=TRUE)    
#  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
#  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
#  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
#  if(method=="fixed") {  
#     T.agg <- sum.wiTi/sum.wi              
#     var.T.agg <- 1/sum.wi                  
#     se.T.agg <- sqrt(var.T.agg)          
#     Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                        
#     k <- sum(!is.na(meta$g))                                   
#     df <- k-1  
#     omnibus <- data.frame(id="Omnibus", g=T.agg)
#     meta$l.ci95 <- meta$g-1.96*sqrt(meta$var.g)     #create fixed ci for each study
#     meta$u.ci95 <- meta$g + 1.96*sqrt(meta$var.g)
#     forest <- ggplot(meta,  ggplot2::aes(y = id,x=g))+       # aes(y = factor(id, levels=rev(levels(id))),  x = g))  +  
#                    geom_vline(xintercept=0) + 
#                    ggplot2::geom_point(data=omnibus, colour="red", size=8, shape=23) + 
#                    #ggplot2::geom_point(aes(size=wi)) + 
#                    opts(title=title,  legend.position="none") + 
#                    geom_errorbarh(data=meta,aes(xmin = l.ci95,  xmax=u.ci95), size=.3, alpha=.6) + 
#                    geom_vline(colour="red", linetype=2,  xintercept=T.agg) + 
#                    #xlim(-1, 1) + 
#                    xlab("Effect Size") + 
#                    #scale_y_discrete(breaks = NA, labels=NA)+  # supress y-labels
#                    ylab(NULL) 
#  }
#  if(method == "random") {
#    comp <- sum.wi-sum.wi2/sum.wi
#    Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
#    k <- sum(!is.na(meta$g))                                  
#    df <- k-1      
#    tau <- (Q-k + 1)/comp 				#random effects variance
#    meta$var.tau <- meta$var.g + tau
#    meta$wi.tau <- 1/meta$var.tau
#    meta$wiTi.tau <- meta$wi.tau*meta$g
#    meta$wiTi2.tau <- meta$wi.tau*(meta$g)^2
#    sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
#    sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE) 
#    sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
#    T.agg.tau <- sum.wiTi.tau/sum.wi.tau
#    var.T.agg.tau <-  1/sum.wi.tau                         #the following is inaccurate 14.23:  (Q - df)/comp
#    se.T.agg.tau <- sqrt(var.T.agg.tau)
#    omnibus.tau  <- data.frame(id="Omnibus", g=T.agg.tau)
#    meta$l.ci95 <- meta$g-1.96*sqrt(meta$var.g)     #create random ci for each study
#    meta$u.ci95 <- meta$g + 1.96*sqrt(meta$var.g)
#    forest <- ggplot(meta,  aes(y = id,x=g))+    #factor(id, levels=rev(levels(id))),  x = g))  +  
#                  geom_vline(xintercept=0) + 
#                  ggplot2::geom_point(data=omnibus.tau, colour="red",  size=8,  shape=23) + 
#                  #ggplot2::geom_point(aes(size=wi.tau)) + 
#                  opts(title=title,  legend.position="none") + 
#                  #ggplot2::geom_point(data=meta,aes(size=wi.tau)) + 
#                  geom_errorbarh(data=meta,aes(xmin = l.ci95,  xmax=u.ci95), size=.3, alpha=.6) + 
#                  geom_vline(colour="red", linetype=2,  xintercept=T.agg.tau) + 
#                  #xlim(-1, 1) + 
#                  xlab("Effect Size") + 
#                  #scale_y_discrete(breaks = NA, labels=NA)+  # supress y-labels
#                  ylab(NULL) 
#  }
#  return(forest)
#}

##=== Funnel Plot ===## 

#plotfunnel <- function(g, var, mod, data, method = "random", title=NULL) {
  # Outputs a funnel plot from a fixed or random effects omnibus analysis to assess for
  # publication bias in the meta-analysis.
  # Computations derived from chapter 14, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with ig, g (standardized mean diff),var.g (variance of g).
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Funnel plot with omnibus effect size (fixed or random), point for each study 
  #   where size of point is based on the study's precision (based primarily on sample 
  #   size) and standard error lines to assess for publication bias. 
  #   The ggplot2 package outputs the rich graphics.     
#  require('ggplot2')
#  call <- match.call()
#  mf <- match.call(expand.dots = FALSE)
#  args <- match(c("g", "var","mod", "data","method"),
#  names(mf), 0)
#  mf <- mf[c(1, args)]
#  mf$drop.unused.levels <- TRUE
#  mf[[1]] <- as.name("model.frame")
#  mf.g <- mf[[match("g", names(mf))]]
#  m <- data
#  m$g <- eval(mf.g, data, enclos = sys.frame(sys.parent()))
#  mf.var.g <- mf[[match("var", names(mf))]]
#  m$var.g <- eval(mf.var.g, data, enclos = sys.frame(sys.parent()))
#  mf.mod <- mf[[match("mod", names(mf))]]
#  m$mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
#  m$mod <- as.character(m$mod)
#  compl<-!is.na(m$mod)
#  meta<-m[compl,]
#  meta$wi <-  1/meta$var.g
#  meta$wi2 <- (1/meta$var.g)^2
#  meta$wiTi <- meta$wi*meta$g
#  meta$wiTi2 <- meta$wi*(meta$g)^2
#  sum.wi <- sum(meta$wi, na.rm=TRUE)       
#  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
#  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)   
#  if(method=="fixed") {  
#    meta$se <- sqrt(meta$var.g)           
##    T.agg <- sum.wiTi/sum.wi                
 #   var.T.agg <- 1/sum.wi                  
#    se.T.agg <- sqrt(var.T.agg)          
#    Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                   
#    k <- sum(!is.na(meta$g))                                  
#    df <- k-1  
#    omnibus <- T.agg
#    funnel <- ggplot(meta,  aes(y = se,  x = g))  +  
#                    geom_vline(colour="black", linetype=1, 
#                    xintercept=omnibus) + 
#                    #ggplot2::geom_point(aes(size=wi)) + 
#                    opts(title=title,  legend.position="none") + 
#                    xlim(-1.7, 1.7) + 
#                    ylim(.028, .5) + 
#                    xlab("g") + 
#                    ylab("Standard Error") + 
#                    stat_abline(intercept=omnibus/1.96, slope=(-1/1.96)) + 
#                    stat_abline(intercept=(-omnibus/1.96), slope=1/1.96) + 
#                    scale_y_continuous(trans="reverse")
#
#  }
#  if(method == "random") {
#    comp <- sum.wi-sum.wi2/sum.wi
#    Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
#    k <- sum(!is.na(meta$g))                                  
#    df <- k-1      
#    tau <- (Q-k + 1)/comp 				#random effects variance
#    meta$var.tau <- meta$var.g + tau
#    meta$se.tau <- sqrt(meta$var.tau)
#    sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)
#    comp <- sum.wi-sum.wi2/sum.wi
#    sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
#    sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE)
#    sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
#    T.agg.tau <- sum.wiTi.tau/sum.wi.tau
#    var.T.agg.tau <-  1/sum.wi.tau                        
#    se.T.agg.tau <- sqrt(var.T.agg.tau)
#    omnibus.tau  <- T.agg.tau
#    funnel <- ggplot(meta,  aes(y = se.tau,  x = g))  +  
#                    geom_vline(colour="black", linetype=1, 
#                    xintercept=omnibus.tau) + 
#                    #ggplot2::geom_point(aes(size=wi.tau)) + 
#                    opts(title=title,  legend.position="none") + 
#                    #xlim(-2, 2) + 
#                    #ylim(.028, .5) + 
#                    xlab("Fisher's z") + 
#                    ylab("Standard Error") + 
#                    stat_abline(intercept=omnibus.tau/1.96, slope=(-1/1.96)) + 
#                    stat_abline(intercept=(-omnibus.tau/1.96), slope=1/1.96) + 
#                    scale_y_continuous(trans="reverse")
#  }
#  return(funnel)
#}

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
  rater1 <- as.factor(rater1)
  rater2 <- as.factor(rater2)
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
  cat("\n")
  cat("strong agreement:    kappa > .75", "", "\n")
  cat("moderate agreement:  .40 < kappa < .75", "", "\n")
  cat("weak agreement:      kappa < .40", "", "\n","\n")
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



##====== Additional Functions ========##

# will not use this function any more--if want independent data, 
# people can use the function to do so (which I will display)

# Function to aggregate data set with complete data for 1 predictor 
# for data having been aggregated by 'agg' function with one mo
#ComplData <- function(meta, mod, type= "independent", cor = .50) {   
  # Outputs an aggregated data.frame that will remove any missing data from the data 
  # set. This is particularly useful to output non-missing data based on a specified
  # number of variables (generally in conjunction with the multivariate moderator
  # functions above)
  # Args:
  #   meta: data.frame with id, g (standardized mean diff),var.g (variance of g).
  #   mod1: Moderator variable wanting to be kept for further analysis.    
  # Returns:
  #   Reduced data.frame (with complete data) for the moderator entered into the 
  #   function while aggregating based on Cooper et al. recommended procedured (2009).  
#  m <- meta
#  m$mod <- mod
#  compl <- !is.na(m$mod)
#  m <- m[compl, ]
#  if(type == "independent") {
#    meta <- agg_g2(m, m$id, m$g, m$var.g, m$n.1, m$n.2, m$mod,cor)
#    meta <- do.call(rbind, lapply(split(meta, meta$id), 
#          function(.data) .data[sample(nrow(.data), 1),]))
#  }
#  if(type == "dependent") {
#     meta <- agg_g2(m, m$id, m$g, m$var.g, m$n.1, m$n.2,m$mod, cor)
#  }
#  meta$wi <-  1/meta$var.g
#  meta$wi2 <- (1/meta$var.g)^2
#  meta$wiTi <- meta$wi*meta$g
#  meta$wiTi2 <- meta$wi*(meta$g)^2
#  mod <- meta$mod
#  sum.wi <- sum(meta$wi, na.rm=TRUE)    
#  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
#  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
#  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
#  comp <- sum.wi-sum.wi2/sum.wi
#  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
#  k <- sum(!is.na(meta$g))                                  
#  df <- k-1      
#  tau <- (Q-(k - 1))/comp 				
#  meta$var.tau <- meta$var.g + tau
#  meta$wi.tau <- 1/meta$var.tau
#  meta$wiTi.tau <- meta$wi.tau*meta$g
#  meta$wiTi2.tau <- meta$wi.tau*(meta$g)^2
#  return(meta)
#}
         
##====== Multivariate Moderator Graphs ======##

#    multimod <- ggplot(meta, aes(conmod, g, weight=wi), na.rm=TRUE) + 
#                       opts(title=title, legend.position="none", na.rm=TRUE) + 
#                       facet_wrap(~catmod)  + 
#                       ggplot2::geom_point( aes(size=wi, shape=catmod)) + 
#                       geom_smooth(aes(group=1, weight=wi),
#                                   method= lm, se=FALSE, na.rm=TRUE) +
#                       ylab("Effect Size") + 
#                       xlab(conmod.name)
  


# Correction for Attenuation

atten <- function (g, xx, yy, data) {
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    args <- match(c("g", "xx", "yy", "data"),
    names(mf), 0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    meta <- data 
    mf.g <- mf[[match("g", names(mf))]]
    g <- eval(mf.g, data, enclos = sys.frame(sys.parent()))
    mf.xx <- mf[[match("xx", names(mf))]]
    xx <- eval(mf.xx, data, enclos = sys.frame(sys.parent()))
    mf.yy <- mf[[match("yy", names(mf))]]
    yy <- eval(mf.yy, data, enclos = sys.frame(sys.parent()))
    meta$g.corrected <-  g/(sqrt(xx)*sqrt(yy))
    return(meta)
}



#---- old functions (for now, will keep for GUI)-----------
#
##=================              MAd               ================##      
##================= Meta-Analysis Mean Differences ================##



# Package created by AC Del Re & William T. Hoyt
# This package contains all the relevant functions to conduct a
# mean differences (d and g) meta-analysis using standard procedures
# as described in Cooper, Hedges, & Valentine's Handbook of
# Research Synthesis and Meta-Analysis (2009).

# suggests('ggplot2')

# a formula for attenuation already created:  correct.cor(x, y)

##=== New Functions ===##

# 3.18.10 internal for GUI: convert to factor

facts <- function(meta, mod) {
  meta[,mod] <- factor(meta[,mod])
  return(meta)
}


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
# df$g <- as.numeric(as.character(df$g))

# set categorical variables to factors or character,  e.g.:
# df$id <- as.character(df$id)

# fix data with errors in factor names, requires car package,e.g.:
# library(car)
# df$outcome3 <- recode(df$outcome2, 'c("?",  "adherence", "compliance", "depression", 
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

prop_to_or <- function(p1, p2, n.ab, n.cd) {
  or <-(p1*(1-p2))/(p2*(1-p1))
  lor <- log(or)
  var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  out <- cbind(or,lor,var.lor)
  return(out)
}

prop_to_d <-function(p1, p2, n.ab, n.cd) {
  or <-(p1*(1-p2))/(p2*(1-p1))
  lor <- log(or)
  var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  d <- lor*sqrt(3)/pi
  var.d <- 3*var.lor/pi^2
  out <- cbind(or,lor,var.lor, d, var.d)
  return(out)
}

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

##============ COMPUTING EFFECT SIZE ESTIMATE ==============##

MeanDiffd <- function(meta, n.1= meta$n.1, m.1= meta$m.1, sd.1= meta$sd.1,
                       n.2= meta$n.2, m.2= meta$m.2, sd.2= meta$sd.2, 
                       denom = "pooled.sd") {
  if(denom == "pooled.sd"){
    meta$s.within <- sqrt(((n.1-1)*sd.1^2+
                      (n.2-1)*sd.2^2)/(n.1 + n.2-2))# pooled.sd
    meta$d <-(m.1-m.2)/meta$s.within
    meta$var.d <- ((n.1+n.2)/(n.1*n.2))+ 
                  ((meta$d^2)/(2*(n.1+n.2)))
    meta$se.d <- sqrt(meta$var.d)
  }
  if(denom == "control.sd") {
    meta$d <-(m.1-m.2)/sd.2  # control.sd in denominator
    meta$var.d <- ((n.1+n.2)/(n.1*n.2))+ 
                  ((meta$d^2)/(2*(n.1+n.2)))
    meta$se.d <- sqrt(meta$var.d)
  }
  return(meta)
}

# if already have d and var.d and want to derive g and var.g: 

MeanDiffg <- function(meta, d = meta$d, var.d = meta$var.d, n.1= meta$n.1, 
                   n.2 =  meta$n.2) {
  meta$df <- (n.1+n.2)-2
  meta$j <- 1-(3/(4*meta$df-1))
  meta$g <- meta$j*d
  meta$var.g <- meta$j^2*var.d
  meta$se.g <- sqrt(meta$var.d)
  meta$df <-NULL 
  meta$j <-NULL  
  return(meta)
}

# will compute both d and g
MeanDiff <- function(meta, n.1= meta$n.1, m.1= meta$m.1, sd.1= meta$sd.1,
                       n.2= meta$n.2, m.2= meta$m.2, sd.2= meta$sd.2, 
                       denom = "pooled.sd") {
  if(denom == "pooled.sd"){
    meta$s.within <- sqrt(((meta$n.1-1)*meta$sd.1^2+
                      (meta$n.2-1)*meta$sd.2^2)/(meta$n.1+meta$n.2-2))# pooled.sd
    meta$d <-(meta$m.1-meta$m.2)/meta$s.within
    meta$var.d <- ((meta$n.1+meta$n.2)/(meta$n.1*meta$n.2))+ 
                  ((meta$d^2)/(2*(meta$n.1+meta$n.2)))
    meta$df <- meta$n.1+meta$n.2-2
    meta$j <- 1-(3/(4*meta$df-1))
    meta$g <- meta$j*meta$d
    meta$var.g <- meta$j^2*meta$var.d
    meta$se.g <- sqrt(meta$var.g)
  }
  if(denom == "control.sd") {
    meta$d <-(meta$m.1-meta$m.2)/meta$sd.2  # control.sd in denominator
    meta$var.d <- (meta$n.1+meta$n.2)/(meta$n.1*meta$n.2)+ 
                  (meta$d^2)/(2*(meta$n.1+meta$n.2))
    meta$df <- meta$n.1+meta$n.2-2
    meta$j <- 1-(3/(4*meta$df-1))
    meta$g <- meta$j*meta$d
    meta$var.g <- meta$j^2*meta$var.d
    meta$se.g <- sqrt(meta$var.g)
  }
  meta$df <-NULL 
  meta$j <-NULL
  return(meta)
}


##============ WITHIN STUDY AGGREGATION OF EFFECT SIZES =============##

# Automated within-study effect size aggregation function accounting for 
# dependencies among effect sizes g (unbiased estimate of d).
# Required inputs are n.1 (treatment sample size), n.2 (comparison sample size)
# id (study id), g (unbiased effect size).

aggs <- function(g,  n.1, n.2, cor = .50) {
  n.1 <- mean(n.1)   
  n.2 <- mean(n.2)    
  #nT = nC = 0.5*nT
  N_ES <- length(g)
  corr.mat <- matrix (rep(cor, N_ES^2), nrow=N_ES)
  diag(corr.mat) <- 1
  g1g2 <- cbind(g) %*% g
  PSI <- (8*corr.mat + g1g2*corr.mat^2)/(2*(n.1+n.2))
  #PSI <- (1/nT + 1/nC)*corr.mat - (0.5*g1g2*corr.mat^2)/(nC+nT)
  PSI.inv <- solve(PSI)
  a <- rowSums(PSI.inv)/sum(PSI.inv)
  var.g <- 1/sum(PSI.inv)
  g <- sum(g*a)
  out<-cbind(g,var.g, n.1, n.2)
  return(out)
  }


# automated agg 
agg_g <- function(meta,id, g, var.g, n.1, n.2, cor = .50) {
  meta$var.g <- var.g
  st <- unique(id)       
  out <- data.frame(id=st)
    for(i in 1:length(st)) { 
    out$id[i] <- st[i]
    out$g[i] <- aggs(g=g[id==st[i]], n.1= n.1[id==st[i]],
                       n.2 = n.2[id==st[i]], cor)[1]
    lid <- sum(match(id, st[i], nomatch = 0))
    out$var.g[i]<-ifelse(lid ==1, meta$var.g[id==st[i]], 
                     aggs(g=g[id==st[i]], 
                          n.1= n.1[id==st[i]],
                          n.2 = n.2[id==st[i]], cor)[2])
    out$n.1[i] <- round(mean(n.1[id==st[i]]),0)
    out$n.2[i] <- round(mean(n.2[id==st[i]]),0)
  }
  return(out)
}

#additional agg functions for mods, etc
agg_g2 <- function(meta, id, g, var.g, n.1, n.2, mod, cor = .50) {
  meta$var.g <- var.g     # What's this for?  (Why would they be different?)
  st <- unique(id)   
  um <- unique(mod)
  # Initialize id, mod for all possible combinations.  Delete NAs at end.
  out <- data.frame(id=rep(st, rep(length(um), length(st))))
  out$mod <- rep(um, length(st))
  for(i in 1:length(st)) {   
    for(j in 1:length(um)) { 
      # row of df to fill
      ro <- (i-1)*length(um) + j
      # Are there any rows in meta where id==i and mod==j?
       m1<-match(id,st[i],nomatch=0)  # (rows where id==i)
       m2<-match(mod,um[j],nomatch=0)  # (rows where mod==j)
       #  m1*m2 will = 1 for each row in which both are true.
       # sum(m1*m2) gives the number of rows for which both are true.
       num <- sum(m1*m2)
       out$g[ro] <- ifelse(num==0,NA, #NA,
                      aggs(g=g[id==st[i]&mod==um[j]],
                      n.1=n.1[id==st[i]&mod==um[j]],
                      n.2=n.2[id==st[i]&mod==um[j]], cor)[1])
       # lid <- sum(match(id, st[i], nomatch = 0))  [num takes place of lid.]
       out$var.g[ro]<-ifelse(num==0, NA, 
                             meta$var.g[id==st[i]&mod==um[j]])
       out$var.g[ro]<-ifelse(num > 1, 
                             aggs(g=g[id==st[i]&mod==um[j]],
                             n.1= n.1[id==st[i]&mod==um[j]],
                             n.2 = n.2[id==st[i]&mod==um[j]],
                               cor)[2],
                               out$var.g[ro])
       out$n.1[ro] <- round(mean(n.1[id==st[i]&mod==um[j]]),0)
       out$n.2[ro] <- round(mean(n.2[id==st[i]&mod==um[j]]),0)
    }
  }
  # Strip out rows with no data.
  out2 <- out[is.na(out$g)==0,]
  return(out2)
}


##=== Add Fixed and Random Effects Weights ===##
# Required input is a data.frame with column names id (study id), 
# g (unbiased standardized mean diff ES),  and n.1 (group 1 sample size),
# n.2 (group 2 sample size).
 
MetaG <-  function(meta, cor = .50) {  
  meta <- agg_g(meta, meta$id, meta$g, meta$var.g, meta$n.1, meta$n.2, cor)
  meta$l.ci95 <- meta$g-1.96*sqrt(meta$var.g)     #create random ci for each study
  meta$u.ci95 <- meta$g + 1.96*sqrt(meta$var.g)
  meta$z.score <- meta$g/sqrt(meta$var.g)
  meta$p.value <- 2*(1-pt(abs(meta$z.score),  (meta$n.1 + meta$n.2) -1))
  meta$wi <-  1/meta$var.g  # computing weight for each study
  meta$wiTi <- meta$wi*meta$g  # used to calculate omnibus
  meta$wiTi2 <- meta$wi*(meta$g)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(meta$g))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  meta$tau <- (Q-(k - 1))/comp  # Level 2 variance
  meta$tau <- ifelse(meta$tau <= 0, 0, meta$tau)
  meta$var.tau <- meta$tau + meta$var.g  # Random effects variance (within study var + between var) 
  meta$wi.tau <- 1/meta$var.tau  # Random effects weights
  meta$wiTi.tau <- meta$wi.tau*meta$g
  meta$wiTi2.tau <- meta$wi.tau*(meta$g)^2
  return(meta)
}

# required inputs are id, g, var.g and will output weights, ci, etc
Wifun <-  function(meta) {  
  meta$l.ci95 <- meta$g-1.96*sqrt(meta$var.g)     #create random ci for each study
  meta$u.ci95 <- meta$g + 1.96*sqrt(meta$var.g)
  meta$z.score <- meta$g/sqrt(meta$var.g)
  #meta$p.value <- 2*(1-pt(abs(meta$z.score),  (meta$n.1+meta$n.2)-1))
  meta$wi <-  1/meta$var.g  # computing weight for each study
  meta$wiTi <- meta$wi*meta$g  # used to calculate omnibus
  meta$wiTi2 <- meta$wi*(meta$g)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(meta$g))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  meta$tau <- (Q-(k - 1))/comp  # Level 2 variance
  meta$var.tau <- meta$tau + meta$var.g  # Random effects variance (within study var + between var) 
  meta$wi.tau <- 1/meta$var.tau  # Random effects weights
  meta$wiTi.tau <- meta$wi.tau*meta$g
  meta$wiTi2.tau <- meta$wi.tau*(meta$g)^2
  return(meta)
}


##================= FIXED AND RANDOM EFFECTS OMNIBUS ===============##
# Function to calculate fixed and random effects omnibus effect size for g,  
# outputing omnibus effect size,  variance,  standard error,  upper and lower 
# confidence intervals,  and heterogeneity test.
# Required input is a data.frame with column names id (study id), 
# g (unbiased standardized mean diff ES),  and n.1 (group 1 sample size),
# n.2 (group 2 sample size).

OmnibusES<-  function(meta,  var="weighted", cor = .50 ) {
  # Computes fixed and random effects omnibus effect size for correlations.  
  # Args:
  #   meta: data.frame with g (standardized mean diff) and n (sample size)for each study.
  #   var:  "weighted" or "unweighted". "weighted" is the default. Use the 
  #   unweighted variance method only if Q is rejected and is very large relative to k.   
  # Returns:
  #   Fixed and random effects omnibus effect size, variance, standard error, 
  #   upper and lower confidence intervals, p-value, Q (heterogeneity test), I2
  #   (I-squared--proportion of total variation in tmt effects due to heterogeneity 
  #   rather than chance). 
  meta <- MetaG(meta, cor)
  k <- length(!is.na(meta$g)) # number of studies
  df <- k-1 
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate omnibus
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate omnibus
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate omnibus
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
  p.homog <- pchisq(Q, df, lower.tail=FALSE)  # <.05 = sig. heterogeneity  
  # random effects #
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
  sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE)
  sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
  T.agg.tau <- sum.wiTi.tau/sum.wi.tau
  if(var == "weighted") {
    var.T.agg.tau <-  1/sum.wi.tau 
    se.T.agg.tau <- sqrt(var.T.agg.tau)
    z.valueR  <-  T.agg.tau/se.T.agg.tau
    p.valueR <-  2*pnorm(abs(z.valueR), lower.tail=FALSE)
    lower.ci.tau <- T.agg.tau-1.96*se.T.agg.tau
    upper.ci.tau <- T.agg.tau + 1.96*se.T.agg.tau
  }
  if(var == "unweighted") {  # unweighted variance method
    var.agg <- (sum(meta$g^2)-sum(meta$g)^2/k)/(k-1) #14.20
    #q.num <- (1/k)*sum(meta$var.g)                                   
    unwgtvar.T.agg.tau <- var.agg-(1/k)*sum(meta$var.g)  #14.22
    var.T.agg.tau <-  ifelse(unwgtvar.T.agg.tau <= 0, 0, unwgtvar.T.agg.tau)  #if var < 0,  its set to 0
    se.T.agg.tau <- sqrt(var.T.agg.tau)
    z.valueR  <-  T.agg.tau/se.T.agg.tau
    p.valueR <-  2*pnorm(abs(z.valueR), lower.tail=FALSE)

    lower.ci.tau <- T.agg.tau-1.96*se.T.agg.tau
    upper.ci.tau <- T.agg.tau + 1.96*se.T.agg.tau
  } 
  Fixed <- list(FixedEffects=c(k=k, r=T.agg,  var.g=var.T.agg,  se=se.T.agg, 
                l.ci=lower.ci,  u.ci=upper.ci,  z.value=z.value,  p.value=p.value,
                Q=Q, df.Q=df, p_homog=p.homog, I2=I2))
  Random <- list(RandomEffects=c(k=k, r=T.agg.tau, var.g=var.T.agg.tau,  
                 se=se.T.agg.tau,  l.ci=lower.ci.tau, u.ci=upper.ci.tau, 
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

# 03-26-10 update: fixed & random effects

CatModf <-  function(meta,  mod) {
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
  meta$k <- 1
  meta$TW <- 1/meta$var.g
  meta$TWD <- meta$TW*meta$g
  meta$TWDS <- meta$TWD*meta$g
  out <- aggregate(meta[c("k","TW","TWD", "TWDS")], by=meta["mod"], FUN=sum)
  out$mod <- as.character(out$mod)
  lastrow <- dim(out)[1] + 1
  out[lastrow, 2:5] <- apply(out[,-1], 2, FUN=sum)
  out$mod[lastrow] <- "Overall"
  out$g <- out$TWD/out$TW
  out$var.g <- 1/out$TW
  out$se.g <- sqrt(out$var.g)
  out$Q <- out$TWDS - out$TWD^2/out$TW
  out$df <- out$k - 1
  out$z.value <- out$g/out$se.g
  out$p.value <- 2*pnorm(abs(out$z.value), lower.tail=FALSE)
  out$L.95ci <- out$g-1.96*out$se.g
  out$U.95ci  <- out$g+1.96*out$se.g
  out$p_homog <- ifelse(out$df==0,  1,  pchisq(out$Q, out$df, lower.tail=FALSE)) 
  out$I2 <- (out$Q-(out$df))/out$Q   #I-squared  
  out$I2 <- ifelse(out$I2<0, 0, out$I2)     
  out$I2 <- paste(round(out$I2*100, 4),  "%",  sep="")
  out$TW <- NULL
  out$TWD <- NULL
  out$TWDS <- NULL
  colnames(out) <- c("Mod", "k", "ES", "var.ES", "SE","Q", "df", "Z",
                     "p", "CI.lower", "CI.upper", "p.h",  "I2")
  out <- out[, c(1, 2, 3, 5, 4, 10, 11, 8, 9, 6,7, 12, 13 )] 
  return(out)
}


# Fixed effect single predictor categorical moderator Q-statistic

CatModfQ <- function(meta,  mod) {
  # Computes fixed effect Q-statistic (homogeneity test) for single predictor categorical 
  # moderator analysis. Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with g (standardized mean diff),var.g (variance of g),
  #         n.1 (grp 1 sample size), n.2 (grp 2 sample size).
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
  modsig <- CatModf(meta, mod)
  modsig$Mod <- as.factor(modsig$Mod)
  com1 <- levels(modsig$Mod)[x1]  # first level of moderator
  com2 <- levels(modsig$Mod)[x2]  # second level of moderator
  x1.es <- modsig[modsig$Mod==com1, "ES"]
  x2.es <- modsig[modsig$Mod==com2, "ES"]
  x1.var <- modsig[modsig$Mod==com1, "var.ES"]
  x2.var <- modsig[modsig$Mod==com2, "var.ES"]
  g <- (-1)*x1.es + 1*x2.es  # pg 288 (Cooper et al.,  2009)
  var <- (-1)^2*x1.var + (1)^2*x2.var
  df <- 2-1
  chi.sqr <- g^2/var
  m <- meta 
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  if(method == "post.hoc1") {  # post-hoc comparison (Tukey HSD method)
    fit<-aov(meta$g~meta$mod,weights=meta$wi)
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
  meta$k <- 1
  meta$TW <- 1/meta$var.g
  meta$TWS <- meta$TW^2
  meta$TWD <- meta$TW*meta$g
  meta$TWDS <- meta$TWD*meta$g
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
  meta$var.tau <- meta$var.g + tau
  meta$TW.tau <- 1/meta$var.tau
  meta$TWD.tau <- meta$TW.tau*meta$g
  meta$TWDS.tau <- meta$TWD.tau*meta$g
  out <- aggregate(meta[c("k","TW.tau","TWD.tau", 
                   "TWDS.tau")], by=meta["mod"], FUN=sum)
  lastrow <- dim(out)[1] + 1
  out[lastrow, 2:5] <- apply(out[,-1], 2, FUN=sum)
  out$mod[lastrow] <- "Overall"
  out$g <- out$TWD.tau/out$TW.tau
  out$var.g <- 1/out$TW.tau
  out$g <- out$TWD.tau/out$TW.tau
  out$se.g <- sqrt(out$var.g)
  out$Q <- out$TWDS.tau - out$TWD.tau^2/out$TW.tau
  out$df <- out$k - 1
  out$z.value <- out$g/out$se.g
  out$p.value <-    2*pnorm(abs(out$z.value), lower.tail=FALSE)
  out$L.95ci <- out$g-1.96*out$se.g
  out$U.95ci  <- out$g+1.96*out$se.g
  out$p_homog <- ifelse(out$df==0,  1,  pchisq(out$Q, out$df, lower.tail=FALSE)) 
  out$I2 <- (out$Q-(out$df))/out$Q   # these values are not to be used 
  out$I2 <- ifelse(out$I2<0, 0, out$I2)     
  out$I2 <- paste(round(out$I2*100, 4),  "%",  sep="")
  out$TW.tau <- NULL
  out$TWD.tau <- NULL
  out$TWDS.tau <- NULL
  colnames(out) <- c("Mod", "k", "ES", "var.ES", "SE","Q", "df", "Z",
                     "p", "CI.lower", "CI.upper", "p.h",  "I2")
  out <- out[, c(1, 2, 3, 5, 4, 10, 11, 8, 9, 6,7, 12, 13 )] 
  return(out)
}

#Q-statistic function (random effects)

CatModrQ <-  function(meta,  mod) {
  # Computes random effect Q-statistic (homogeneity test) for single predictor categorical 
  # moderator analysis. Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with g (standardized mean diff),var.g (variance of g),
  #         n.1 (grp 1 sample size), n.2 (grp 2 sample size).
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
  Qb <- Q-Qw  # overall between-group heterogeneity
  Q <- Qw + Qb
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
  #   meta: data.frame with g (standardized mean diff),var.g (variance of g),
  #         n.1 (grp 1 sample size), n.2 (grp 2 sample size).
  #   mod: Categorical moderator variable used for moderator analysis.    
  #   x1: One level of categorical moderator
  #   x2: Other level (comparison group) of same categorical moderator
  #   method: "post.hoc1" assumes the comparision was not planned prior to conducting
  #           the meta-analysis. The other option "planned" assumes you have planned 
  #           a priori to compare these levels of the categorical moderator. 
  #           Default is "post.hoc1". 
  # Returns:
  #   Random effects moderator means per group, mean difference (d), variance of difference,
  #   p-value, and 95% confidence intervals. 
  modsig <- CatModr(meta, mod)
  modsig$Mod <- as.factor(modsig$Mod)
  com1 <- levels(modsig$Mod)[x1]  # first level of moderator
  com2 <- levels(modsig$Mod)[x2]  # second level of moderator
  x1.es <- modsig[modsig$Mod==com1, "ES"]
  x2.es <- modsig[modsig$Mod==com2, "ES"]
  x1.var <- modsig[modsig$Mod==com1, "var.ES"]
  x2.var <- modsig[modsig$Mod==com2, "var.ES"]
  g <- (-1)*x1.es + 1*x2.es  # pg 288 (Cooper et al.,  2009)
  var <- (-1)^2*x1.var + (1)^2*x2.var
  df <- 2-1
  chi.sqr <- g^2/var
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  meta$wi <-  1/meta$var.g
  meta$wi2 <- (1/meta$var.g)^2
  meta$wiTi <- meta$wi*meta$g
  meta$wiTi2 <- meta$wi*(meta$g)^2
  meta$k <- meta$mod
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$g))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.g + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$g
  meta$wiTi2.tau <- meta$wi.tau*(meta$g)^2
  if(method == "post.hoc1") {  # post-hoc comparison (Tukey HSD method)
    fit<-aov(meta$g~meta$mod,weights=meta$wi.tau)
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

##==== META-REGRESSION FUNCTIONS (for continuous & categorical moderators)====##


# Assess the model fit from a meta-regression
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



       
##====== Multiple Moderator Graphs ======##


# MultiModGraph <- function(meta, conmod,  catmod, method="random",  
#                           conmod.name=NULL,  title=NULL) {
#   # Outputs a scatterplot and boxplot faceted by the categorical moderator from a 
#   # fixed or random effects moderator analysis. Computations derived from chapter 
#   # 14 and 15, Cooper et al. (2009). 
#   # Args:
#   #   meta: data.frame with id, g (standardized mean diff),var.g (variance of g).
#   #   conmod: Continuous moderator variable used for analysis.
#   #   catmod: Categorical moderator variable used for analysis.    
#   #   method: Model used, either "random" or "fixed" effects. Default is "random".
#   #   conmod.name: Name of continuous moderator variable to appear on x-axis of plot. 
#   #   Default is NULL.
#   #   title: Plot title. Default is NULL.
#   # Returns:
#   #   Multivariate moderator scatterplot graph faceted by categorical moderator levels. Also
#   #   places a weighted regression line based on either a fixed or random effects analysis. 
#   #   The ggplot2 packages outputs the rich graphics.  
#   m <- meta
#   m$conmod <- conmod
#   m$catmod <- catmod
#   compl <- !is.na(m$conmod)& !is.na(m$catmod)
#   meta <- m[compl, ]
#   meta$wi <-  1/meta$var.g
#   meta$wi2 <- (1/meta$var.g)^2
#   meta$wiTi <- meta$wi*meta$g
#   meta$wiTi2 <- meta$wi*(meta$g)^2
#   mod <- meta$mod
#   sum.wi <- sum(meta$wi, na.rm=TRUE)    
#   sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
#   sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
#   sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
#   comp <- sum.wi-sum.wi2/sum.wi
#   Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
#   k <- sum(!is.na(meta$g))                                  
#   df <- k-1      
#   tau <- (Q-(k - 1))/comp 				
#   meta$var.tau <- meta$var.g + tau
#   meta$wi.tau <- 1/meta$var.tau
#   meta$wiTi.tau <- meta$wi.tau*meta$g
#   meta$wiTi2.tau <- meta$wi.tau*(meta$g)^2
#   if(method=="fixed") {
#     multimod <- ggplot(meta, aes(conmod, g, weight=wi), na.rm=TRUE) + 
#                        opts(title=title, legend.position="none", na.rm=TRUE) + 
#                        facet_wrap(~catmod)  + 
#                        ggplot2::geom_point( aes(size=wi, shape=catmod)) + 
#                        geom_smooth(aes(group=1, weight=wi),
#                                    method= lm, se=FALSE, na.rm=TRUE) +
#                        ylab("Effect Size") + 
#                        xlab(conmod.name)
#   }  
#   if(method=="random") {
#     multimod <- ggplot(meta, aes(conmod, g, weight=wi.tau), na.rm=TRUE) + 
#                               opts(title=title, legend.position="none", na.rm=TRUE) + 
#                               facet_wrap(~catmod)  + 
#                               ggplot2::geom_point(aes(size=wi.tau, shape=catmod)) + 
#                               geom_smooth(aes(group=1, weight=wi.tau), 
#                                           method = lm, se = FALSE,  na.rm=TRUE) + 
#                               ylab("Effect Size") + 
#                               xlab(conmod.name)
#   }
#   return(multimod)
# }


# Correction for Attenuation

Rho_TU<- function(g,xx,yy) {
  g.corrected<-g/(sqrt(xx)*sqrt(yy))
  return(g.corrected)
   }


CorAtten <- function (meta, xx, yy) {
    m <- meta
    meta$xx <- xx
    meta$yy <- yy
    meta$g <- ifelse(is.na(meta$xx & meta$yy), meta$g, Rho_TU(meta$g, 
        meta$xx, meta$yy))
    meta$wi <- 1/meta$var.g
    meta$wiTi <- meta$wi * meta$g
    meta$wiTi2 <- meta$wi * (meta$g)^2
    sum.wi <- sum(meta$wi, na.rm = TRUE)
    sum.wi2 <- sum(meta$wi2, na.rm = TRUE)
    sum.wiTi <- sum(meta$wiTi, na.rm = TRUE)
    sum.wiTi2 <- sum(meta$wiTi2, na.rm = TRUE)
    comp <- sum.wi - sum.wi2/sum.wi
    Q <- sum.wiTi2 - (sum.wiTi^2)/sum.wi
    k <- sum(!is.na(meta$g))
    df <- k - 1
    tau <- (Q - k + 1)/comp
    meta$var.tau <- meta$var.g + tau
    meta$wi.tau <- 1/meta$var.tau
    meta$wiTi.tau <- meta$wi.tau * meta$g
    meta$wiTi2.tau <- meta$wi.tau * (meta$g)^2
    return(meta)
}

