wilcox.signtest <- function (x,...) {
  UseMethod("wilcox.signtest")
}

wilcox.signtest.formula <- function(formula,data,subset,...) {
  if (missing(formula) || (length(formula)!=3L) || (length(attr(terms(formula[-2L]), 
      "term.labels"))!=1L)) {stop("'formula' missing or incorrect")}
  m <- match.call(expand.dots=FALSE)
  if (is.matrix(eval(m$data,parent.frame())))  {m$data <- as.data.frame(data)}
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m,parent.frame())
  DNAME <- paste(names(mf),collapse=" by ")
  names(mf) <- NULL
  response <- attr(attr(mf,"terms"),"response")
  g <- factor(mf[[-response]])
  if (nlevels(g)!=2L) {stop("grouping factor must have exactly 2 levels")}
  DATA <- setNames(split(mf[[response]],g),c("x","y"))
  res <- do.call("wilcox.signtest",c(DATA,list(...)))
  res$data.name <- DNAME
  return(res)
}

wilcox.signtest.default <- function(x,y=NULL,mu=0,conf.level=0.95,...) {
  if (!is.numeric(x)) {stop("'x' must be numeric")}
  if (!is.numeric(mu)) {stop("'mu' must be numeric")}
  if (!(conf.level>0 & conf.level<1)) {stop("'conf.level' must be a number between 0 and 1")}
  names(mu) <- est.name <- ifelse(is.null(y),"location","location shift")
  if (is.null(y)) {
    dname <- deparse(substitute(x))
  } else {
    dname <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
    OK <- complete.cases(x,y)
    if (length(x[OK])!=length(y[OK])) {stop("'x' and 'y' lengths differ")}
    x <- x[OK]-y[OK]
    y <- NULL
  }
  signs <- x-mu
  signs <- signs[signs!=0]
  if (length(signs)<length(x)) {warning("zeroes are present")}
  pval <- binom.test(length(signs[signs>0]),length(signs),p=0.5)$p.value
  x <- signs+mu
  n <- length(x)
  ci.all <- 1-2*pbinom(0:(n/2),n,0.5)
  ci.ranks <- which.min(abs(ci.all-conf.level))
  ci <- sort(x)[c(ci.ranks,n-ci.ranks+1)]
  achieved.ci <- ci.all[ci.ranks]
  if (abs(conf.level-achieved.ci)>0.1) {
    warning("requested conf.level not achievable")
    conf.level <- signif(achieved.ci, 2)
  }
  attr(ci,"conf.level") <- conf.level
  estimate <- median(x)
  names(estimate) <- est.name
  res <- list(method="Wilcoxon sign test",data.name=dname,null.value=mu,p.value=pval,
    alternative="two.sided",estimate=estimate,conf.int=ci)
  class(res) <- "htest"
  return(res)
}
