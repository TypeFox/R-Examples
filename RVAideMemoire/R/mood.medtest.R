mood.medtest <- function(x,...) {
  UseMethod("mood.medtest")
}

mood.medtest.formula <- function(formula,data,subset,...) {
  if (missing(formula) || (length(formula)!=3L)) {stop("'formula' missing or incorrect")}
  m <- match.call(expand.dots=FALSE)
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(data)}
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m,parent.frame())
  if (length(mf)>2L) {stop("'formula' should be of the form response ~ group")}
  DNAME <- paste(names(mf),collapse=" by ")
  names(mf) <- NULL
  response <- attr(attr(mf,"terms"),"response")
  x <- mf[[response]]
  g <- factor(mf[[-response]])
  res <- do.call("mood.medtest",c(list(x,g),list(...)))
  res$data.name <- DNAME
  return(res)
}

mood.medtest.default <- function(x,g,exact=NULL,...) {
  if (!is.factor(g)) {g <- factor(g)}
  OK <- complete.cases(x,g)
  if (length(x[OK])!=length(g[OK])) {stop("'x' and 'g' lengths differ")}
  dname <- paste(deparse(substitute(x)),"by",deparse(substitute(g)))
  x <- x[OK]
  g <- g[OK]
  ng <- table(g)
  if (any(ng<2)) {
    ng.to.del <- names(ng)[which(ng<2)]
    x <- x[g!=ng.to.del]
    g <- droplevels(g[g!=ng.to.del])
  }
  med <- median(x)
  cont <- table(x>med,g)
  if(is.null(exact)) {exact <- sum(cont)<200}
  res <- list(method="Mood's median test",data.name=dname)
  if (exact) {
    test <- fisher.test(cont)
    res$p.value <- test$p.value
  } else {
    test <- suppressWarnings(chisq.test(cont))
    res$statistic <- test$statistic
    res$parameter <- test$parameter
    res$p.value <- test$p.value
  }
  class(res) <- "htest"
  return(res)
}
