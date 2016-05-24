hmftest <- function(x,...){
  UseMethod("hmftest")
}

hmftest.formula <- function(x, alt.subset, ...){
  formula <- x
  x <- mlogit(formula,...)
  x$call$data <- match.call()$data
  xs <- mlogit(formula, alt.subset=alt.subset, ...)
  hmftest(x,xs)
}

hmftest.mlogit <- function(x, z, ...){
  if (is.character(z)) xs <- update(x,alt.subset=z)
  if (class(z)=="mlogit") xs <- z
  coef.x <- coef(x)
  coef.s <- coef(xs)
  un <- names(coef.x) %in% names(coef.s)
  diff.coef <- coef.s-coef.x[un]
  diff.var <- vcov(xs)-vcov(x)[un,un]
  hmf <- as.numeric(diff.coef%*%solve(diff.var)%*%diff.coef)
  names(hmf) <- "chisq"
  df <- sum(un)
  names(df) <- "df"
  pv <- pchisq(hmf,df=df,lower.tail=FALSE)
  res <- list(data.name = x$call$data,
              statistic = hmf,
              p.value =pv,
              parameter = df,
              method = "Hausman-McFadden test",
              alternative = "IIA is rejected")
  class(res) <- "htest"
  res
}  



mfR2 <- function(x){
##   ll <- logLik(x)
##   data.name <- x$call$data
##   choice.name <- as.character(x$call$formula[[2]])
##   data <- eval(data.name,envir=parent.frame())
##   alt <- data[[2]]
##   choice <- data[[choice.name]]
##   eff <- table(alt[choice])
##   n <- sum(eff)
##   llo <- sum(eff*log(eff/n))
##   1-ll/llo
  logLik0 <- attr(x$logLik, 'null')
  1-x$logLik/logLik0
}

lratio <- function(object){
  freq <- object$freq
  llo <- sum(freq*log(prop.table(freq)))
  data.name <- object$call$data
  stat <- -2*(llo-logLik(object))
  names(stat) <- "chisq"
  parameter <- length(coef(object))-length(freq)+1
  names(parameter) <- "df"
  pval <- pchisq(stat,df=parameter,lower.tail=FALSE)
  lrtest <- list(statistic = stat,
                 data.name = data.name,
                 p.value = pval,
                 parameter = parameter,
                 method = "likelihood ratio test")
  class(lrtest) <- "htest"
  lrtest
}
