confint.mcprofile <-
function(object, parm, level=0.95, adjust=c("single-step","none","bonferroni"), alternative=c("two.sided","less","greater"), ...){
  pam <- c("bonferroni", "none", "single-step")
  if (!(adjust[1] %in% pam)) stop(paste("adjust has to be one of:", paste(pam, collapse=", ")))
  CM <- object$CM
  df <- object$df
  cr <- NULL
  sdlist <- object$srdp  
  spl <- lapply(sdlist, function(x){
    x <- na.omit(x)
    try(interpSpline(x[,1], x[,2]))
  })
  
  if (adjust[1] == "none" | nrow(CM) == 1){
    if (alternative[1] == "two.sided") alpha <- (1-level)/2 else alpha <- 1-level
    if (is.null(df)) quant <- qnorm(1-alpha) else quant <- qt(1-alpha, df=df)
  }
  if (adjust[1] == "bonferroni"){
    if (alternative[1] == "two.sided") alpha <- (1-level)/2 else alpha <- 1-level
    if (is.null(df)) quant <- qnorm(1-alpha/nrow(CM)) else quant <- qt(1-alpha/nrow(CM), df=df)
  }
  if (adjust[1] == "single-step" & nrow(CM) > 1){
    vc <- vcov(object$object)
    VC <- CM %*% vc %*% t(CM)
    d <- 1/sqrt(diag(VC))
    dd <- diag(d)
    cr <- dd %*% VC %*% dd
    if (alternative[1] == "two.sided"){
      if (is.null(df)) quant <- qmvnorm(level, corr=cr, tail="both.tails")$quantile else quant <- qmvt(level, df=df, corr=cr, tail="both.tails")$quantile
    } else {
      if (is.null(df)) quant <- qmvnorm(level, corr=cr, tail="lower.tail")$quantile else quant <- qmvt(level, df=df, corr=cr, tail="lower.tail")$quantile
    }
  }
   
  if (alternative[1] == "two.sided"){
    ci <- data.frame(t(sapply(spl, function(x, quant){
      pfun <- function(xc, obj, quant) predict(obj, xc)$y-quant
      upper <- try(uniroot(pfun, range(predict(x)$x), obj=x, quant=quant)$root, silent=TRUE)
      if (class(upper)[1] == "try-error") upper <- NA
      lower <- try(uniroot(pfun, range(predict(x)$x), obj=x, quant=-quant)$root, silent=TRUE)
      if (class(lower)[1] == "try-error") lower <- NA
      c(lower, upper)
    }, quant=quant)))
    names(ci) <- c("lower", "upper")
  }
  if (alternative[1] == "less"){
    ci <- data.frame(sapply(spl, function(x, quant){
      pfun <- function(xc, obj, quant) predict(obj, xc)$y-quant
      upper <- try(uniroot(pfun, range(predict(x)$x), obj=x, quant=quant)$root, silent=TRUE)
      if (class(upper)[1] == "try-error") upper <- NA
      cbind(c(upper))
    }, quant=quant))
    names(ci) <- "upper"
  }
  if (alternative[1] == "greater"){
    ci <- data.frame(sapply(spl, function(x, quant){
      pfun <- function(xc, obj, quant) predict(obj, xc)$y-quant
      lower <- try(uniroot(pfun, range(predict(x)$x), obj=x, quant=-quant)$root, silent=TRUE)
      if (class(lower)[1] == "try-error") lower <- NA
      cbind(c(lower))
    }, quant=quant))
    names(ci) <- "lower"
  }  

  out <- list()
  out$estimate <- data.frame(Estimate=CM %*% coefficients(object$object))
  out$confint <- ci
  out$cr <- cr
  out$CM <- CM
  out$quant <- quant
  out$alternative <- alternative[1]
  out$level <- level
  out$adjust <- adjust[1]
  class(out) <- "mcpCI"
  return(out)
}

