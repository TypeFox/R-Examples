exp.mcpCI <-
function(x){
  x$estimate <- exp(x$estimate)
  x$confint <- exp(x$confint)
  return(x)
}
expit <- function(x) 1/(1 + exp(-x))
expit.mcpCI <-
function(x){
  expit <- function(x) 1/(1 + exp(-x))
  x$estimate <- expit(x$estimate)
  x$confint <- expit(x$confint)
  return(x)
}
plot.mcpCI <-
function(x, ...){
  Estimate <- x$estimate
  CM <- x$CM
  grp <- factor(rownames(CM), levels=rownames(CM)[nrow(CM):1])
  cidat <- data.frame(grp, Estimate)
  switch(x$alternative,
         two.sided = {
           lower <- x$confint$lower
           upper <- x$confint$upper
         }, greater = {
           lower <- x$confint$lower
           upper <- Inf
         }, lower = {
           lower <- x$confint$lower
           upper <- Inf
         })
  lower[is.na(lower)] <- -Inf
  upper[is.na(upper)] <- Inf
  cidat$lower <- lower
  cidat$upper <- upper
  ggplot(cidat, aes(x=Estimate, y=grp, xmin=lower, xmax=upper)) + geom_errorbarh(height=0.3) + geom_point() + ylab("") + xlab("")
}

plot.mcprofile <- function(x, ...){
  sdlist <- x$srdp  
  CM <- x$CM
  nr <- sapply(sdlist, nrow)
  grp <- factor(rep(rownames(CM), nr), levels=rownames(CM))
  b1 <- unlist(lapply(sdlist, function(x) x[,1]))
  z1 <- unlist(lapply(sdlist, function(x) x[,2]))
  sdd <- data.frame(b=b1, z=z1, grp)

  spl <- lapply(sdlist, function(x){
    x <- na.omit(x)
    sp <- try(interpSpline(x[,1], x[,2]))
    bcoord <- seq(min(x[,1]), max(x[,1]), length=200)
    preds <- if (class(sp)[1] == "try-error") y <- rep(NA, 200) else predict(sp, x=bcoord)$y
    data.frame(b=bcoord, z=preds)
  })
  pgrp <- factor(rep(rownames(CM), each=200), levels=rownames(CM))
  b <- unlist(lapply(spl, function(x) x[,1]))
  z <- unlist(lapply(spl, function(x) x[,2]))
  spd <- data.frame(b, z, grp=pgrp)
  
  ggplot(spd, aes(x=b, y=z)) + geom_line() + facet_wrap(~grp, scales = "free_x") + geom_hline(yintercept=0, linetype=2) + geom_point(data=sdd, shape="|")
}

print.mcpCI <-
function(x, ...){
  cat("\n   mcprofile - Confidence Intervals \n\n")
  cat("level: \t\t", x$level,"\n")
  cat("adjustment:\t", x$adjust,"\n\n")
  print(data.frame(Estimate=x$estimate, x$confint), digits=3)
  cat("\n")
}
print.mcprofile <-
function(x, ...){
  cat("\n   Multiple Contrast Profiles\n\n")
  est <- x$CM %*% coefficients(x$object)
  vest <- x$CM %*% vcov(x$object) %*% t(x$CM)
  sdest <- sqrt(diag(vest))
  if (is.null(x$adjestimates)){
    print(data.frame(Estimate=est, Std.err=sdest), digits=3)
  } else {
    print(data.frame(Estimate=est, Std.err=sdest, adjEstimate=x$adjestimates), digits=3)
  }
  cat("\n")
}
print.mcpSummary <-
function(x, ...){
  cat("\n   mcprofile - Multiple Testing\n\n")
  cat("Adjustment:\t", x$adjust, "\n")
  cat("Margin:  \t", x$margin, "\n")
  cat("Alternative:\t", x$alternative, "\n\n")
  dat <- data.frame(round(x$CM %*% x$estimate,2), round(x$statistic,2), round(x$p.values,3))
  pname <- switch(x$alternative, less = paste("Pr(<", ifelse(is.null(x$df), "z", "t"), ")", sep = ""), greater = paste("Pr(>", ifelse(is.null(x$df), "z", "t"), ")", sep = ""), two.sided = paste("Pr(>|", ifelse(is.null(x$df), "z", "t"), "|)", sep = ""))
  names(dat) <- c("Estimate", "Statistic", pname)
  alt <- switch(x$alternative, two.sided = "==", less = ">=", greater = "<=")
  rownames(dat) <- paste(rownames(dat), alt, x$margin)
  error <- attr(x$p.values, "error")
  if (!is.null(error) && error > .Machine$double.eps) {
    sig <- which.min(abs(1/error - (10^(1:10))))
    sig <- 1/(10^sig)
  } else {
    sig <- .Machine$double.eps
  }
  printCoefmat(dat, digits = max(3, getOption("digits") - 3), has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
  cat("\n")
}
wald <-
function(object){
  srdp <- object$srdp
  est <- object$CM %*% coefficients(object$object)
  sde <- sqrt(diag(object$CM %*% vcov(object$object) %*% t(object$CM)))
  wsrdp <- lapply(1:length(srdp), function(i){
    srdpi <- srdp[[i]]
    b <- srdpi[,1]
    srdpi[,2] <- (b-est[i])/sde[i]
    srdpi 
  })
  object$srdp <- wsrdp
  return(object)
}
