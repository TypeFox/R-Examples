"rateratio.wald" <-
function(x, y = NULL,
         conf.level = 0.95,
         rev = c("neither", "rows", "columns", "both"),
         verbose = FALSE){
  if(is.matrix(x) && !is.null(y)){stop("y argument should be NULL")}
  if(is.null(y)){
    x <- ratetable(x, rev = rev)
  } else {
    xn <- substitute(x)
    yn <- substitute(y)
    x <- ratetable(x, y, rev = rev)
    colnames(x) <- c(xn, yn)
  }
  tmx <- table.margins(x)[,-3]
  Z <- qnorm(0.5*(1 + conf.level))
  nr <- nrow(x)
  wald <- matrix(NA, nr, 3)
  pval <- matrix(NA, nr, 2)
  wald[1,1] <- 1  
  for(i in 2:nr){
    aa <- x[i,1]; bb <- x[1,1]; pt1 <- x[i,2]; pt0 <- x[1,2]
    est <- (aa/pt1)/(bb/pt0)
    logRR <- log(est)
    SElogRR <- sqrt((1/aa) + (1/bb))
    ci <- exp(logRR + c(-1, 1)*Z*SElogRR)
    wald[i,] <- c(est, ci)
  }
  pval <- rate2by2.test(x)$p.value
  colnames(wald) <- c("estimate", "lower", "upper")
  rownames(wald) <- rownames(x)
  cn2 <- paste("rate ratio with",
               paste(100*conf.level, "%", sep=""),
               "C.I.")  
  names(dimnames(wald)) <- c(names(dimnames(x))[1], cn2)
  rr <- list(x = x,
             data = tmx,
             measure = wald,
             conf.level = conf.level,
             p.value = pval
             )
  rrs <- list(data = tmx,
              measure = wald,
              p.value = pval
              )
  attr(rr, "method") <- "Unconditional MLE & normal approximation (Wald) CI"
  attr(rrs, "method") <- "Unconditional MLE & normal approximation (Wald) CI"
  if(verbose==FALSE) {
    rrs
  } else rr 
}
