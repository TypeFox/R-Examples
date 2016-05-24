"rateratio.midp" <-
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
  alpha <- 1 - conf.level
  nr <- nrow(x)
  est <- matrix(NA, nr, 3)
  est[1,1] <- 1
  for(i in 2:nr){
    aa <- x[i,1]; bb <- x[1,1]; pt1 <- x[i,2]; pt0 <- x[1,2]
    pt <- pt0 + pt1
    mm <- aa + bb
    s.irr <- uniroot(function(s) (dbinom(aa, mm, s)/2 + pbinom(aa-1, mm, s)) - 0.5, c(0, 1))$root
    irr <- (s.irr / (1 - s.irr)) * (pt0 / pt1)
    s.lower <- uniroot(function(s) 1 - (dbinom(aa, mm, s)/2 + pbinom(aa-1, mm, s)) - alpha/2, c(0,1))$root
    irr.lower <- (s.lower/ (1 - s.lower)) * (pt0 / pt1)
    s.upper <- uniroot(function(s) dbinom(aa, mm, s)/2 + pbinom(aa-1, mm, s) - alpha/2, c(0,1))$root
    irr.upper <- (s.upper/ (1 - s.upper)) * (pt0 / pt1)
    est[i,] <- c(irr, irr.lower, irr.upper)
  }
  pval <- rate2by2.test(x)$p.value
  colnames(est) <- c("estimate", "lower", "upper")
  rownames(est) <- rownames(x)
  cn2 <- paste("rate ratio with",
               paste(100*conf.level, "%", sep=""),
               "C.I.")  
  names(dimnames(est)) <- c(names(dimnames(x))[1], cn2)
  rr <- list(x = x,
             data = tmx,
             measure = est,
             conf.level = conf.level,
             p.value = pval
             )
  rrs <- list(data = tmx,
              measure = est,
              p.value = pval
              )
  attr(rr, "method") <- "Median unbiased estimate & mid-p exact CI"
  attr(rrs, "method") <- "Median unbiased estimate & mid-p exact CI"
  if(verbose==FALSE) {
    rrs
  } else rr 
}
