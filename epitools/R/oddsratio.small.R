"oddsratio.small" <-
function(x, y = NULL,
         conf.level = 0.95,
         rev = c("neither", "rows", "columns", "both"),
         correction = FALSE,
         verbose = FALSE){
  if(is.matrix(x) && !is.null(y)){stop("y argument should be NULL")}
  if(is.null(y)){
    x <- epitable(x, rev = rev)
  } else {
    x <- epitable(x, y, rev = rev)
  }
  tmx <- table.margins(x)
  p.exposed <- sweep(tmx,2,tmx["Total",],"/")
  p.outcome <- sweep(tmx,1,tmx[,"Total"],"/")  
  Z <- qnorm(0.5*(1 + conf.level))
  nr <- nrow(x)
  small <- matrix(NA, nr, 3)
  or <- rep(NA, nr)
  small[1,1] <- 1
  for(i in 2:nr){
    a0<-x[1,2]; b0<-x[1,1]; a1<-x[i,2]; b1<-x[i,1]
    or[i] <- (b0*a1)/(a0*b1)
    est <- (b0*a1)/((a0+1)*(b1+1))
    logORss <- log(((b0+0.5)*(a1+0.5))/((a0+0.5)*(b1+0.5)))
    SElogORss <- sqrt((1/(b0+0.5))+(1/(a0+0.5))+(1/(b1+0.5))+(1/(a1+0.5)))
    ci <- exp(logORss + c(-1, 1)*Z*SElogORss)    
    small[i,] <- c(est, ci)
  }
  if(any(or, na.rm=TRUE)<1){
    cat("CAUTION: At least one unadjusted odds ratio < 1.
Do not use small sample-adjusted OR to esimate 1/OR.",fill=1)
  }
  pv <- tab2by2.test(x, correction = correction)
  colnames(small) <- c("estimate", "lower", "upper")
  rownames(small) <- rownames(x)
  cn2 <- paste("odds ratio with",
               paste(100*conf.level, "%", sep=""),
               "C.I.")  
  names(dimnames(small)) <- c(names(dimnames(x))[1], cn2)
  rr <- list(x = x,
             data = tmx,
             p.exposed = p.exposed,
             p.outcome = p.outcome,
             measure = small,
             conf.level = conf.level,
             p.value = pv$p.value,
             correction = pv$correction             
             )
  rrs <- list(data = tmx,
               measure = small,
               p.value = pv$p.value,
               correction = pv$correction
               )  
  attr(rr, "method") <- "small sample-adjusted UMLE & normal approx (Wald) CI"
  attr(rrs, "method") <- "small sample-adjusted UMLE & normal approx (Wald) CI"
  if(verbose==FALSE){
    rrs
  } else rr
}
