"oddsratio.wald" <-
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
  wald <- matrix(NA, nr, 3)
  wald[1,1] <- 1  
  for(i in 2:nr){
    a0<-x[1,2]; b0<-x[1,1]; a1<-x[i,2]; b1<-x[i,1]
    est <- (b0*a1)/(a0*b1)
    logOR <- log(est)
    SElogOR <- sqrt((1/b0)+(1/a0)+(1/b1)+(1/a1))
    ci <- exp(logOR + c(-1, 1)*Z*SElogOR)
    wald[i,] <- c(est, ci)
  }
  pv <- tab2by2.test(x, correction = correction)
  colnames(wald) <- c("estimate", "lower", "upper")
  rownames(wald) <- rownames(x)
  cn2 <- paste("odds ratio with",
               paste(100*conf.level, "%", sep=""),
               "C.I.")  
  names(dimnames(wald)) <- c(names(dimnames(x))[1], cn2)
  rr <- list(x = x,
             data = tmx,
             p.exposed = p.exposed,
             p.outcome = p.outcome,
             measure = wald,
             conf.level = conf.level,
             p.value = pv$p.value,
             correction = pv$correction             
             )
  rrs <- list(data = tmx,
               measure = wald,
               p.value = pv$p.value,
               correction = pv$correction
               )  
  attr(rr, "method") <- "Unconditional MLE & normal approximation (Wald) CI"
  attr(rrs, "method") <- "Unconditional MLE & normal approximation (Wald) CI"
  if(verbose==FALSE){
    rrs
  } else rr
}
