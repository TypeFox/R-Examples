"riskratio.boot" <-
function(x, y = NULL,
         conf.level = 0.95,
         rev = c("neither", "rows", "columns", "both"),
         correction = FALSE,
         verbose = FALSE,
         replicates = 5000){
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
  boot <- matrix(NA, nr, 3)
  boot[1,1] <- 1  
  rr.boot <-
    function(a1, a0, b1, b0, conf.level = 0.95,
             replicates = 5000){
      alpha <- 1 - conf.level
      n1 <- a1 + b1; n0 <- a0 + b0
      p1 <- a1/n1; p0 <- a0/n0
      r1 <- rbinom(replicates, n1, p1)/n1
      r0 <- rbinom(replicates, n0, p0)/n0
      rrboot <- r1/r0
      rrbar <- mean(rrboot)
      ci <- quantile(rrboot, c(alpha/2, 1-alpha/2))
      list(p0 = p0, p1 = p1, rr = p1/p0, rr.mean = rrbar,
           conf.level = conf.level, conf.int = unname(ci),
           replicates = replicates)
    }  
  for(i in 2:nr){
    a0<-x[1,2]; b0<-x[1,1]; a1<-x[i,2]; b1<-x[i,1]
    rrb <- rr.boot(a1 = a1, a0 = a0, b1 = b1, b0 = b0,
                   conf.level = conf.level,
                   replicates = replicates)
    boot[i,] <- c(rrb$rr, rrb$conf.int)
  }
  pv <- tab2by2.test(x, correction = correction)
  colnames(boot) <- c("estimate", "lower", "upper")
  rownames(boot) <- rownames(x)
  cn2 <- paste("risk ratio with",
               paste(100*conf.level, "%", sep=""),
               "C.I.")  
  names(dimnames(boot)) <- c(names(dimnames(x))[1], cn2)
  rr <- list(x = x,
             data = tmx,
             p.exposed = p.exposed,
             p.outcome = p.outcome,
             measure = boot,
             replicates = rrb$replicates,
             p.value = pv$p.value,
             correction = pv$correction             
             )
  rrs <- list(data = tmx,
              measure = boot,
              p.value = pv$p.value,
              correction = pv$correction
              )
  attr(rr, "method") <- "Unconditional MLE & bootstrap CI"
  attr(rrs, "method") <- "Unconditional MLE & bootstrap CI"
  if(verbose==FALSE) {
    rrs
  } else rr 
}
