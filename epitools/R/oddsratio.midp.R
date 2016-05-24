"oddsratio.midp" <-
function(x, y = NULL,
         conf.level = 0.95,
         rev = c("neither", "rows", "columns", "both"),
         correction = FALSE,
         verbose = FALSE,
         interval =  c(0, 1000)){
  if(is.matrix(x) && !is.null(y)){stop("y argument should be NULL")}
  if(is.null(y)){
    x <- epitable(x, rev = rev)
  } else {
    x <- epitable(x, y, rev = rev)
  }
  tmx <- table.margins(x)
  p.exposed <- sweep(tmx,2,tmx["Total",],"/")
  p.outcome <- sweep(tmx,1,tmx[,"Total"],"/")
  nr <- nrow(x)
  midp <- matrix(NA, nr, 3)
  midp[1,1] <- 1  
  for(i in 2:nr){
    a0<-x[1,2]; b0<-x[1,1]; a1<-x[i,2]; b1<-x[i,1]
    tmpx <- matrix(c(a1,a0,b1,b0),2,2, byrow=TRUE)
    OR <- or.midp(tmpx, conf.level = conf.level, interval = interval)
    midp[i,] <- c(OR$estimate, OR$conf.int)
  }
  pv <- tab2by2.test(x, correction = correction)
  colnames(midp) <- c("estimate", "lower", "upper")
  rownames(midp) <- rownames(x)
  cn2 <- paste("odds ratio with",
               paste(100*conf.level, "%", sep=""),
               "C.I.")  
  names(dimnames(midp)) <- c(names(dimnames(x))[1], cn2)
  rr <- list(x = x,
             data = tmx,
             p.exposed = p.exposed,
             p.outcome = p.outcome,
             measure = midp,
             conf.level = conf.level,
             p.value = pv$p.value,
             correction = pv$correction             
             )
  rrs <- list(data = tmx,
              measure = midp,
              p.value = pv$p.value,
              correction = pv$correction
              )  
  attr(rr, "method") <- "median-unbiased estimate & mid-p exact CI"
  attr(rrs, "method") <- "median-unbiased estimate & mid-p exact CI"
  if(verbose==FALSE){
    rrs
  } else rr
}
