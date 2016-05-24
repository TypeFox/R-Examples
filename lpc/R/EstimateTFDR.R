EstimateTFDR <- function(x,y,type,censoring.status=NULL){
  call <- match.call()
  dat <- list(x=x,y=y,censoring.status=censoring.status)
  CheckEstimateTFDRFormat(dat,type)
  ttstar <- NULL
  if(type=="regression"){
    tt <- quantitative.func(dat$x, dat$y, .05)$tt
    for (i in 1:100) ttstar <- c(ttstar, quantitative.func(dat$x, sample(dat$y), .05)$tt)
  } else if (type=="multiclass"){
    tt <- multiclass.func(dat$x, dat$y, .05)$tt
    for(i in 1:100) ttstar <- c(ttstar, multiclass.func(dat$x, sample(dat$y), .05)$tt)
  } else if(type=="two class"){
    tt <- ttest.func(dat$x, dat$y, .05)$tt
    for (i in 1:100) ttstar <- c(ttstar, ttest.func(dat$x, sample(dat$y), .05)$tt)
  } else if(type=="survival"){
    tt <- cox.func(dat$x, dat$y, dat$censoring.status, .05)$tt
    for (i in 1:100){
      oo <- sample(ncol(dat$x))
      ttstar <- c(ttstar, cox.func(dat$x, dat$y[oo], dat$censoring.status[oo], .05)$tt)
    }  
  }
  pi0 <- min(1, sum(abs(tt)>quantile(abs(ttstar), 0) & abs(tt)<quantile(abs(ttstar), .5))/(0.5*length(tt)))
  fdr <- NULL
  for(q in 1:nrow(dat$x)){
    fdr <- c(fdr, sum(abs(ttstar) >= abs(tt[q]))/((length(ttstar)/length(tt))*sum(abs(tt) >= abs(tt[q]))))
  }
  tfdrobj <- list(fdrt=pmin(1,pi0*fdr), pi0=pi0,call=call)
  class(tfdrobj) <- "tfdrobj"
  return(tfdrobj)
}  

