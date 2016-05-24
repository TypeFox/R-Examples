nca.complete <- function(conc, time, n.tail=3, dose=0, method=c("z", "boott"), conf.level=0.95, 
     nsample=1000, data) {

  if(!missing(data)){
    cnames <- colnames(data)
    if(!any(cnames=='conc')){stop("data does not contain a variable conc")}
    if(!any(cnames=='time')){stop("data does not contain a variable time")}
    conc <- data$conc
    time <- data$time
  }

  method <- "z"
  dat <- cbind(conc,time)
  dat <- dat[order(time),]
  n <- as.numeric(table(dat[,"time"])[1])
  tp <- nrow(dat)/n
  conc.list <-  split(dat[,"conc"],rep(1:n,tp))
  
  ests <- matrix(NA,nrow=7,ncol=n)
  for (i in 1:n) {

    res <- nca.ssd(conc=rep(conc.list[[i]],2), time=rep(unique(dat[,"time"]),2),  n.tail=n.tail, dose=dose, method=method, conf.level=conf.level)
    ests[,i] <- res$est
  }
  

  ## return geometric mean of parameters + message
  if(n>1) {
    message("Note that the geometric mean of the individual parameters is reported.\n") 
  }
  res$est <- matrix(apply(ests,1,prod)^(1/n),ncol=1)
  rownames(res$est) <- c('AUC to tlast', 'AUC to infinity', 'AUMC to infinity', 'Mean residence time', 'non-compartmental half-life', 'Clearance', 'Volume of distribution at steady state') 
  res$CIs[,"est"]<-  apply(ests,1,prod)^(1/n)
  ## set CIs NA as complete data design is used
  res$CIs[,c("lower","upper")]<- NA
  res$CIs[,"stderr"] <- NA
  res$design <- "complete"
  return(res)
}
