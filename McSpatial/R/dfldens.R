dfldens <- function(y,lgtform,window=0,bandwidth=0,kern="tcub",probit=FALSE,graph=TRUE,yname="y",alldata=FALSE,data=NULL) {

  xmat <- model.frame(lgtform,data=data)
  datevar <- xmat[,1]
  n = length(datevar)
  n1 = sum(datevar)
  n0 = n-n1
  p0 = n0/n
  linktype = "logit"
  if (probit==TRUE) {linktype=="probit"}
  fit <- glm(lgtform,data=data,family=binomial(link=linktype))
  px0 <- 1-fitted(fit)
  tau <- (px0/(1-px0))/(p0/(1-p0))
  tau <- tau[datevar==1]
  tau <- tau/mean(tau)

  if (kern=="rect")  { wgt <- function(psi) {ifelse(abs(psi)>=0,.5,0) } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { .75*(1-psi^2) } }
  if (kern=="bisq")  { wgt <- function(psi) { (15/16)*((1-psi^2)^2) } }
  if (kern=="tcub")  { wgt <- function(psi) { (70/81)*((1 - abs(psi)^3)^3) } }
  if (kern=="trwt")  { wgt <- function(psi) { (35/32)*((1 - psi^2)^3) } }
  if (kern=="gauss") { wgt <- function(psi) { dnorm(psi) } }

  y1 <- y[datevar==1]

  if (alldata==TRUE){target <- y1}
  if (alldata==FALSE){
    if (window==0&bandwidth==0) {bandwidth = (.9*(quantile(y1,.75)-quantile(y1,.25))/1.34)*(n1^(-.20)) }
    if (window>0)    {target <- lfeval(locfit(~lp(y1,nn=window,deg=0),kern=kern))$xev }
    if (bandwidth>0) {target <- lfeval(locfit(~lp(y1,h=bandwidth,deg=0),kern=kern))$xev }
  }

  nt = length(target)
  dtarget1 <- array(0,dim=nt)
  dtarget10 <- array(0,dim=nt)
  for (i in seq(1,nt)) {
    dist <- abs(y1 - target[i])
    if (window>0) {h = quantile(dist,window)}
    if (bandwidth>0) {h = bandwidth}
    sampvar <- dist<=h
    if (kern=="gauss") {sampvar <- dist<=max(dist)}
    dtarget1[i] <-  sum( wgt(dist[sampvar]/h) )/(h*n1)
    dtarget10[i] <- sum( wgt(dist[sampvar]/h)*tau[sampvar] )/(h*n1)
  }
  if (alldata==TRUE) {
    dhat1 <- dtarget1
    dhat10 <- dtarget10
  }
  if (alldata==FALSE) {
    dhat1 <- smooth12(target,dtarget1,y1)
    dhat10 <- smooth12(target,dtarget10,y1)
  }
 
  if (graph==TRUE) {
    o <- order(y1)
    ymin = min(dhat1,dhat10)
    ymax = max(dhat1,dhat10)
    plot(y1[o],dhat1[o],xlab=yname,ylab="Density",type="l",ylim=c(ymin,ymax))
    lines(y1[o],dhat10[o],col="red")
    legend("topright",c("Actual","Counterfactual"),col=c("black","red"),lwd=1)
  }

  out <- list(target,dtarget1,dtarget10,dhat1,dhat10)
  names(out) <- c("target","dtarget1","dtarget10","dhat1","dhat10")
  return(out)
}

