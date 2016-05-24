ksim <- function(long1,lat1,long2,lat2,kilometer=FALSE,noplot=FALSE,dmin=0,dmax=0,dlength=512,h=0,kern="gaussian",
  nsim=2000,nsamp=0,pval=.05,cglobal=FALSE) {

  if (kern=="gauss") {kern = "gaussian"}
  if (kern=="epan") {kern = "epanechnikov"}
  if (kern=="rect") {kern = "rectangular"}
  if (kern=="tria") {kern = "triangular"}
  if (kern=="bisq") {kern = "biweight"}
  if (kern=="cosine"|kern=="cos") {kern = "optcosine"}

  if (kern=="gaussian")     {k2 = 1/(2*sqrt(pi))}
  if (kern=="epanechnikov") {k2 = .6}
  if (kern=="rectangular")  {k2 = .5}
  if (kern=="triangular")   {k2 = 2/3}
  if (kern=="biweight")     {k2 = 5/7}
  if (kern=="optcosine")    {k2 = (pi^2)/16}

  lat1  <- 2*pi*lat1/360
  long1 <- 2*pi*long1/360
  lat2 <- 2*pi*lat2/360
  long2 <- 2*pi*long2/360
  n1 = length(long1)
  n2 = length(long2)
  if (nsamp >0) {sampvar <- sample(seq(1:n1),nsamp) }
  if (nsamp==0) {sampvar <- seq(1:n1)}
  la1 <- lat1[sampvar]
  lo1 <- long1[sampvar]
  radj = ifelse(kilometer==FALSE,3958,6371)

  maked <- function(lo,la) {
    coslat  <- cos(la)
    sinlat  <- sin(la)
    dvect <- sinlat%o%sinlat + (coslat%o%coslat)*cos(-outer(lo,lo,"-"))
    dvect <- dvect[lower.tri(dvect)]
    dvect <- ifelse(dvect>1,1,dvect)
    dvect <- radj*acos(dvect)
    return(dvect) 
  }

  dvect <- maked(lo1,la1)
  ndist = 2*length(dvect) 

  dmax = ifelse(dmax>0,dmax,max(dvect))
  h = ifelse(h>0,h, (.9*(quantile(dvect,.75)-quantile(dvect,.25))/1.34)*(ndist^(-.20)))
  distance <- seq(from=dmin,to=dmax,length=dlength)
  fit1 <- density(dvect,from=dmin,to=dmax,n=dlength,bw=h,kernel=kern) 
  fit2 <- density(-dvect,from=dmin,to=dmax,n=dlength,bw=h,kernel=kern)
  dhat <- fit1$y + fit2$y

  se <- sqrt(dhat*k2/(ndist*h))
  lo <- dhat + qnorm(pval/2)*se
  hi <- dhat + qnorm(1-pval/2)*se

  n = length(la1)
  n1 = n
  n2 = length(lat2) 
  coslat  <- cos(lat2)
  sinlat  <- sin(lat2)
  dbase <- sinlat%o%sinlat + (coslat%o%coslat)*cos(-outer(long2,long2,"-"))
  dbase <- ifelse(dbase>1,1,dbase)
  dbase <- radj*acos(dbase)

  rp = FALSE
  if (n1>=n2) {
    cat("n1>=n2.  Sample with replacement from n2 observations","\n")
    rp = TRUE
  }
  dmat <- array(0,dim=c(nsim,dlength))
  for (i in seq(1:nsim)) {
    obs <- sample(seq(1:n2),n1,replace=rp)
    o <- as.matrix(expand.grid(obs,obs))
    o <- o[o[,1]<o[,2],]
    dvect <- dbase[o]
    fit1 <- density(dvect,from=dmin,to=dmax,n=dlength,bw=h,kernel=kern) 
    fit2 <- density(-dvect,from=dmin,to=dmax,n=dlength,bw=h,kernel=kern) 
    dmat[i,] <- fit1$y + fit2$y
  }
  a0 = pval/2
  a1 = 1-a0
  qfunc0 <- function(x) {quantile(x,a0)}
  qfunc1 <- function(x) {quantile(x,a1)}
  local.lo <- apply(dmat,2,qfunc0)
  local.hi <- apply(dmat,2,qfunc1)

  global.lo = NULL
  global.hi = NULL

  if (cglobal==TRUE) {
    dvect <- apply(dmat,2,sort)
    tmat <- dvect[1,]
    f <- function(y) {y<tmat}
    growlo = 0
    gnum = a0*nsim
    for (i in seq(1:gnum)) {
      tmat <- dvect[i,]
      num <- sum(rowSums(t(apply(dmat,1,f))>0)>0)
      growlo = ifelse((num>gnum)&(growlo==0),i-1,growlo)
      if (growlo>0) {break}
    }
    growlo = ifelse(growlo==0,1,growlo)
    cat("Lower value for global CI = entry number",growlo,"\n")
    global.lo <- dvect[growlo,]

    f <- function(y) {y>tmat}
    growhi = 0
    for (i in seq(nsim,a1*nsim,-1)) {
      tmat <- dvect[i,]
      num <- sum(rowSums(t(apply(dmat,1,f))>0)>0)
      growhi = ifelse((num>gnum)&(growhi==0),i+1,growhi)
      if (growhi>0) {break}
    }
    growhi = ifelse(growhi==0,nsim,growhi)
    cat("Upper value for global CI = entry number",growhi,"\n")
    global.hi <- dvect[growhi,]
  }

  if (noplot==FALSE) {
    if (cglobal==TRUE) {
      lo <- global.lo
      hi <- global.hi
      yvect <- c(min(lo), max(hi))
      header <- "Global Confidence Intervals"
    }
    if (cglobal==FALSE) {
      lo <- local.lo
      hi <- local.hi
      yvect <- c(min(lo), max(hi))
      header <- "Local Confidence Intervals"
    }
    plot(distance,lo,xlab="Distance",ylab="Density",type="l",lty="dashed",ylim=yvect,main=header) 
    lines(distance,hi,lty="dashed")
  }

  out <- list(distance,dhat,h,local.lo,local.hi,global.lo,global.hi)
  names(out) <- c("distance","dhat","h","local.lo","local.hi","global.lo","global.hi")
  return(out)
}



