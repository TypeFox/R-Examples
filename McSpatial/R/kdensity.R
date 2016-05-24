kdensity <- function(longitude,latitude,kilometer=FALSE,noplot=FALSE,dmin=0,dmax=0,dlength=512,h=0,kern="gaussian",nsamp=0,confint=TRUE,pval=.05) {

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

  n = length(longitude)
  if (nsamp>0) {sampvar <- sample(seq(1:n),nsamp) }
  if (nsamp==0) {sampvar <- seq(1:n)}
  lat1  <- 2*pi*latitude[sampvar]/360
  long1 <- 2*pi*longitude[sampvar]/360
  sinlat  <- sin(lat1)
  coslat  <- cos(lat1)
  dvect <- sinlat%o%sinlat + (coslat%o%coslat)*cos(-outer(long1,long1,"-"))
  dvect <- dvect[lower.tri(dvect)]
  dvect <- ifelse(dvect>1,1,dvect)
  radj = ifelse(kilometer==FALSE,3958,6371)
  dvect <- radj*acos(dvect)
  n = 2*length(dvect)

  if (dmin==0&dmax==0) {dmax = max(dvect)} 
  if (h==0) {
    h = (.9*(quantile(dvect,.75)-quantile(dvect,.25))/1.34)*(n^(-.20))
    dfit1 <- density(dvect,from=dmin,to=dmax,n=dlength,kernel=kern,bw=h)
    dfit2 <- density(-dvect,from=dmin,to=dmax,n=dlength,kernel=kern,bw=h)
  }
  if (h>0)  {
    dfit1 <- density(dvect,from=dmin,to=dmax,n=dlength,kernel=kern,bw=h)
    dfit2 <- density(-dvect,from=dmin,to=dmax,n=dlength,kernel=kern,bw=h)
  }
  distance <- dfit1$x
  dhat <- dfit1$y + dfit2$y

  se <- sqrt(dhat*k2/(n*h))
  lo <- dhat + qnorm(pval/2)*se
  hi <- dhat + qnorm(1-pval/2)*se

  if (noplot==FALSE) { 
    ymin = ifelse(confint==TRUE,min(lo),min(dhat))
    ymax = ifelse(confint==TRUE,max(hi),max(dhat))
    plot(distance,dhat,xlab="Distance",ylab="Density",type="l",ylim=c(ymin,ymax)) 
    if (confint==TRUE) {
      lines(distance,lo,lty="dashed")
      lines(distance,hi,lty="dashed")
    }
  }
  out <- list(distance,dhat,dvect,h,se)
  names(out) <- c("distance","dhat","dvect","h","se")
  return(out)
}



