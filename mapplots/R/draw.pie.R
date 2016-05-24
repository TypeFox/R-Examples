draw.pie <-
function(x,y,z,radius,scale=T,labels=NA,silent=TRUE,...){
  nx <- length(x)
  nz <- dim(z)[2]
  if(length(y)!=nx) stop('x and y should be vectors of the same length')
  if(length(dim(z))!=2) stop('z should be a 2-dimensional array')
  if(dim(z)[1]!=nx) stop('the number of rows in of z should match as the length of x and y')
  if(sum(z,na.rm=T)==0) stop('z has no data')
  maxsumz <- max(rowSums(z),na.rm=T)
  pm <- setProgressMsg(1,nx)
  for(i in 1:nx){
    xi <- x[i]
    yi <- y[i]
    zi <- z[i,]
    zi <- ifelse(is.na(zi),0,zi)
    if(length(radius)>1) radiusi <- radius[i] else radiusi = radius
    if(scale & length(radius)==1) radiusi <- radius*sqrt(sum(zi,na.rm=T))/sqrt(maxsumz)
    if(sum(zi)>0) add.pie(zi,xi,yi,labels,radius=radiusi,...)
    if(!silent) pm <- progressMsg(pm,i)
  }
}

