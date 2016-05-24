#' plot1= joint probability plot #
#' @export
#' @param  obj is output object from the binteg function - only input required
#' @param  z is a table of structure function variables
#' @param  px is the annual exceedance probability associated  with rows 
#' @param  py is the annual exceedance probability associated with columns
#' @param  prm is parameter scalar/vector for model being used
#' @param  out is the annual exceedance probability of output
#' @param  prob is how the probabilities are interpreted - assumed AEP 
#' @param  model is density model being used - assumed logistic
#' @param ...  is for other graphical parameters
#' alternative specification instead of obj - you can call this plot without using binteg
#' can specify (px,py,z,model,prob,prm)  OR as obj
#' the former is if you want to check data before integration
#' the latter is for output from binteg function
plot1=function(obj,prob){
  
  model=obj$model
  if(model=="log")    G=Glog
  if(model=="neglog") G=Gneglog   # not tested yet, not working
  #if(model=="beta")   G=Gbeta     # not tested yet, not working
  f=function(x) -1/log(1-10^(-x)) # inverse transform for frechet scale (and log10 scale)
  
  # get input from object
  if(obj$prob=="AEP"){
    px.aep=obj$px;py.aep=obj$py;  # Annual exceedance probabilities
  }else if(obj$prob=="ARI"){
    px.aep=1-exp(-1/obj$px); py.aep=1-exp(-1/obj$py); # Annual exceedance probabilities
  }
  # convert probability scale
  px.anep=1-px.aep;             py.anep=1-py.aep;             # Annual non-exceedance probabilities
  px.ari=-1/log(px.anep);       py.ari=-1/log(py.anep);       # Annual Average Recurrence Intervals
  if (px.anep[1]==0){px.ari[1]=1/365}
  if (py.anep[1]==0){py.ari[1]=1/365} 
  px.dri=px.ari*365;            py.dri=py.ari*365;            # Daily Average Recurrence Intervals, # Note that lowest bound (marginal) is given a daily recurrence here
  lpx.dri=log10(px.dri);        lpy.dri=log10(py.dri);        # Log10 scale for easier handling
 
  # generate dimensions of pdf grid
  tmp_dri <- min(max(lpx.dri),max(lpy.dri))
  lpx = seq(0,tmp_dri,0.1)
  lpy = seq(0,tmp_dri,0.1)
  nx=length(lpx);ny=length(lpy)
  # generate cdf on grid
  G.grid=outer(f(lpy),f(lpx),G,obj$prm) 
  # generate pdf via numerical derivative g=dG/dxdy
  g.grid=(G.grid[2:ny,2:nx] - G.grid[1:(ny-1),2:nx] - G.grid[2:ny,1:(nx-1)] + G.grid[1:(ny-1),1:(nx-1)]) 
  g.grid=rbind(0,cbind(0,g.grid)) # add row/col of 0 to make conformable 

  # deal with the vertical lines
   tmp_fz <- 1E-6;z_z <- obj$oz;k=1
   for (i in 1:nrow(obj$oz)) { 
   for (j in 1:ncol(obj$oz)) {
    z_z[i,j]<- obj$oz[i,j]+k* tmp_fz
    k=k+1}}
    z <- z_z

 if (prob=="ARI"){
  # plot the pdf
  plot(x=NULL,y=NULL, type="n", xlim=c(0,max(lpx)),ylim=c(0,max(lpy)),axes=F, xlab="X variable (ARI)", ylab="Y variable (ARI)")
  uu = par("usr")
  rect(uu[1], uu[3], uu[2], uu[4],  border = "black",lwd=4)
  xlab=c(round(px.ari[which(px.ari<=1)],1),round(px.ari[which(px.ari>1)],1))
  ylab=c(round(py.ari[which(py.ari<=1)],1),round(py.ari[which(py.ari>1)],1))
  axis(1, at=lpx.dri, labels=xlab, tick=1, lty=1)
  axis(2, at=lpy.dri, labels=ylab, tick=1, lty=1)
  contour(lpx,lpy,g.grid, col="blue", levels=c(0.0001,0.00001,0.000001,0.0000001,0.00000001), add=TRUE, lwd=2,drawlabels = F)
  # plot the design variable contours  
  nz.pretty=10 # number of z thresholds to plot
  ninc1=100     # number of intervals along line
  thresh.pretty = (min(z)+(1:nz.pretty)*diff(range(z))/(nz.pretty+1)) # contour levels
  h.pretty=lapply(contourLines(x=lpx.dri,y=lpy.dri,z=z,levels=thresh.pretty),approx,n=ninc1) # interpolated (x,y) coords of contours
  # plot the contour lines  
  for(i in 1:length(h.pretty)) {
    lines(h.pretty[[i]],col="red",lwd=2)
    text(h.pretty[[i]]$x[nx/2],h.pretty[[i]]$y[ny/2],format(thresh.pretty[[i]],digits=2),col="red",adj=0)}
  }

   if (prob=="AEP"){
  # plot the pdf
  plot(x=NULL,y=NULL, type="n", xlim=c(0,max(lpx)),ylim=c(0,max(lpy)),axes=F, xlab="X variable (AEP)", ylab="Y variable (AEP)")
  uu = par("usr")
  rect(uu[1], uu[3], uu[2], uu[4],  border = "black",lwd=4)
  xlab=round(px.aep,4)
  ylab=round(py.aep,4)
  axis(1, at=lpx.dri, labels=xlab, tick=1, lty=1)
  axis(2, at=lpy.dri, labels=ylab, tick=1, lty=1)
  contour(lpx,lpy,g.grid, col="blue", levels=c(0.0001,0.00001,0.000001,0.0000001,0.00000001), add=TRUE, lwd=2,drawlabels = F)

  # plot the design variable contours  
  nz.pretty=10 # number of z thresholds to plot
  ninc1=100     # number of intervals along line
  thresh.pretty = (min(z)+(1:nz.pretty)*diff(range(z))/(nz.pretty+1)) # contour levels
  h.pretty=lapply(contourLines(x=lpx.dri,y=lpy.dri,z=z,levels=thresh.pretty),approx,n=ninc1) # interpolated (x,y) coords of contours
  # plot the contour lines  
  for(i in 1:length(h.pretty)) {
    lines(h.pretty[[i]],col="red",lwd=2)
    text(h.pretty[[i]]$x[nx/2],h.pretty[[i]]$y[ny/2],format(thresh.pretty[[i]],digits=2),col="red",adj=0)}
  }

}
 
