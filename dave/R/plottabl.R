plottabl <-
function(veg,rorder=NULL,sorder=NULL,grr=NULL,grs=NULL,y=0.5) {
# ========================================================================
# plotting vegetation tables based on function image()  vers. 8.4.2012
# rorder and sorder are the orders, typically taken from o.hclust$order
# grr and gss are orders of group labels (factors), resulting from cuttree()
# y is for transformation of gray values
  sp.names<- names(veg)
  rel.names<- rownames(veg)
  sp.names<- strtrim(sp.names, 18)
  nrel <- length(veg[,1])
  nspec <- length(veg[1,])

# default handling
  l.rorder<- is.null(rorder)
  if(l.rorder == TRUE) {
    rorder<- rep(1:nrel,1)
    rorder<- order(rorder)
    grr<- rep(1,nrel) 
  }
  l.sorder<- is.null(sorder)
  if(l.sorder == TRUE) {
    sorder<- rep(1:nspec,1)
    sorder<- order(sorder)
    grs<- rep(1,nspec) 
  }

  l.grs<- is.null(grs)
  if(l.grs == TRUE) grs<- rep(1,nspec)
     
  l.grr<- is.null(grr)
  if(l.grr == TRUE) grr<- rep(1,nrel)

# reverse species order
#  sorder<- order(sorder,decreasing=TRUE)

# pixel size
  largedim<- max(nrel,nspec)
  psize<- 1/(largedim-1)
#  cat("psize", psize,"\n")
  hpsize<-psize/2

# transforming veg for plotting
  veg<- veg^y
  vrange<- range(-veg)

# setting up multiple of a wsz x wsz matrix
  wsz<- max(c(nrel,nspec))
  mrow<- ceiling(nspec/wsz)
  mcol<- ceiling(nrel/wsz)
  
  pmatrix<- matrix(rep(0,wsz*wsz*mrow*mcol),ncol=wsz*mrow)
  rn<- seq(1,nrel,1)
  sn<- seq(1,nspec,1)
  ind <- as.matrix(expand.grid(rn,sn))
  pmatrix[ind]<- veg[ind]
# order within pmatrix
  o.py<- seq(1,mrow*wsz,1)
  o.px<- seq(1,mcol*wsz,1)
  o.px[1:nrel]<- rorder
  o.py[nspec:1]<- sorder
  par(mfrow=c(1,1),mar=c(0,0,0,0),omi=c(0,0,0,0))
# plot matrix
# -----------
  plot(c(-0.10,1.05),c(-0.10,1.05),asp=1,type="n",axes=FALSE)
  image(-pmatrix[o.px,o.py],zlim=vrange,col=gray((0:32)/32),add=TRUE,useRaster = TRUE)
   yt<- c(0,0.3*(nspec-1)*psize)
   xt<- c(-0.04,-0.04)
#   lines(xt,yt,lwd=0.5)
   yt<- c(0.7*(nspec-1)*psize,1.0*(nspec-1)*psize)
#   lines(xt,yt,lwd=0.5)
   text(-0.04,0.5*(nspec-1)*psize,"Species",srt=90,cex=0.6)
   xr<- c(0,0.3*(nrel-1)*psize)
   yr<- c(1.05*(nspec-1)*psize,1.05*(nspec-1)*psize)
#   lines(xr,yr,lwd=0.5)
   xr<- c(0.7*(nrel-1)*psize,1.0*(nrel-1)*psize)
#   lines(xr,yr,lwd=0.5)
   text(0.5*(nrel-1)*psize,1.05*(nspec-1)*psize,"Releves",cex=0.6)
# new lines
  lines(c(0-hpsize,0-hpsize),c(0-hpsize,(nspec*psize)-hpsize+0.001),lwd=0.5)                              # left
  lines(c(0-hpsize,(nrel*psize)-hpsize),c(0-hpsize,0-hpsize),lwd=0.5)                                     # below
  lines(c(0-hpsize,(nrel*psize)-hpsize),c((nspec*psize)-hpsize+0.001,(nspec*psize)-hpsize+0.001),lwd=0.5) # top
  lines(c((nrel*psize)-hpsize,(nrel*psize)-hpsize),c(0-hpsize,(nspec*psize)-hpsize+0.001),lwd=0.5)        # right
# lines separating the species groups
  iposy<- 0
  o.set<- setgroupsize(grs[o.py])
  for(k in 1:(o.set$ngroups)) {
      iposy<- iposy+o.set$groupcounts[k]*psize
      hiposy<- iposy-(o.set$groupcounts[k]*psize*0.5)
      lines(c(0-hpsize,(nrel*psize)-hpsize),c(iposy-hpsize,iposy-hpsize),col=gray(0.5),lwd=0.1)
#      text(nrel*psize,hiposy-hpsize,o.set$grouplabs[k],cex=0.8,pos=4)
  }
# lines separating the releve groups
  iposx<- 0
  o.set<- setgroupsize(grr[o.px])
  for(k in 1:(o.set$ngroups)) {
      iposx<- iposx+o.set$groupcounts[k]*psize
      hiposx<- iposx-(o.set$groupcounts[k]*psize*0.5)
      lines(c(iposx-hpsize,iposx-hpsize),c(0-hpsize,(nspec*psize)-hpsize+0.001),col=gray(0.5),lwd=0.1)
      text(hiposx-hpsize,0-hpsize,o.set$grouplabs[k],cex=0.8,pos=1)
  }
 }
