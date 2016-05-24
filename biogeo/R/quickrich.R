quickrich <-
function(dat,world,ID='ID',Species='Species',x='x',y='y',countries = "",others='',res,msk,ext){
  rs<-res
  dat1<-quickclean(world,dat,ID=ID,Species=Species,x=x,y=y,countries=countries,others=others,res=rs,msk=msk,ext=ext)
  # choose  a raster
  resn<-rs/60
  ra<-raster(xmn=-180, xmx=180, ymn=-90, ymx=90,res=resn,vals=NA)#
  ra[msk]<-msk
  
  if (class(ext)== "character"){
    fact<-0.1
    ex<-getextent(dat1$x,dat1$y, "p")
    mnx<-ex$xlm[1]
    mxx<-ex$xlm[2]
    xrng<-(mxx-mnx)*fact
    mnx2<-ifelse(mnx<0,mnx-xrng,mnx+xrng)
    mxx2<-ifelse(mxx<0,mxx-xrng,mxx+xrng)
    mny<-ex$ylm[1]
    mxy<-ex$ylm[2]
    yrng<-(mxy-mny)*fact
    mny2<-ifelse(mny<0,mny-yrng,mny+yrng)
    mxy2<-ifelse(mxy<0,mxy-yrng,mxy+yrng)
    ext<-c(mnx2,mxx2,mny2,mxy2)
  } else{
    ex<-getextent(dat1$x,dat1$y,ext=ext)
    ext<-c(ex$xlm,ex$ylm)
  }
  rst<-crop(ra,ext)
  rm(ra)
  
  rich <- richnessmap(dat1, rst, option="richness")
  return(rich)
}
