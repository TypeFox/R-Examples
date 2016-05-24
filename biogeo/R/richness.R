richness <-
function(dat,res=10,option="richness",buf=5,ext=""){
    
  if(res<=0){stop("res must be >0 minutes")}
    f1<-ifelse(option=="richness",1,0)
    f2<-ifelse(option=="records",1,0)
    if ((f1+f2)<1){stop("invalid option")}
  resn<-res
  ra<-raster(xmn=-180, xmx=180, ymn=-90, ymx=90,res=resn/60)#
  nrw<-nrow(ra)
  ncl<-ncol(ra)
  cn<-colFromX(ra, dat$x)
  left1<-min(cn)-buf
  left<-ifelse(left1<1,1,left1)
  rt1<-max(cn)+buf
  rt<-ifelse(rt1>ncl,ncl,rt1)
  mnx<-xFromCol(ra,left)
  mxx<-xFromCol(ra,rt)
  rn<-rowFromY(ra, dat$y)
  bot1<-max(rn)+buf
  top1<-min(rn)-buf
  bot<-ifelse(bot1>nrw,nrw,bot1)
  top<-ifelse(top1<1,1,top1)
  mny<-yFromRow(ra,bot)
  mxy<-yFromRow(ra,top)
  r<-raster(xmn=mnx, xmx=mxx, ymn= mny, ymx=mxy,res=resn/60,vals=0)
  cx<-class(ext)
  if (cx[1]=="character"){r<-r}
  if (cx[1]=="numeric"){
    ext1<-extent(ext)
    r<-crop(r,ext1)}
  if (cx[1]=="Extent"){r<-crop(r,ext)}
  rst<-r  
  rst0<-r
    v<-values(rst0)
    xy<-data.frame(dat$x,dat$y)
    vals <- extract(rst0,xy) # extract values from the raster
    f<-which(!is.na(vals)) # first find those that are in the sea
    ce1<-cellFromXY(rst, xy[f,])
    dat2<-data.frame(Species=dat$Species[f],cell=ce1,n=1)
    if (option=="richness"){
      # species richness with only one record per species
      u<-unique(dat2$cell)
      d<-rep(0,length(u))
      #b<-aggregate(Species~cell,length,data=dat2) # count records per species
      for (i in 1:length(u)){
        bi<-dat2[dat2$cell==u[i],]
        d[i]<-length(unique(bi$Species))
      }
      v[u]<-d  # set the values
      rich <-setValues(rst,v)
    } 
    if (option=="records"){
      nrec<-aggregate(n~cell,sum,data=dat2) # no of records per cell
      v[nrec$cell]<-nrec$n  # set the values
      rich <-setValues(rst,v)
    } 
    return(rich)
  }
