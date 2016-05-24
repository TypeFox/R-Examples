stplotGoogleMaps<-function(SPT,
                           zcol=1,
                           stfilename='spacetime.htm',
                           filename='file', # prefix for time data
                           w='100%',
                           h='49.5%',
                           openMap=FALSE,
                           colPalette=NULL,
                           do.bubble=FALSE,
                           at=NULL,
                           bubble= list(max.radius=10000,
                                        key.entries = if(do.bubble){quantile(SPT@data[,zcol],(1:5)/5, na.rm=T)}else{0},
                                         do.sqrt = TRUE),
                           iconMarker="",
                           ...) {
  
  


if(do.bubble){
  max.radius=bubble$max.radius
  key.entries=bubble$key.entries
  do.sqrt=bubble$do.sqrt
  SPT=SPT[!is.na(SPT@data[,zcol]),]
  z=SPT@data[,zcol]

  kkk=c(min(z),key.entries)
  kkk=sapply(2:length(kkk), function(i) mean( c(kkk[i],kkk[i-1]) ) )
  
  ke<-abs(kkk) + mean(abs(kkk))    # no zeros as max for radius vecor
  # creating a vector for subgroups
  if(do.sqrt){
    scale.level<- sqrt(ke/(max(ke)) ) }else{scale.level<-ke/(max(ke))}
  radius.level<-max.radius*scale.level
  # list of radiuses for createSphereCircle
  breakss<-factor(c(min(z),key.entries))
  break_unique<-as.numeric(levels(breakss))
  break_unique[length(break_unique)]<-max(z)
  if(length(unique(z))==length(key.entries)){ zz=factor(z,labels=radius.level)
                                              radius.vector<-floor(as.numeric(as.vector(zz))) 
  }else{ 
    zz=factor(cut(z,break_unique,include.lowest=TRUE ),labels=radius.level)
    radius.vector<-floor(as.numeric(as.vector((zz))))
  }
  
  
  if(is.null(colPalette) & min(key.entries)<0){
    colPalette=rep("#99000D",length(key.entries))
    colPalette[which(key.entries<0)]="#084594"
  }
  
  if(!is.null(colPalette)){
    rgbc<-col2rgb(colPalette)
    colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}
  
  if(length(colPalette)==1){
    colPalette=rep(colPalette,length(key.entries))
    rgbc<-col2rgb(colPalette)
    colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}
  
  if(length(colPalette)==2 & min(key.entries)<0){
    cop=rep(colPalette[2],length(key.entries))
    cop[which(key.entries<0)]=colPalette[1]
    rgbc<-col2rgb(cop)
    colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}
  
  colsb<-PolyCol(factor(zz,labels=key.entries),colPalette)$cols
} # end of do.bubble

if(!is.null(colPalette)){
  rgbc<-col2rgb(colPalette)
  colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}

attribute=SPT@data[,zcol]
cols=PolyCol(attribute,colPalette,at=at )$cols

if(class(SPT@sp)[1]=="SpatialPoints"){
  
  if(length(iconMarker)<length(attribute) ){
    iconMarker=iconlabels(attribute,cols,at,height=10,icon=TRUE,scale=0.6) }else{
      iconMarker=iconMarker[1:attribute] }
}



for(i in 1:length(SPT@data)) {
  if( identical(attribute,SPT@data[,i])){
    attributeName<-names(SPT@data)[i]  }
}

lti=0
wd=getwd()
dd=gsub(".h","_h",stfilename)

dir.create(dd)
setwd(dd)


if(class(SPT)[1]=="STFDF"){

  filename=paste(filename,"_",as.Date(SPT@endTime),sep="")
if(!do.bubble){
     for(i in 1: dim(SPT)[2]){
            m=plotGoogleMaps(SP=SPT[,i],filename=paste(filename[i],'.htm',sep=''),
                             zcol=zcol,
                             layerName = paste(attributeName,index(SPT@time[i]),sep=" "),
                   colPalette=cols[(1:length(SPT[,i]@data[,1])) +lti],openMap=openMap, ...) 
                           lti=lti+length(SPT[,i]@data[,1]) }
}else{
  for(i in 1: dim(SPT)[2]){
    tmp=bubbleSP(SPDF=SPT[,i],
                 zcol=zcol,
                 max.radius=bubble$max.radius,
                 radius.vector=radius.vector[(1:length(SPT[,i]@data[,1])) +lti])
    m=plotGoogleMaps(SP=tmp,
                     filename=paste(filename[i],'.htm',sep=''),
                     zcol=zcol,
                     layerName = paste(attributeName,index(SPT@time[i]),sep=""),
                     colPalette=colsb[(1:length(SPT[,i]@data[,1])) +lti] ,
                     openMap=openMap, ...)  
                    lti=lti+length(SPT[,i]@data[,1]) }
}

}else{
  SP=as(SPT, "Spatial")
  time_uni=unique(SP$time)
  filename=paste(filename,"_",as.Date(time_uni),sep="")
if(!do.bubble){  
          for(i in 1: length(filename)) {
              m=plotGoogleMaps(SP=SP[SP$time==time_uni[i],],zcol=zcol,
                   filename=paste(filename[i],'.htm',sep=''),
                   layerName = paste(attributeName,time_uni[i], sep=""),
                   colPalette=cols[(1:length(SP[SP$time==time_uni[i],]@data[,1])) +lti] ,openMap=openMap, ...) 
                       lti=lti+length(SP[SP$time==time_uni[i],]@data[,1]) }
   }else{

  for(i in 1: length(filename)) {
    tmp=bubbleSP(SPDF=SP[SP$time==time_uni[i],],
                 zcol=zcol,
                 max.radius=bubble$max.radius,
                 radius.vector=radius.vector[(1:length(SP@data[,1])) +lti])
    m=plotGoogleMaps(SP=tmp,
                     zcol=zcol,
                     filename=paste(filename[i],'.htm',sep=''),
                     layerName = paste(attributeName,time_uni[i], sep=""),
                     colPalette=colsb[(1:length(tmp@data[,1])) +lti],openMap=openMap, ...) 
                                    lti=lti+length(tmp@data[,1]) }
}
  
}


ifr=sapply(1:length(filename), function(i){
  paste('<iframe src="',getwd(),'/', filename[i],'.htm', '" width="',w,'" height="',h,'" align="left"  title=',i,' >',i,'</iframe>',sep='')
})

ifr=paste(ifr,collapse="")
ifr=paste(' <html> <head>  </head> <body>',ifr,' </body>  </html>',sep="",collapse=" ")

setwd(wd)
write(ifr, stfilename)
browseURL(stfilename)

}