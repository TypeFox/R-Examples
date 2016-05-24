`extract.val` <-
function(latlon, layers){
m=length(layers)
res=as.data.frame(cbind(latlon$lat,latlon$lon))
names(res)=c("lat","long")

for (k in 1:m){
  alayer=readGDAL(layers[[k]])	
  #alayer=layers[[k]]
  #check format
  if(class(alayer)=="SpatialGridDataFrame") {
	alayer=as.image.SpatialGridDataFrame(alayer)
  }  
  if(k==1){
    idx = latlon2idx(latlon$lat,latlon$lon,alayer)
    res1=as.data.frame(cbind(idx$x,idx$y))
    names(res1)=c("y","x")
    res=cbind(res,res1)
  }
  n=length(idx$x)

  out=matrix(NA,ncol=1,nrow=n)
  for(i in 1:n){
    out[i]=alayer$z[idx$x[i], idx$y[i]]
  }
  res1=as.data.frame(out)
  names(res1)=names(layers[k])
  res=cbind(res,res1)
}
return(res)
}

