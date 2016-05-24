read.rate<-function(data)
 { 
    no_col <- max(count.fields(data, sep=""))
    k<-1:((no_col-17)/2)
    read.fwf(data, skip=3,fill=T,widths=c(5,7,1,6,1,3,3,3,5,5,9,3,5,4,5,6,4, 
             rep(c(3,5),times=(no_col-17)/2)),col.names=c("IMOcode","long","EW",
             "lat","NS","day","month","year","start","stop","sollong","fovRA","fovDEC",
             "Teff","F","lmg","SPO",paste(c("Shw","N"),rep(k,each=2),sep="")))
 }
