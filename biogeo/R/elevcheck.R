elevcheck <-
function(dat, dem, elevc="elevation",diff=50){
    
    fieldsmissing(dat, fields=c("ID","Species","x","y","x_original","y_original","Correction","Modified","Exclude","Reason"))
    n<-nrow(dat)
    x <- coord2numeric(dat$x)
    y <- coord2numeric(dat$y)
    xy <- cbind(x,y)
    dat$demElevation <- extract(dem, xy)
    f <- which(!is.na(dat[,elevc]))
    altdat <- dat[f,]
    x <- coord2numeric(altdat$x)
    y <- coord2numeric(altdat$y)
    xy <- cbind(x,y)
    cid <- cellFromXY(dem, xy)
    dups <- duplicated(cid) * 1
    reslm <- altdat$demElevation-dat[,elevc]
    elevMismatch<-rep(0,n)
    evm<- (abs(reslm)>=diff)*1
    elevMismatch[f]<-evm
    demElevation<-dat[,elevc]
    demElevation[f]<-altdat$demElevation
    z<-data.frame(elevMismatch,demElevation)
    
    return(z)
  }
