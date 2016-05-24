frameAutoSize <-
function(width,height,fdc){ 
  if(is.null(width)==TRUE&&is.null(height)==TRUE){
    width <- 10
    height<-23*width*(bbox(fdc)[2]-bbox(fdc)[4])/(bbox(fdc)[1]-bbox(fdc)[3])/20
  }else{
    if(is.null(width)==FALSE&&is.null(height)==TRUE){
      height<-23*width*(bbox(fdc)[2]-bbox(fdc)[4])/(bbox(fdc)[1]-bbox(fdc)[3])/20
    }else{
      if(is.null(width)==TRUE&&is.null(height)==FALSE){
        width<-17*height*(bbox(fdc)[1]-bbox(fdc)[3])/(bbox(fdc)[2]-bbox(fdc)[4])/20
      }
    }
  }
  dimFrame <- c(width,height)
  return(dimFrame)
}
