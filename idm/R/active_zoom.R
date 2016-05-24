active_zoom<-function(outmo,nframes=25,att=TRUE){
  nc=outmo$nchunk
  out=list()
  
  if(att==TRUE){
    lout=list()
    lout$xr=cbind(seq(from=outmo[[(1)]]$xr[1],to=outmo[[(1)]]$xr[1],length=nframes),
                  seq(from=outmo[[(1)]]$xr[2],to=outmo[[(1)]]$xr[2],length=nframes))*1.01
    
    lout$yr=cbind(seq(from=outmo[[(1)]]$yr[1],to=outmo[[(1)]]$yr[1],length=nframes),
                  seq(from=outmo[[(1)]]$yr[2],to=outmo[[(1)]]$yr[2],length=nframes))*1.01
    
    out[[1]]=lout
    
    for(i in 2:nc){
      lout=list()    
      lout$xr=cbind(seq(from=outmo[[(i-1)]]$xr[1],to=outmo[[(i)]]$xr[1],length=nframes),
                seq(from=outmo[[(i-1)]]$xr[2],to=outmo[[(i)]]$xr[2],length=nframes))*1.01
      
      lout$yr=cbind(seq(from=outmo[[(i-1)]]$yr[1],to=outmo[[(i)]]$yr[1],length=nframes),
                seq(from=outmo[[(i-1)]]$yr[2],to=outmo[[(i)]]$yr[2],length=nframes))*1.01
      
      out[[i]]=lout
    }
  }else{
    lout=list()
    lout$xr=cbind(seq(from=outmo[[(1)]]$uxr[1],to=outmo[[(1)]]$uxr[1],length=nframes),
                  seq(from=outmo[[(1)]]$uxr[2],to=outmo[[(1)]]$uxr[2],length=nframes))*1.01
    
    lout$yr=cbind(seq(from=outmo[[(1)]]$uyr[1],to=outmo[[(1)]]$uyr[1],length=nframes),
                  seq(from=outmo[[(1)]]$uyr[2],to=outmo[[(1)]]$uyr[2],length=nframes))*1.01
    
    out[[1]]=lout
    for(i in 2:nc){
      lout=list()
      lout$xr=cbind(seq(from=outmo[[(i-1)]]$uxr[1],to=outmo[[(i)]]$uxr[1],length=nframes),
                seq(from=outmo[[(i-1)]]$uxr[2],to=outmo[[(i)]]$uxr[2],length=nframes))*1.01
      
      lout$yr=cbind(seq(from=outmo[[(i-1)]]$uyr[1],to=outmo[[(i)]]$uyr[1],length=nframes),
                seq(from=outmo[[(i-1)]]$uyr[2],to=outmo[[(i)]]$uyr[2],length=nframes))*1.01
      
      out[[i]]=lout
    } 
  }
  
  
return(out)
}