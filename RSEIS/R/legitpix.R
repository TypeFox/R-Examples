legitpix<-function(sel, zloc, zenclick)
  {
    ############   sel is the selected index of traces
    ############   zloc are the clicks on the screen
    ############   zenclick is the number of clicks
    
    legpick = length(sel)-floor(length(sel)*zloc$y[1:(zenclick-1)])
    ileg = which(legpick>=1 & legpick<=length(sel))
    
    ypick  = legpick[ ileg ]
    ppick  = zloc$x[ ileg ]
    
    return( list(ypick=ypick, ppick=ppick, ileg=ileg) )
    
  }

