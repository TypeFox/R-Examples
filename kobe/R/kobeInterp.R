
kobeInterp=function(x,levels=seq(0.0,1.0,0.05),
                  col   =c(colorRampPalette(c("red4","red"))(12),colorRampPalette(c("yellowgreen","darkgreen"))(8)),
                  nIterp=101){

  x=x[!is.na(x[,1]) & !is.na(x[,2]) & !is.na(x[,3]),]
  
  return(NULL)
  ##### smooth
  ## uses akima
#   t.<-interp(x[,1],x[,2],x[,3],
#                     xo=seq(min(x[,1]),   max(x[,1]), length=nIterp),
#                     yo=seq(min(x[,2]),   max(x[,2]), length=nIterp),
#                     duplicate="mean")
#   
#   
#   res=cbind(expand.grid(x=t.$x,y=t.$y),z=cut(t.$z,levels,include.lowest=T),w=c(t.$z))
#   res$col=col[as.numeric(res$z)]
#   
#  res
}

