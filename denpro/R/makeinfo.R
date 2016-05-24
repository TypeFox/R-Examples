makeinfo<-function(left,right,mean,low,upp)
{
lehdet<-findleafs(left,right)

d<-dim(low)[2]
nodenum<-length(lehdet)         #length(left)

value<-matrix(0,nodenum,1)
infolow<-matrix(0,nodenum,d)
infoupp<-matrix(0,nodenum,d)
nodefinder<-matrix(0,nodenum,1)
infopointer<-matrix(0,nodenum,1)

runner<-1
leafnum<-0
while (runner<=nodenum){
  if ((!is.na(lehdet[runner])) && (lehdet[runner]==1) && (mean[runner]>0)){  
      # we are in leaf where the value is positive
      leafnum<-leafnum+1
      value[leafnum]<-mean[runner]
      nodefinder[leafnum]<-runner
      infolow[leafnum,]<-low[runner,]
      infoupp[leafnum,]<-upp[runner,]

      infopointer[runner]<-leafnum
  }
  runner<-runner+1
}
value<-value[1:leafnum]
nodefinder<-nodefinder[1:leafnum]
infolow<-infolow[1:leafnum,]
infoupp<-infoupp[1:leafnum,]

return(list(value=value,low=infolow,upp=infoupp,nodefinder=nodefinder,
infopointer=infopointer,
terminalnum=leafnum))
}
