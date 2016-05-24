plotCI.yhat<-function(sampStat, upperCI, lowerCI, pid=1:ncol(sampStat),nr=2,nc=2){

par(mfrow=c(nr,nc))
for (j in 1:length(pid)){
  i<-pid[j]
  k<-nrow(sampStat)
  plotCI(x=sampStat[1:k,i],uiw=upperCI[1:k,i]-sampStat[1:k,i],
       liw=sampStat[1:k,i]-lowerCI[1:k,i],ylab="", xlab="",
       col="black",xaxt="n",xlim=c(0,(k+1)), ylim=c(min(na.omit(lowerCI[1:k,i])),max(na.omit(upperCI[1:k,i]))),pch=21,pt.bg=par("bg"))
  axis(side=1,at=1:k,labels=rownames(sampStat),las=2)
  title(main=colnames(sampStat)[[i]])
}
return()
}
