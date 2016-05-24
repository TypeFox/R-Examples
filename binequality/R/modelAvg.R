modelAvg <-
function(fitList,ID,nonCon=TRUE){
	fitList$distribution<-as.character(fitList$distribution)
  dists<-unique(fitList[,ID])
  outAIC<-c()
  outBIC<-c()
  outAICavg<-c()
  outBICavg<-c()
  outWeightsAIC<-c()
  outWeightsBIC<-c()
  outWeightsAIC<-rep(NA, nrow(fitList))
  outWeightsBIC<-rep(NA, nrow(fitList))
  for(i in dists){
    use.i<-which(fitList[,ID]==i)
    dat.i<-fitList[use.i,]
    if(nonCon==TRUE){
      con.i<-which(dat.i$didConverge==TRUE&is.na(dat.i$estMean)==FALSE)
    }else{
      con.i<-which(is.na(dat.i$est)==FALSE)
    }
    if(length(con.i)==0){
      bestMod<-'nothingConverged'
      avgAIC.i<-'nothingConverged'
      avgBIC.i<-'nothingConverged'
    }else{
      bestMod.aic<-which(dat.i$aic[con.i]==min(dat.i$aic[con.i],na.rm=TRUE))
      bestMod.bic<-which(dat.i$bic[con.i]==min(dat.i$bic[con.i],na.rm=TRUE))
      use.AIC.i<-which(fitList[,ID]==i&fitList$distribution==dat.i$distribution[con.i][bestMod.aic])[1]
      use.BIC.i<-which(fitList[,ID]==i&fitList$distribution==dat.i$distribution[con.i][bestMod.bic])[1]
      wAIC.i<-makeWeightsAIC(dat.i$aic[con.i])
      wBIC.i<-makeWeightsAIC(dat.i$bic[con.i])
      avgAIC.i<-mAvg(dat.i$estMean[con.i],wAIC.i)
      avgBIC.i<-mAvg(dat.i$estMean[con.i],wBIC.i)
      outAIC<-c(outAIC, use.AIC.i)
      outBIC<-c(outBIC, use.BIC.i)
      outWeightsAIC[use.i[con.i]]<-wAIC.i
      outWeightsBIC[use.i[con.i]]<-wBIC.i
      #outWeightsAIC<-c(outWeightsAIC,wAIC.i)
      #outWeightsBIC<-c(outWeightsBIC,wBIC.i)
    }#end if/else con.i
    outAICavg<-c(outAICavg,avgAIC.i)
    outBICavg<-c(outBICavg,avgBIC.i)
  }#end for i in dists
  aicOut<-matrix(NA, ncol=ncol(fitList), nrow=length(dists))
  aicOut<-as.data.frame(aicOut)
  use.aic<-which(outAICavg!='nothingConverged')
  if(length(use.aic)>0){
  	aicOut[use.aic,]<-fitList[outAIC,]
  }
  aicOut[,1]<-dists
  colnames(aicOut)<-colnames(fitList)
  bicOut<-matrix(NA, ncol=ncol(fitList), nrow=length(dists))
  bicOut<-as.data.frame(bicOut)
  use.bic<-which(outBICavg!='nothingConverged')
  if(length(use.bic)>0){
  	 bicOut[use.bic,]<- fitList[outBIC,]
  }
  bicOut[,1]<-dists
  colnames(bicOut)<-colnames(fitList)
  out<-list('aic'= aicOut,'bic'= bicOut,'aicAvg'=outAICavg,'bicAvg'=outBICavg,'waic'=outWeightsAIC,'wbic'=outWeightsBIC)
  return(out)
}
