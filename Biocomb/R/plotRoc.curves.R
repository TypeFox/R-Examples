plotRoc.curves <- function(dattable,file.name=NULL,colours=NULL,ltys=NULL,add.legend=F,
                           curve.names=NULL,include.auc=F,xaxis="",yaxis="",line.width=2,headline="",ispercent=F)
  {
  labs <- dattable[,ncol(dattable)]
  scores<-dattable[,-ncol(dattable)]

  if (is.null(ncol(scores))){scores <- data.frame(scores)}
  if (is.null(colours)){colours <- 1:ncol(scores)}
  if (is.null(ltys)){
  ltys <- rep(1,ncol(scores))
  if(ncol(scores)>8){
    N<-ncol(scores)%/%8
    for(i in 1:N){
      ltys[(8*i+1):length(ltys)]<-ltys[(8*i+1):length(ltys)]+1

    }
  }}#####if ncol(scores)>8 ltys should be different.ncol(scores)
  aucvals <- rep(0,ncol(scores))
  pred <- prediction(scores[,1],labs)
  aucv <- performance(pred,"tpr", "fpr",measure="auc")
  aucval <- attr(aucv,"y.values")[[1]]
  if (aucval<0.5){
    aucval <- 1-aucval
    pred <- prediction(-scores[,1],labs)
  }
  aucvals[1] <- round(1000*aucval)/1000
  perf <- performance(pred,"tpr", "fpr")

 if (!is.null(file.name)){pdf(file=file.name)}

    if (xaxis=="" & yaxis==""){
      xaxis = "False positive rate"
      yaxis = "True positive rate"
      if (ispercent){
        xaxis <- paste(xaxis,"(%)",sep=" ")
        yaxis <- paste(yaxis,"(%)",sep=" ")
      }
    }
    if (ispercent){
      attr(perf,"x.values")[[1]] <- attr(perf,"x.values")[[1]]*100
      attr(perf,"y.values")[[1]] <- attr(perf,"y.values")[[1]]*100
    }

    plot(perf,lwd=line.width,col=colours[1],lty=ltys[1],xlab=xaxis,ylab=yaxis,main=headline)

    if (ncol(scores)>1){
      for (i in 2:ncol(scores)){
        pred <- prediction(scores[,i],labs)
        aucv <- performance(pred,"tpr", "fpr",measure="auc")
        aucval <- attr(aucv,"y.values")[[1]]
        if (aucval<0.5){
          aucval <- 1-aucval
          pred <- prediction(-scores[,i],labs)
        }
        aucvals[i] <- round(1000*aucval)/1000
        perf <- performance(pred,"tpr", "fpr")

        if (ispercent){
          attr(perf,"x.values")[[1]] <- attr(perf,"x.values")[[1]]*100
          attr(perf,"y.values")[[1]] <- attr(perf,"y.values")[[1]]*100
        }
        plot(perf,lwd=line.width,col=colours[i],lty=ltys[i],add=T)
      }
    }


    if (add.legend){

      if (is.null(curve.names)){curve.names=names(scores)}
      leg.text <- curve.names
      if (include.auc){
        for (i in 1:ncol(scores)){
          leg.text[i] <- paste(curve.names[i],",  AUC=",aucvals[i],sep="")
        }
      }
      legend("bottomright",leg.text,lwd=line.width,lty=ltys,col=colours,cex =(0.3+1.4/(ncol(scores)%/%15+2)))
    }

  if (!is.null(file.name)){dev.off()}

}
