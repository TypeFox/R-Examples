##' @export
plotpredictsmcure <-
function(object, type="S", xlab="Time",ylab="Predicted Survival Probability",model=c("ph","aft"), ...)
{
  pred <- object$prediction
  if(model=="ph"){ 
    pdsort <- pred[order(pred[,"Time"]),]
    if(length(object$newuncureprob)==1) plot(pdsort[,"Time"],pdsort[,1], type="S")
    else 
      matplot(pdsort[,"Time"],pdsort[,1:(ncol(pred)-1)],col=1,type="S",lty=1:(ncol(pred)-1),xlab=xlab,ylab=ylab)
    }
  if(model=="aft"){
    nplot=ncol(pred)/2
    pdsort <- pred[order(pred[,1+nplot]),c(1,1+nplot)]
    plot(pdsort[,2],pdsort[,1],xlab=xlab,ylab=ylab,col=1,type="S",ylim=c(0,1))
    if(nplot>1){
      for(i in 2:nplot){
        pdsort<- pred[order(pred[,i+nplot]),c(i,i+nplot)]
        lines(pdsort[,2],pdsort[,1],lty=i,type="S")
        }
      }
    }
  }

