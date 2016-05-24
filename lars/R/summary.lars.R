summary.lars <-
function(object,sigma2=NULL,...){
  heading=c(paste("LARS/",object$type,sep=""),paste("Call:",format(object$call)))
  
  df=object$df
  rss=object$RSS
  K=length(object$actions)
  stepno=c(0,seq(K))
  Cp=object$Cp  
  if(!missing(sigma2)){
    n=attr(Cp,"n")
    Cp=rss/sigma2 -n +2*df
  }
sumob=data.frame(Step=stepno,Df=df,Rss=rss,Cp=Cp,row.names=1)
  structure(sumob,heading=heading,class=c("anova","data.frame"))
}

