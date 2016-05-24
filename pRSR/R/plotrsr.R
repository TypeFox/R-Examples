plotrsr <-
function(n=20:50, q=c(90,95), ...){
 # n: sequence of series size
 # q: vector of quantile values
 # ... : argument for function lines()
 probs<-TPout$probs 
 quant<-probs*100
 nn<-TPout$n  # seq(6,70,2) : the original sequence of  n that was run for RSR
 if(any(!(q%in%quant)))stop("One of the quantiles is out of bound. Type 'TPout$probs*100'
     to see possible quantile values.")
   MV=TPout$MeanVar
  m<-length(quant)
  w<-quant%in%q
  ind<-(1:m)[w]
  LM=~I(1/nn)+I(1/nn^2)
  Qts=paste(c("Qt"),probs[ind]*100,sep="")
  Wt=paste(c("1/Var"),probs[ind]*100,sep="")
  Regressor=paste(as.character(LM),sep="")
  Regressor=sapply(Qts, function(a) paste(a,"~", Regressor[2], sep=""))
  AllMod=function(X){
        apply(cbind(Regressor, Wt),1,  function(Regressor){
        V=eval(parse(text=Regressor[2]),X, parent.frame())
        do.call("lm", list(as.formula(Regressor[1]),  weights=quote(V), data=quote(X)))})
    }
 Res=lapply(MV, AllMod)[[1]]
 pred.val<-sapply(Res, function(x,n) predict(x, data.frame(nn=n)),n)
 y.min=min(pred.val)
 y.max=max(pred.val)
 plot(c(min(n),max(n)),c(y.min,y.max),type="n",xlab="n",ylab="Quantiles")
 for( i in 1:length(ind)){
    lines(n, pred.val[,i],...)
    text(min(n) ,pred.val[1,i], label=names(as.data.frame(pred.val))[i], cex=0.8)
     }
}