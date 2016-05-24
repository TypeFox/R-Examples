vuongtest<-function(model1,model2,selection="AIC"){
    ll1<-model1$ll
    ll2<-model2$ll
    m<-ll1-ll2
    t.stat<-sqrt(length(m))*mean(m)/sd(m)
    if (selection=="AIC"){
	t.stat=t.stat -(model1$npar -model2$npar)
	}
    if (selection=="BIC"){
	n=length(ll1)
	t.stat=t.stat -(model1$npar*log(n)/2 - model2$npar*log(n)/2)
	}
    return(t.stat)



}
