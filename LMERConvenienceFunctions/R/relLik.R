relLik<-function(x,y,method=c("AIC","BIC"),ndigits=6,...){
        if(method[1]=="AIC"){
                tmp<-c(AIC(x),AIC(y),AIC(x)-AIC(y),exp(abs(AIC(x)-AIC(y))/2))
                names(tmp)<-c("AIC(x)","AIC(y)","diff","relLik")
                return(round(tmp,ndigits))
        }else if(method[1]=="BIC"){
                tmp<-c(BIC(x),BIC(y),BIC(x)-BIC(y),exp(abs(BIC(x)-BIC(y))/2))
                names(tmp)<-c("BIC(x)","BIC(y)","diff","relLik")
                return(round(tmp,ndigits))

        }else{
                stop("unrecognized method.\n")
        }
}

