regsubsets2aic <-
function(x,y,z){
      # Calculate AIC and cAIC for every model stored in "z", the
      # output of "regsubsets". This version sorts the output by AICc.
      # Requires the same x and y objects inputted to "regsubsets"
      
      #If intercept
      if (colnames(as.matrix(summary(z)$which))[1]=="(Intercept)") {
      rows<-1:length(summary(z)$bic)
      for(i in 1:nrow(summary(z)$which)){
            
            vars<-which(summary(z)$which[i,-1])
            newdat<-cbind.data.frame(y,x[,vars])
            znew<-lm(y~.,data=newdat)
            L<-logLik(znew)
            aL<-attributes(L)
            z$model[i]<-paste(colnames(summary(z)$which)[-1][vars],collapse=" + ")
            z$p[i]<-as.numeric(rownames(summary(z)$which)[i])+1  # terms in model, including intercept
            z$k[i]<-aL$df      # parameters (including sigma^2)
            z$AIC[i]<- -2*L[1] + 2*aL$df
            z$AICc[i]<- -2*L[1] + 2*aL$df + 
                  (2*aL$df*(aL$df+1))/(aL$nobs-aL$df-1)
            z$BIC[i]<- -2*L[1] + log(aL$nobs)*aL$df
      }
      z$AIC<-z$AIC-min(z$AIC)
      z$AICc<-z$AICc-min(z$AICc)
      z$BIC<-z$BIC-min(z$BIC)
      result<-do.call("cbind.data.frame",
                      list(z[c("model","p","k","AIC","AICc","BIC")],
                           stringsAsFactors=FALSE))
      result<-result[order(result$AICc),]
      rownames(result)<-1:nrow(result)
      newdat<-cbind.data.frame(y,x)
      return(result)
}

else {
      rows<-1:length(summary(z)$bic)
      for(i in 1:nrow(summary(z)$which)){
            vars<-which(summary(z)$which[i,])
            newdat<-cbind.data.frame(y,x[,vars])
            znew<-lm(y~0+.,data=newdat)
            L<-logLik(znew)
            aL<-attributes(L)
            z$model[i]<-paste(colnames(summary(z)$which)[vars],collapse=" + ")
            z$p[i]<-as.numeric(rownames(summary(z)$which)[i])  # terms in model, no intercept
            z$k[i]<-aL$df      # parameters (including sigma^2)
            z$AIC[i]<- -2*L[1] + 2*aL$df
            z$AICc[i]<- -2*L[1] + 2*aL$df + 
                  (2*aL$df*(aL$df+1))/(aL$nobs-aL$df-1)
            z$BIC[i]<- -2*L[1] + log(aL$nobs)*aL$df
      }
      z$AIC<-z$AIC-min(z$AIC)
      z$AICc<-z$AICc-min(z$AICc)
      z$BIC<-z$BIC-min(z$BIC)
      result<-do.call("cbind.data.frame",
                      list(z[c("model","p","k","AIC","AICc","BIC")],
                           stringsAsFactors=FALSE))
      result<-result[order(result$AICc),]
      rownames(result)<-1:nrow(result)
      newdat<-cbind.data.frame(y,x)
      return(result)
}
}
