###### within 2sls
ivplm.w2sls <- function(Y,X,H = NULL, endog = NULL, lag=FALSE, listw, lag.instruments,T,N,NT, twow = FALSE, listw2 = NULL){

indic <- rep(1:N,T)
##transform y and X	
ywithin <-panel.transformations(Y,indic, type= "within")
Xwithin <- panel.transformations(X, indic, type= "within")

colnames(Xwithin)<-colnames(X)
del <- which(diag(var(Xwithin)) == 0)
Xwithin <- Xwithin[,-del]
# print(Xwithin[1:5,])	


if(!lag){
##transform the instruments H
Hwithin <-panel.transformations(H, indic, type= "within")

if(lag.instruments) {
	
	L.Hwithin <- as.matrix(listw %*% Hwithin)
	L2.Hwithin <- as.matrix(listw %*% L.Hwithin)
	Hwithin <- cbind(Hwithin, as.matrix(L.Hwithin), as.matrix(L2.Hwithin))
	
}


##transform the endogenous variables endog
endogwithin <-panel.transformations(endog, indic, type= "within")
colnames(endogwithin)<-colnames(endog)

res<-spgm.tsls(as.matrix(ywithin), as.matrix(endogwithin), Xwithin, as.matrix(Hwithin))
varb<-res$var *res$df /((N * (T -1)) - ncol(as.matrix(Xwithin)) - ncol(endogwithin)) 
res$var<-varb
sigma2v1<- res$sse/ ((N * (T -1)) - ncol(as.matrix(Xwithin)) - ncol(endogwithin)) 
res$sigmav<- sigma2v1	
res$Hwithin <- Hwithin


	}
	
else{

   wywithin <- listw %*% as.matrix(ywithin)
   wywithin <- as.matrix(wywithin)
   colnames(wywithin)<-"lambda"

if(is.null(endog)){

            if(twow){
            	
	WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
	W2Xwithin <- listw2 %*%  Xwithin
    W2WXwithin <- listw2 %*% WXwithin
    W2WWXwithin <- listw2 %*% WWXwithin
            	
 	Hwithin <-cbind(as.matrix(WXwithin), as.matrix(WWXwithin), as.matrix(W2Xwithin), as.matrix(W2WXwithin), as.matrix(W2WWXwithin))            	
    
            }
else{            
	
	WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
 	Hwithin <-cbind(as.matrix(WXwithin), as.matrix(WWXwithin))
 	
 	}


res<-spgm.tsls(ywithin, wywithin, Xwithin, Hwithin)
varb<-res$var *res$df / ((N * (T -1)) - ncol(as.matrix(Xwithin)) - 1) 
res$var<-varb
sigma2v1<- res$sse / ((N * (T -1)) - ncol(as.matrix(Xwithin)) - 1) 
res$sigmav <- sigma2v1
res$Hwithin <- Hwithin
		}
		
else{
			
			Hwithin <-panel.transformations(H, indic, type= "within")
			
if(lag.instruments ) {
	
	L.Hwithin <- listw %*% Hwithin
	L2.Hwithin <- listw %*% L.Hwithin

	if(twow){
		
		w2.Hwithin <- as.matrix(listw2 %*% Hwithin)
		w2w.Hwithin <- as.matrix(listw2 %*% L.Hwithin)
		w2ww.Hwithin <- as.matrix(listw2 %*% L2.Hwithin)
	Hwithin <- cbind(Hwithin, as.matrix(L.Hwithin), as.matrix(L2.Hwithin), w2.Hwithin, w2w.Hwithin, w2ww.Hwithin)		
	
	}
	
	else 	Hwithin <- cbind(Hwithin, as.matrix(L.Hwithin), as.matrix(L2.Hwithin))
}


            if(twow){
            	
	WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
	W2Xwithin <- listw2 %*%  Xwithin
    W2WXwithin <- listw2 %*% WXwithin
    W2WWXwithin <- listw2 %*% WWXwithin
            	
 	Hwithin <-cbind(Hwithin, as.matrix(WXwithin), as.matrix(WWXwithin), as.matrix(W2Xwithin), as.matrix(W2WXwithin), as.matrix(W2WWXwithin))            	
    
            }
else{            
	
	WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
 	Hwithin <-cbind(Hwithin, as.matrix(WXwithin), as.matrix(WWXwithin))
 	
 	}


endogwithin <- panel.transformations(endog, indic, type= "within")

endogwithin <-cbind(endogwithin, wywithin)
colnames(endogwithin)<-c(colnames(endog), "lambda")
# colnames(Xwithin)<-colnames(X)[-del]

res<-spgm.tsls(ywithin, endogwithin, Xwithin, Hwithin)

varb<-res$var *res$df / ((N * (T -1)) - ncol(as.matrix(Xwithin)) - ncol(endogwithin)) 
res$var<-varb
sigma2v1<- res$sse / ((N * (T -1)) - ncol(as.matrix(Xwithin)) - ncol(endogwithin)) 
res$sigmav <- sigma2v1
res$Hwithin <- Hwithin
	}		
	
	}	
res
}




