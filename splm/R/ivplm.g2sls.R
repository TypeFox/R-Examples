# # # # # g2sls

ivplm.g2sls <-function(Y,X,H = NULL,endog = NULL, lag=FALSE, listw, lag.instruments, T = T, N = N, NT = NT){

indic <- rep(1:N,T)
listwnn <- listw[1:N,1:N]
##transform y	
transy<-panel.transformations(Y,indic, type= "both")
ybetween<-transy[[2]]
ywithin<-transy[[1]]
ybetweennt<- rep(ybetween, T)	

##transform X	
transx<-panel.transformations(X,indic, type= "both")
Xbetween<-transx[[2]]
Xwithin<-transx[[1]]
colnames(Xwithin)<-colnames(X)
colnames(Xbetween)<-colnames(X)
Xbetweennt<-matrix(,NT, ncol(Xbetween))
for (i in 1:ncol(Xbetween)) Xbetweennt[,i]<-rep(Xbetween[,i],T)
del<- which(diag(var(Xwithin))==0)
colnames(Xbetweennt)<-colnames(X)

if(!lag){

##transform the instruments H
transH<-panel.transformations(H, indic, type= "both")
Hbetween<-transH[[2]]
Hwithin<-transH[[1]]

if(lag.instruments ) {
	
	L.Hwithin <- listw %*% Hwithin
	L2.Hwithin <- listw %*% L.Hwithin
	Hwithin <- cbind(Hwithin, as.matrix(L.Hwithin), as.matrix(L2.Hwithin))

	L.Hbetween <- listwnn %*% Hbetween
	L2.Hbetween <- listwnn %*% L.Hbetween
	Hbetween <- cbind(Hbetween, as.matrix(L.Hbetween), as.matrix(L2.Hbetween))
	
}

Hbetweennt<-matrix(,NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i], T)

##transform the endogenous variables endog
transendog<-panel.transformations(endog, indic, type= "both")
endogbetween<-transendog[[2]]
endogwithin<-transendog[[1]]
endogbetweennt<-matrix(,NT, ncol(endogbetween))
for (i in 1:ncol(endogbetween)) endogbetweennt[,i] <- rep(endogbetween[,i], T)
colnames(endogbetweennt)<-colnames(endog)
colnames(endogwithin)<-colnames(endog)

#W2SLS
resw<-spgm.tsls(as.matrix(ywithin), endogwithin, Xwithin, Hwithin )

sigma2v1<-resw$sse / ((N * (T -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 

#B2SLS
resb<-spgm.tsls(sqrt(T)*as.matrix(ybetween), sqrt(T)*as.matrix(endogbetween), sqrt(T)*Xbetween, sqrt(T)*as.matrix(Hbetween) )

sigma21<-resb$sse /  resb$df

ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
endogstar<-endogwithin/sqrt(sigma2v1) + endogbetweennt/sqrt(sigma21)


Hstar<-Hwithin/sqrt(sigma2v1) + Hbetweennt/sqrt(sigma21)

res<-spgm.tsls(ystar, endogstar, xstar, Hstar )
res$sigma1<-sigma21
res$sigmav<-sigma2v1

}

else{
       wy<-listw %*% Y
       wy <- as.matrix(wy)
     colnames(wy)<-"lambda"  
	  wywithin <- listw %*% ywithin
     wywithin <- as.matrix(wywithin)
     colnames(wywithin)<-"lambda"
  	  wybetween <- listwnn %*%  as.matrix(ybetween)
     colnames(wybetween) <- "lambda"

           WXwithin <- as.matrix(listw %*% Xwithin)
           WWXwithin <- as.matrix(listw %*% WXwithin)

            WXbetween <- as.matrix(listwnn %*% Xbetween)
            WWXbetween <- as.matrix(listwnn %*% WXbetween)
	
	if(is.null(endog)){

Hwithin<-cbind(WXwithin, WWXwithin)    
    
resw<-spgm.tsls(ywithin, wywithin, Xwithin, Hwithin)

sigma2v1<- resw$sse / ((N * (T -1)) - ncol(as.matrix(Xwithin[,-del])) - 1) 

        
Hbetween<-cbind(WXbetween, WWXbetween)        

resb<-spgm.tsls(sqrt(T)*as.matrix(ybetween), sqrt(T)*as.matrix(wybetween), sqrt(T)*Xbetween, sqrt(T)*as.matrix(Hbetween) )
sigma21<-resb$sse /  resb$df


ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
endogstar<-wywithin/sqrt(sigma2v1) + rep(as.matrix(wybetween), T)/sqrt(as.numeric(sigma21))
endogstar<-as.matrix(endogstar)
colnames(endogstar)<-"lambda"

Hbetweennt<-matrix(,NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i], T)

Hstar<-Hwithin/sqrt(sigma2v1) + Hbetweennt/sqrt(sigma21)

res <- spgm.tsls(ystar, endogstar, xstar, Hstar)

res$sigma1 <- sigma21
res$sigmav <- sigma2v1

}


else{

transH<-panel.transformations(H, indic, type= "both")
Hbetween<-transH[[2]]
Hwithin<-transH[[1]]

if(lag.instruments ) {
	
	L.Hwithin <- as.matrix(listw %*% Hwithin)
	L2.Hwithin <- as.matrix(listw %*% L.Hwithin)
	Hwithin <- cbind(Hwithin, L.Hwithin, L2.Hwithin)

	L.Hbetween <- as.matrix(listwnn %*% Hbetween)
	L2.Hbetween <- as.matrix(listwnn %*% L.Hbetween)
	Hbetween <- cbind(Hbetween, L.Hbetween, L2.Hbetween)
	
}


Hbetweennt<-matrix(, NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i], T)

Hwithin<-cbind(Hwithin, WXwithin, WWXwithin)


transendog<-panel.transformations(endog, indic, type= "both")
endogbetween<-transendog[[2]]
endogwithin<-transendog[[1]]

endogwithin<-cbind(endogwithin, wywithin)
    
resw<-spgm.tsls(as.matrix(ywithin), as.matrix(endogwithin), Xwithin, Hwithin )
sigma2v1<-resw$sse / ((N * (T -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 


Hbetween<-cbind(Hbetween, as.matrix(WXbetween), as.matrix(WWXbetween))
endogbetween<-cbind(endogbetween, as.matrix(wybetween))
endogbetweennt<-matrix(,NT, ncol(endogbetween))
for (i in 1:ncol(endogbetween)) endogbetweennt[,i]<-rep(endogbetween[,i], T)

resb<-spgm.tsls(sqrt(T)*as.matrix(ybetween), sqrt(T)*as.matrix(endogbetween), sqrt(T)*Xbetween, sqrt(T)*as.matrix(Hbetween))
sigma21<-resb$sse / resb$df

ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
endogstar<-endogwithin/sqrt(sigma2v1) + endogbetweennt/sqrt(sigma21)

Hbetweennt<-matrix(,NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i], T)

Hstar<- Hwithin/sqrt(sigma2v1) + Hbetweennt/sqrt(sigma21)

res <- spgm.tsls(ystar, endogstar, xstar, Hstar)

res$sigma1<- sigma21
res$sigma1<- sigma2v1
	}
}


res
}


