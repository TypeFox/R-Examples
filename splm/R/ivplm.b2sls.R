###### between 2sls

ivplm.b2sls <- function(Y,X,H = NULL, endog = NULL, lag=FALSE, listw, lag.instruments, T = T, N = N, NT = NT, twow = FALSE, listw2 = NULL){

indic <- rep(1:N,T)

##transform y	
ybetween<-panel.transformations(Y,indic, type= "between")
ndim <- length(ybetween)
listwnn <- listw[1:ndim, 1:ndim]

Xbetween<-panel.transformations(X,indic, type= "between")
colnames(Xbetween)<-colnames(X)

if (colnames(Xbetween)[1] == "(Intercept)") Xbetween<-Xbetween[,-1]
delb<-as.numeric(which(diag(var(Xbetween))==0))
if(length(delb)==0) Xbetween<-Xbetween
else Xbetween<-Xbetween[,-delb]

if (colnames(X)[1] == "(Intercept)") Xbetween<-cbind(1,Xbetween)
colnames(Xbetween)[1]<-"(Intercept)"

if(!lag){
##transform the instruments H and the endogenous variable
	Hbetween<-panel.transformations(H,indic, type= "between")

	if(lag.instruments ) {
	
	L.Hbetween <- listwnn %*% Hbetween
	L2.Hbetween <- listwnn %*% L.Hbetween
	Hbetween <- cbind(Hbetween, as.matrix(L.Hbetween), as.matrix(L2.Hbetween))
}

	endogbetween<-panel.transformations(endog,indic, type= "between")
   colnames(endogbetween)<-colnames(endog)

res <-spgm.tsls(sqrt(T)*as.matrix(ybetween), sqrt(T)*endogbetween, sqrt(T)*Xbetween, sqrt(T)*as.matrix(Hbetween) )
res$Hbetween <- Hbetween
}

else{
	
	wybetween <- listwnn %*% as.matrix(ybetween)
	wybetween <- as.matrix(wybetween)
    colnames(wybetween) <- ("lambda")
	
	if(is.null(endog)){
		
		            if(twow){
		            	
     listw2nn <- listw2[1:ndim, 1:ndim]       	
	WXbetween <- listwnn %*%  Xbetween
    WWXbetween <- listwnn %*% WXbetween
	W2Xbetween <- listw2nn %*%  Xbetween
    W2WXbetween <- listw2nn %*% WXbetween
    W2WWXbetween <- listw2nn %*% WWXbetween
            	
 	Hbetween <-cbind(as.matrix(WXbetween), as.matrix(WWXbetween), as.matrix(W2Xbetween), as.matrix(W2WXbetween), as.matrix(W2WWXbetween))            	
    
            }
else{            
	
            WXbetween <- as.matrix(listwnn %*% Xbetween)
            WWXbetween <- as.matrix(listwnn %*% WXbetween)
        
Hbetween<-cbind(WXbetween, WWXbetween)        
 	
 	}


res<-spgm.tsls(sqrt(T)*as.matrix(ybetween), sqrt(T)*as.matrix(wybetween), sqrt(T)*Xbetween, sqrt(T)*as.matrix(Hbetween) )
res$Hbetween <- Hbetween
		}
		
else{
	
		Hbetween <- panel.transformations(H,indic, type= "between") 
			
if(lag.instruments ) {
	
	L.Hbetween <- listwnn %*% Hbetween
	L2.Hbetween <- listwnn %*% L.Hbetween

	if(twow){
		listw2nn <- listw2[1:ndim, 1:ndim] 
		w2.Hbetween <- as.matrix(listw2nn %*% Hbetween)
		w2w.Hbetween <- as.matrix(listw2nn %*% L.Hbetween)
		w2ww.Hbetween <- as.matrix(listw2nn %*% L2.Hbetween)
	Hbetween <- cbind(Hbetween, as.matrix(L.Hbetween), as.matrix(L2.Hbetween), w2.Hbetween, w2w.Hbetween, w2ww.Hbetween)		
	
	}
	
	else 		Hbetween <- cbind(Hbetween, as.matrix(L.Hbetween), as.matrix(L2.Hbetween))
}


            if(twow){
    listw2nn <- listw2[1:ndim, 1:ndim]        	
	WXbetween <- listwnn %*%  Xbetween
    WWXbetween <- listwnn %*% WXbetween
	W2Xbetween <- listw2nn %*%  Xbetween
    W2WXbetween <- listw2nn %*% WXbetween
    W2WWXbetween <- listw2nn %*% WWXbetween
            	
 	Hbetween <-cbind(Hbetween, as.matrix(WXbetween), as.matrix(WWXbetween), as.matrix(W2Xbetween), as.matrix(W2WXbetween), as.matrix(W2WWXbetween))            	
    
            }
else{            
	
	WXbetween <- listwnn %*%  Xbetween
    WWXbetween <- listwnn %*% WXbetween
 	Hbetween <-cbind(Hbetween, as.matrix(WXbetween), as.matrix(WWXbetween))
 	
 	}


	##transform the endogenous variables endog
	endogbetween<-panel.transformations(endog,indic, type= "between")
	endogbetween<-cbind(endogbetween, wybetween)

colnames(endogbetween)<-c(colnames(endog), "lambda")


	
res<-spgm.tsls(sqrt(T)*as.matrix(ybetween), sqrt(T)*endogbetween, sqrt(T)*Xbetween, sqrt(T)*as.matrix(Hbetween) )
res$Hbetween <- Hbetween
	}			
	}	

res
}




