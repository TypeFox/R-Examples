StackingXY <-
function(YList,XList,RespDist,Beta,VT,B=NULL,ZZMatrix,LinkList) {
    
    for (i in 1:length(YList)){
        if (i==1) {
            Y<-YList[[1]]
            X<-XList[[1]]
        }    
        else {
            Y<-rbind(Y,YList[[i]])
            X<-dbind(X,XList[[i]])
        }
   }
   if (is.null(B)) B<-rep(1,length(Y))
   nModels <- length(YList)
   ModelsDims <- sapply(YList,nrow)
   cModelsDims <- cumsum(c(0,ModelsDims))
   mu<-0
   Wvec<-0
   dWdmu<-0
   dmudeta<-0
   d2Wdmu2<-0
   
   # Create the TT matrix #
   TT <- cbind(X,ZZMatrix)
   nrand <- ncol(ZZMatrix)
   ptot <- ncol(X)
   ntot <- nrow(Y)
   TT <- rbind(TT, cbind(matrix(0,nrand,ptot),diag(nrand)))
   # Create the output #
   eta<-TT[1:ntot,]%*%as.matrix(c(Beta,unlist(VT)))
   for (i in 1:nModels){
        mu[(cModelsDims[i]+1):cModelsDims[i+1]]<-B[(cModelsDims[i]+1):cModelsDims[i+1]]*InvLinkY(eta[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]])
        Wvec[(cModelsDims[i]+1):cModelsDims[i+1]]<-Wmatgen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
        dWdmu[(cModelsDims[i]+1):cModelsDims[i+1]]<-dWdmugen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
        d2Wdmu2[(cModelsDims[i]+1):cModelsDims[i+1]]<-d2Wdmu2gen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
        dmudeta[(cModelsDims[i]+1):cModelsDims[i+1]]<-dmudetagen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
    }   
   
    return(list(TT=TT),mu=mu,Wvec=Wvec,dWdmu=dWdmu,d2Wdmu2,dmudeta=dmudeta)
}
