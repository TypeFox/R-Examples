"histfreq.t" <- function(X,dfreq)
{
    nl <- dim(X)[1]
    t <- if(dfreq) dim(X)[2]-1 else dim(X)[2]
    vecli <- rep(0,nl)      # vecteur du numero de l'historique (selon l'ordre de histpos.t) de chaque individu
    for (i in (1:nl))  vecli[i] <- 1+sum(na.rm=TRUE,(1-X[i,1:t])*(2^{(t-1):0}))
    X <- X[vecli<2^t,]      # pour mettre de cote les lignes de X comprenant uniquement des zeros s'il y en a
    vecli <- vecli[vecli<2^t] 
    Y <- rep(0,2^t-1)       # on cree un vecteur de 0 de la taille du nbre des differents historiques de capture possibles
    for (i in (1:length(vecli)))  Y[vecli[i]] <- Y[vecli[i]] + if(dfreq) X[i,t+1] else 1
    return(Y)
}


"histfreq.0" <- function(X,dfreq,dtype="hist",vt)
{
  # Note : Il faut absolument que les valeurs de vt soient >= au nombre maximum de captures
  # dans chacune des périodes. Sinon il y a une sortie générée, mais elle n'est pas bonne.  
    I <- length(vt) # nombre de periodes primaires
    t <- sum(vt)
    if (I==1) {
         if (dtype=="hist") { 
               Y <- getfi(X=X,dfreq=dfreq,t=t) 
         } else {
               fi <- if (dfreq) tapply(X[,2],X[,1],sum,na.rm=TRUE) else table(X[,1])
               i <- as.numeric(names(fi))
               Y <- rep(0,t)
               Y[i] <- fi
         } 
         Y <- rev(Y) # On veut la fréquence pour nbcap = t à 1
    } else {
         Xh<-apply(X[,1:vt[1]],1,sum)
         for (i in 2:I)  Xh<-cbind(Xh,apply(X[,c((sum(na.rm=TRUE,vt[1:(i-1)])+1):sum(na.rm=TRUE,vt[1:i]))],1,sum))
         nl<-dim(Xh)[1]
         vecli<-rep(0,nl)    # vecteur du numero de l'historique (selon l'ordre de histpos.0) de chaque individu
         cvt<-c(cumprod((vt+1)[I:2])[(I-1):1],1)
         for (i in (1:nl))  vecli[i]<-1+sum(na.rm=TRUE,(vt-Xh[i,])*cvt)
         X <- X[vecli<prod(vt+1),]      # pour mettre de cote les lignes de X comprenant uniquement des zeros s'il y en a
         vecli<-vecli[vecli<prod(vt+1)]
         Y <- rep(0,prod(vt+1)-1)       # on cree un vecteur de 0 de la taille du nbre  des differents historiques de capture possibles
         for (i in (1:length(vecli)))  Y[vecli[i]] <- Y[vecli[i]] + if(dfreq) X[i,t+1] else 1
    }
    return(Y)
}
