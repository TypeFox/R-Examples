# -------------------------------------------------
# Funcion para convertir una TC con DP a dataframe
# Cada celda es una fila con cuatro factores
# L J I K
# entra tab, rbl,cbl y iden
# iden vector con el numero de caracteres por
#      cada factor para identifcar la fila
# factor L names(rbl)
# factor J names(cbl)
# factor I rownames(tab)
# factor K colnames(tab)
# Campo Elias Pardo
# julio 24 de 2010
#-------------------------------------------------
ctdp2df<-function(tab,rbl,cbl,iden=rep(3,4))
  {
    n <- nrow(tab)
    m <- ncol(tab)
    l <- length(rbl)
    j <- length(cbl)
    freq <- numeric(n*m)
    freq <- unlist(tab)
    K.fac <- factor(rep(1:m,each=n),labels=colnames(tab))
    L.fac <- factor(rep(rep(1:l,rbl,each=TRUE),m),labels=names(rbl))
    J.fac <- factor(rep(rep(1:j,cbl,each=TRUE),each=n),labels=names(cbl))
    I.fac <- factor(rep(1:n,m),labels=rownames(tab))
    df <- data.frame(J=J.fac,K=K.fac,L=L.fac,I=I.fac,freq)
    rn <-paste(substr(df[,1],1,iden[1]),substr(df[,2],1,iden[2]),sep=".")
    rn <-paste(rn,substr(df[,3],1,iden[3]),sep=".")
    rownames(df) <- paste(rn,substr(df[,4],1,iden[4]),sep=".")
    return(df)
  }
# fin de funcion
