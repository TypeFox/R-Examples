consistencyplot <- function (t,cod,cex.axis=.8,mark.code=TRUE,...) {

 # INPUT:
 # t = years
 # cod = code

 cod <- as.factor(cod)
 n <- nlevels(cod)
 codici <- levels(cod)
 ni <- tapply(t,cod,length)

 dimensioni=cex.axis
  
 time <- c(min(t),max(t))
 code <- c(1,n)
 plot(time,code,type="n",axes=FALSE,...)
 axis(1)
 if (mark.code==TRUE) axis(2,at=c(n:1),labels=codici,las=1,cex.axis=dimensioni)
 box(); grid()
 for (i in n:1) {
  ti <- t[cod==codici[i]]
  points(ti,rep(n-i+1,ni[i]),pch=45)
 }
}


# --------------------------------------------------------------------------------- #

serieplot <- function (x,t,lim.x=c(min(x),max(x)),lim.t=c(min(t),max(t)),...) {

 X=x
 T=t
 x=lim.x
 t=lim.t
 plot(t, x, type="n", ...)
 lines(T,X, type="h")
}
