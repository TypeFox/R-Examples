hP <- function(x,theta) theta^x - 1
hD <- function(x) (x^2)/2
hG <- function(x,theta) -log(theta + x) + log(theta)

Xclosedp <- function(t, m, h, theta, histpos, nbcap, mX = NULL)
{
  if (missing(histpos))   histpos <- histpos.t(t)
  if (missing(nbcap))     nbcap <- rowSums(histpos)
  
  if (is.null(mX)) {
    if (m=="none") {
      mX <- as.matrix(ifelse(rowSums(histpos)>0,1,0))
      cnames <- "beta"
    }
    if (m %in% c("M0", "Mh")) {
      mX <- as.matrix(nbcap)
      cnames <- "beta"
    }
    if (m %in% c("Mt", "Mth")) {
      mX <- histpos
      cnames <- paste("beta",1:t,sep="")
    }
    if (m=="Mb") {
      nbcap_av <- rep(t-t:1,2^(t:1-1))  # nbre d'occasions de capture avant la premiere capture
      nbcap_ap <- (nbcap-1) # nbre de capture apres la premiere capture
      mX <- cbind(nbcap_av,nbcap_ap)
      cnames <- c("beta1","beta2")
    }
    if (m=="Mbh") {
      inv_c1<-1-histpos[,1]
      nbcap_av <- rep(t-t:1,2^(t:1-1))  # nbre d'occasions de capture avant la premiere capture
      nbcap_ap <- (nbcap-1) # nbre de capture après la premiere capture
      mX <- cbind(inv_c1,nbcap_av,nbcap_ap)
      cnames <- c("eta","beta1","beta2")
    }
  } else {
    # On ne touche pas à mX, c'est celui donné en entrée que l'on prend
    cnames <- if(ncol(mX) > 0 && (is.null(colnames(mX)) || any(colnames(mX)=="")))  {
                 paste("mX.",1:ncol(mX),sep="") 
              } else {
                colnames(mX)
              }
  }
  if (is.function(h) || !is.null(h) && h != "Normal") { ## Colonnes pour l'hétérogénéité au besoin
    if (is.function(h)) {
      mX2 <- h(nbcap)
      cnames2 <- "tau"
    } else {
      if (h %in% c("Chao", "LB")) {
        if (t < 3) {
          mX2 <- rep(0, nrow(mX))
          cnames2 <- "eta"
        } else {
          mX2 <- matrix(0,dim(histpos)[1],t-2)
          for (j in (3:t)) { mX2[,j-2]<-pmax(nbcap-j+1,0) }
          cnames2 <- paste("eta",3:t,sep="")
        }
      } else {
        if (h=="Poisson")  mX2 <- hP(nbcap,theta) else
        if (h=="Darroch")  mX2 <- hD(nbcap) else
        if (h=="Gamma")    mX2 <- hG(nbcap,theta)
        cnames2 <- "tau"
      }
    }
    mX <- cbind(mX,mX2)
    cnames <- c(cnames,cnames2)
  }
  
  dimnames(mX)<-NULL
  
  list(mat=mX, coeffnames=cnames, nbcoeff=ncol(mX))
}


Xomega <- function(vt,vm,vh,vtheta,fct.call,typet,histpos)
{
  if (missing(histpos)) {
    # matrice des historiques de captures possibles pour le nombre d'occasions de capture total
    histpos <- if (typet)  histpos.t(sum(na.rm=TRUE,vt)) else histpos.0(vt)
  }
  I <- length(vt) # nombre de periodes primaires
  M <- cnames <- NULL
  nbcoeff <- rep(0,I)
  models <- rep(0,I)
  
  # on cree la matrice periode par periode en respectant les models demandes en entree
  for (i in (1:I))
  {
    # selection des colonnes correspondantes a la periode etudiee dans cette boucle
    if (typet) {
      if (i==1) { histposp <- histpos[,c(1:vt[i])]   
      } else { histposp <- histpos[,c((sum(na.rm=TRUE,vt[1:(i-1)])+1):sum(na.rm=TRUE,vt[1:i]))] }
    } else  histposp <- histpos[,i,drop=FALSE]
    
    Xclosedp.out<-Xclosedp(vt[i],vm[i],vh[[i]],vtheta[i],histposp)
    mXp<-Xclosedp.out$mat
    M <- cbind(M,mXp)
    cnames<-c(cnames,paste(Xclosedp.out$coeffnames,".",i,sep=""))
    nbcoeff[i] <- dim(mXp)[2]
    models[i] <- if (vm[i]%in%c("none","M0","Mt")) vm[i] else
        if (is.function(vh[[i]])) {
          if(length(fct.call$vh)==1)  paste(vm[i],deparse(fct.call$vh))
          else paste(vm[i],fct.call$vh[min(i+1,length(fct.call$vh))])
        } else if(vh[[i]]=="Poisson"||vh[[i]]=="Gamma") paste(vm[i],paste(vh[[i]],vtheta[i],sep="")) else paste(vm[i],vh[[i]])
  }
  
  list(mat=M,models=models,coeffnames=cnames,nbcoeff=nbcoeff)
}



Zdelta <- function (Xdelta)
{
  Xdelta <- as.matrix(Xdelta)
  I <-dim(Xdelta)[2]
  Z <- matrix(0,dim(Xdelta)[1],2*(I-1))
  Z[,1] <- (1-Xdelta[,1])
  if(I>2) {
    i <- 2
    for  (j in (2:(I-1))) {
      Z[,j] <- Z[,(j-1)]*(1-Xdelta[,i])
      i<- i+1
    }
  }
  Z[,I] <- (1-Xdelta[,I])
  if(I>2) {
    i <- 1
    for  (j in ((I+1):(2*(I-1)))) {
      Z[,j] <- Z[,(j-1)]*(1-Xdelta[,(I-i)])
      i <- i+1
    }
  }
  return(Z)
}
