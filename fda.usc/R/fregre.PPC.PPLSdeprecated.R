
################
################
fdata2ppls<-function(fdataobj,y,ncomp = 2,lambda=0,P=c(0,0,1),norm=TRUE,...) {
    if (!is.fdata(fdataobj)) fdataobj<-fdataobj(fdataobj)
    C <- match.call()
    X<-fdataobj$data
    tt<-fdataobj[["argvals"]]
    rtt<-fdataobj[["rangeval"]]
    nam<-fdataobj[["names"]]
    J <- ncol(X);n <- nrow(X)
    Jmin<-min(c(J,n))
    Jmax = min(J+1,n-1)
    Beta <- matrix(, J,ncomp)
    W <- V <- Beta
    dW <- dBeta <- dV <- array(dim = c(ncomp, J, n))
    X0 <- X;  y0 <- y
    mean.y <- mean(y)
    y <- scale(y, scale = FALSE)
    center<-fdata.cen(fdataobj)
    mean.X<-center$meanX
    repnX<-rep(1, n)
    X<-center$Xcen$data
    
if (norm)    {
    sd.X <- sqrt(apply(X, 2, var))
     X <- X/(rep(1, nrow(X)) %*% t(sd.X))
    X2<-fdata(X,tt,rtt,nam)
    dcoefficients = NULL
    tX<-t(X)    
    A <- crossprod(X)
    b <- crossprod(X,y)
    if (lambda>0) {
     if (is.vector(P))  {     P<-P.penalty(tt,P)          }
#     else {
          dimp<-dim(P)
          if (!(dimp[1]==dimp[2] & dimp[1]==J))
              stop("Incorrect matrix dimension P")
 #         }
     M <- solve( diag(J) + lambda*P)
     W[, 1]<- M %*%b
     }
     else    {
            M<-NULL
            W[, 1] <- b
            }
    dV[1, , ] <- dW[1, , ] <- dA(W[, 1], A, tX)
    W[, 1] <- W[, 1]/sqrt((sum((W[, 1]) * (A %*% W[,1]))))
    V[, 1] <- W[, 1]
    Beta[, 1] <- sum(V[, 1] * b) * V[, 1]
    dBeta[1, , ] <-dvvtz(V[, 1], b, dV[1, , ],tX)
    if (ncomp>1) {
    for (i in 2:ncomp) {
           vsi <-b - A %*% Beta[, i - 1]
           if (!is.null(M))  vsi<- M %*% vsi
            W[, i]<-vsi
            dW[i, , ] <- t(X) - A %*% dBeta[i - 1, , ]
            V[, i] <- W[, i] - vvtz(V[, 1:(i - 1), drop = FALSE],A %*% W[, i])
            dV[i, , ] = dW[i, , ] - dvvtz(V[, 1:(i - 1),drop = FALSE],
              A %*% W[, i], dV[1:(i - 1), , , drop = FALSE], A %*% dW[i, , ])
            dV[i, , ] <- dA(V[, i], A, dV[i, , ])
            V[, i] <- V[, i]/sqrt((sum(t(V[, i]) %*% A %*% V[,i])))
            Beta[, i] = Beta[, i - 1] + sum(V[, i] * b) * V[,i]
            dBeta[i,,]<-dBeta[i-1,,]+dvvtz(V[,i],b,dV[i,,],tX)
            }
    }
    dcoefficients <- NULL
    dcoefficients <- array(0, dim = c(ncomp + 1, J, n))
    dcoefficients[2:(ncomp + 1), , ] = dBeta
    sigmahat <- RSS <- yhat <- vector(length = ncomp + 1)
    DoF <- 1:(ncomp + 1)
    Yhat <- matrix(, n, ncomp + 1)
    dYhat <- array(dim = c(ncomp + 1, n, n))
    coefficients <- matrix(0, J, ncomp + 1)
    coefficients[, 2:(ncomp + 1)] = Beta/(sd.X %*% t(rep(1, ncomp)))
    intercept <- rep(mean.y, ncomp + 1) - t(coefficients) %*% t(mean.X$data)
    covariance <- array(0, dim = c(ncomp + 1, J, J))
    DD <- diag(1/sd.X)
    for (i in 1:(ncomp + 1)) {
        Yhat[, i] <- X0 %*% coefficients[, i] + intercept[i]
        res <- y0 - Yhat[, i]
        yhat[i] <- sum((Yhat[, i])^2)
        RSS[i] <- sum(res^2)
        dYhat[i, , ]<-X %*%dcoefficients[i, , ] + matrix(1,n, n)/n
        DoF[i] <- sum(diag(dYhat[i, , ]))
#        dummy <- (diag(n) - dYhat[i, , ]) %*% (diag(n) - t(dYhat[i, , ]))
#        sigmahat[i] <- sqrt(RSS[i]/sum(diag(dummy)))
#        if (i > 1) {
#                covariance[i, , ] <- sigmahat[i]^2 * DD %*% dcoefficients[i,
#                  , ] %*% t(dcoefficients[i, , ]) %*% DD
#            }
    }
    V2<- fdata(t(V)*(rep(1, nrow(t(V))) %*% t(sd.X)),tt,rtt,nam)
    V2$data<-sweep(V2$data,1,norm.fdata(V2),"/")
#   V2<-fdata(t(V),tt,rtt,nam)
#   V2$data<-sweep(V2$data,1,norm.fdata(V2),"/")
#    W2<-fdata(t(W),tt,rtt,nam)
#    X3<-fdata(X,tt,rtt,nam)
#    beta.est<-fdata(t(Beta),tt,rtt,nam)
     DoF[DoF > Jmax] = Jmax
#    intercept <- as.vector(intercept)    
      }
    else {
       plsr<-mplsr(X,y,ncomp=ncomp,lambda=lambda,P=P,...)   
       V2<-fdata(t(plsr$loading.weights),tt,rtt,nam)
       X2<-fdata(X,tt,rtt,nam)
       DoF<-1:ncomp
       Yhat<-plsr$fitted.values    
       yhat<-sum(Yhat)^2          
   }    
   scores<-inprod.fdata(X2,V2,...)
 #   TT = X %*% V  ## son los scores
    l<-1:ncomp
    colnames(scores) <- paste("PLS", l, sep = "")
    outlist = list(call=C,df = DoF, rotation=V2,x=scores,lambda=lambda,P=P,
    fdataobj=fdataobj,y=y0,l=l,fdataobj.cen=center$Xcen,mean=mean.X,Yhat = Yhat,yhat = yhat)
#    Yhat = Yhat, coyhat = yhat,efficients = coefficients,intercept = intercept,
#     RSS = RSS,TT=TT, sigmahat = sigmahat,covariance = covariance,
#     W2=W2,beta.est=beta.est,
    class(outlist)<-"fdata.comp"
    return(outlist)
}

###############
################################################################
################################################################
fdata2ppc<-function (fdataobj,  ncomp = 2,norm = TRUE,lambda=0,P=c(0,0,1),...)
{
    C <- match.call()
    if (!is.fdata(fdataobj))
        stop("No fdata class")
    nas1 <- apply(fdataobj$data, 1, count.na)
    if (any(nas1))
        stop("fdataobj contain ", sum(nas1), " curves with some NA value \n")
    X <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    nam <- fdataobj[["names"]]
    mm <- fdata.cen(fdataobj)
    xmean <- mm$meanX
    Xcen.fdata <- mm$Xcen
    dimx<-dim(X)
    n <- dimx[1]
    J <- dimx[2]
    Jmin <- min(c(J, n))
    if (lambda>0) {
     if (is.vector(P))  {     P<-P.penalty(tt,P)          }
#     else {
          dimp<-dim(P)
           if (!(dimp[1]==dimp[2] & dimp[1]==J))
              stop("Incorrect matrix dimension P")
#          }
         M <- solve( diag(J) + lambda*P)
         Xcen.fdata$data<- Xcen.fdata$data %*%t(M)
     }
    eigenres <- svd(Xcen.fdata$data)
    v <- eigenres$v
    u <- eigenres$u
    d <- eigenres$d
    D <- diag(d)
    vs <- fdata(t(v), tt, rtt, list(main = "fdata2pc", xlab = "t",
        ylab = "rotation"))
    scores <- matrix(0, ncol = J, nrow = n)
    if (norm) {
        dtt <- diff(tt)
        drtt <- diff(rtt)
        eps <- as.double(.Machine[[1]] * 10)
        inf <- dtt - eps
        sup <- dtt + eps
        if (all(dtt > inf) & all(dtt < sup))
            delta <- sqrt(drtt/(J - 1))
        else delta <- 1/sqrt(mean(1/dtt))
        no <- norm.fdata(vs)
        vs <- vs/delta
        newd <- d * delta
        scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs, ...)
    }
    else {
        scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs, ...)
        newd <- d
    }
    colnames(scores) <- paste("PC", 1:J, sep = "")
    l <- 1:ncomp
    out <- list(call = C,d = newd, rotation = vs[1:ncomp],x = scores,
    lambda = lambda,P=P, fdataobj.cen = Xcen.fdata,
    mean = xmean, fdataobj = fdataobj,l=l,u=u[,1:ncomp,drop=FALSE])
    class(out) = "fdata.comp"
    return(out)
}


#################################################################
#################################################################
fregre.ppc.cv=function (fdataobj, y, kmax=8,lambda=0,P=c(0,0,1),criteria = "SIC",...) {
  if (class(fdataobj)=="fdata.comp") {
  fdataobj<-fdataobj$fdataobj
  if (min(lambda)!=0 | !is.null(P)) warning("The arguments lambda and P are not used")
  }
  else {if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  tt<-fdataobj[["argvals"]]
  x<-fdataobj[["data"]]
  X<-fdata.cen(x)[[1]]$data
  np<-ncol(x)
  if (min(lambda)>0) {
   if (is.vector(P)) {    P<-P.penalty(tt,P=P)   }
   if (is.logical(lambda[1]))   {
      normx<-sqrt(sum(abs((t(x)%*%x)^2)))
      normp<-sqrt(sum(abs(P^2)))
      normip<-sqrt(sum(abs(diag(np)+P)^2))
      pii <-traza(P)
      lambda0<-(-2*pii+sqrt(4*pii^2-4*(np-normx^2)*normp^2))/(2*normp^2)
      #  lambda1<-normx/normp
      #print(lambda0)
      #print(lambda1)
      lambda<-seq(0,sqrt(lambda0),len=10)
    }
  }
  }
#  pc<-fdata2ppc(fdataobj,ncomp=kmax,lambda=lambda,P=P,...)
if (is.null(names(y))) names(y)<-1:length(y)
rtt<-fdataobj[["rangeval"]]
n <- nrow(x);    #nc <- ncol(x)
cv.opt1 = Inf;    pc.opt1 = NA
c1 = matrix(1:kmax, nrow = 1)
num.pc = nrow(c1)
l = l2 = list()
max.c = length(c1)
c0 = 1:kmax
use = rep(FALSE, kmax)
tab = list("AIC", "AICc","SIC", "SICc","rho","CV")
type.i = pmatch(criteria, tab)
#pc2<-pc
lenrn<-length(lambda)
MSC3<-list()
pc.opt2 <- matrix(NA,nrow=lenrn,ncol=kmax)
rownames(pc.opt2)<-paste("lambda=",zapsmall(signif(lambda)),sep="")
colnames(pc.opt2)<-paste("PC(",1:kmax,")",sep="")
MSC2<-pc.opt2
MSC.min<-Inf
min.rn<-lambda[1]
if (is.na(type.i))     stop("Error: incorrect criteria")
else {
   if (type.i < 6) {
   for (r in 1:lenrn) {
       pc<-fdata2ppc(fdataobj,ncomp=kmax,lambda=lambda[r],P=P,...)
       cv.opt1 = Inf
       pc.opt1 = NA
       l = l2 = list()
       c1 = matrix(1:kmax, nrow = 1)
       num.pc = nrow(c1)
       max.c = length(c1)
       c0 = 1:kmax
       use = rep(FALSE, kmax)
       pc2<-pc
    for (k in 1:kmax) {
                cv.AIC <- rep(NA, max.c)
                for (j in 1:max.c) {
                  pc2$rotation <- pc$rotation#[c1[, j]]
                  pc2$l <- pc$l[c1[, j]]
                  out = fregre.ppc(pc2, y,l=c1[, j],...)
                  ck<-out$df
                  s2 <- sum(out$residuals^2)/n
                  cv.AIC[j]<-switch(criteria,
              "AIC"=log(s2) + 2 * (ck)/n,
              "AICc"=log(s2) + 2 * (ck)/(n - ck - 2),
              "SIC"=log(s2) + log(n) * ck/n,
              "SICc"=log(s2) + log(n) * ck/(n-ck-2),
              "HQIC"=log(s2) + 2*log(log(n)) * ck/n,
              "rho"={A<-out$residuals;B<-1-diag(out$H)/n; D1<-(A/B)^2;sum(D1)})
                }
                min.AIC = min(cv.AIC)
                pc.opt1 <- c1[, which.min(cv.AIC)]
                l[[k]] = pc.opt1[k]
                l2[[k]] = min.AIC
                use[pc.opt1[k]] = TRUE
                l[[k + 1]] = c0[use == FALSE]
                c1 = t(expand.grid(l))
                ck = nrow(c1) + 1
                max.c = ncol(c1)
            }
     mn = which.min(l2)
     MSC = as.numeric(l2)
     if ( MSC.min>MSC[mn]) {
       min.rn<-r
       MSC.min = MSC[mn]
       pc.opt3<-pc.opt1[1:mn]
       }
    pc.opt = pc.opt1[1:mn]
    MSC2[r,]<-MSC
    pc.opt2[r,]<-pc.opt1
    }
    }
#### CV criteria
    else {
    pb=txtProgressBar(min=0,max=lenrn,width=50,style=3)
    pcl<-list()
    for (r in 1:lenrn) {
     setTxtProgressBar(pb,r-0.5)
     for (i in 1:n) {
      pcl[[i]]<-fdata2ppc(fdataobj[-i,],ncomp=kmax,lambda=lambda[r],P=P,...)
      }
     cv.opt1 = Inf
     pc.opt1 = NA
     l = l2 = list()
     c1 = matrix(1:kmax, nrow = 1)
     num.pc = nrow(c1)
     max.c = length(c1)
     c0 = 1:kmax
     use = rep(FALSE, kmax)
     for (k in 1:kmax) {
      cv.AIC <- rep(NA, max.c)
      cv.AIC2 <- matrix(NA,nrow=max.c,ncol=lenrn)
      rownames(cv.AIC2)<-1:max.c
      for (j in 1:max.c) {
        residuals2<- rep(NA, n)
        maxk<-max(c1[, j])
          for (i in 1:n){
            pc2<-pcl[[i]]
            pc2$rotation<-pcl[[i]]$rotation#[c1[,j]]
            pc2$l<-pcl[[i]]$l[c1[,j]]
            out = fregre.pc(pc2,y[-i],l=c1[,j],...) #####
#            out = fregre.pc(fdataobj[-i,],y[-i],l=c1[,j],rn=rn[r],...)
            ck<-out$df
            residuals2[i] <- ((y[i] - predict(out,fdataobj[i,]))/(n-ck))^2
            }
          cv.AIC[j] <-sum(residuals2)/n
          }
                min.AIC = min(cv.AIC)
                pc.opt1 <- c1[, which.min(cv.AIC)]
                l[[k]] = pc.opt1[k]
                l2[[k]] = min.AIC
                use[pc.opt1[k]] = TRUE
                l[[k + 1]] = c0[use == FALSE]
                c1 = t(expand.grid(l))
                ck = nrow(c1) + 1
                max.c = ncol(c1)
            }
     mn = which.min(l2)
     MSC = as.numeric(l2)
     if ( MSC.min>MSC[mn]) {
       min.rn<-r
       MSC.min = MSC[mn]
       pc.opt3<-pc.opt1[1:mn]
       }
    pc.order<-names(MSC)
    pc.opt = pc.opt1[1:mn]
    MSC2[r,]<-MSC
    pc.opt2[r,]<-pc.opt1
    setTxtProgressBar(pb,r)
    }
    close(pb)
    }
 mn = which.min(l2)
 MSC = as.numeric(l2)
 names(pc.opt3)<-paste("PC", pc.opt3, sep = "")
 rn.opt<-lambda[min.rn]
 fregre=fregre.ppc(fdataobj,y,l=drop(pc.opt),lambda=rn.opt,P=P,...)
 return(list("fregre.ppc"=fregre,pc.opt = pc.opt3,lambda.opt=rn.opt,
 PC.order=pc.opt2,MSC.order=MSC2))
}
}
#################################################################
#################################################################

####################################################################
####################################################################
fregre.ppc=function (fdataobj, y, l =NULL,lambda=0,P=c(0,0,1),...){
if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    if (is.null(l))    {
       l<-pc$l
       }
    else if (length(l)>nrow(pc$rotation)) stop("Incorrect value for  argument l")
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    if (lambda!=0 | !is.null(P)) warning("The arguments lambda and P are not used")
   }
else {
 if (is.null(l)) l<- 1:3
 if (!is.fdata(fdataobj))    fdataobj = fdata(fdataobj)
# omit<-omit.fdata(fdataobj,y)
# fdataobj<-omit[[1]]
# y<-omit[[2]]
  tt<-fdataobj[["argvals"]]
  x<-fdataobj[["data"]]
  X<-fdata.cen(x)[[1]]$data
  np<-ncol(X)
  pc<-fdata2ppc(fdataobj,ncomp=max(l),lambda=lambda,P=P,...)
  }
########################3
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x); np <- ncol(x);lenl = length(l)
    if (is.null(rownames(x)))        rownames(x) <- 1:n
    X <-xcen<-pc$fdataobj.cen
    if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
#    if (length(l) == 1)     vs<-pc$rotation$data[l,]
#    else                    vs<-(pc$rotation$data[l,])
    scores<-Z<-(pc$x[,l])
    cnames<-colnames(pc$x)[l]
#    df<-lenl+1
    J<-min(np,lenl)
    ymean<-mean(y)
    ycen<- y - ymean
###########################
    response = "y"
    dataf<-data.frame(y,Z)
    colnames(dataf)<-c("y",cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf,"+",cnames[i],sep="")
    object.lm = lm(formula = pf, data =dataf , x = TRUE,y = TRUE)
    beta.est<-object.lm$coefficients[2:(lenl+1)]*pc$rotation[l]
    beta.est$data<-apply(beta.est$data,2,sum)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1) 
    Z=cbind(rep(1,len=n),Z)
    S=solve(t(Z)%*%Z)
    H<-Z%*%S%*%t(Z)
    e<-object.lm$residuals
    df<-traza(H)
    rdf<-n-df
    sr2 <- sum(e^2)/rdf
    Vp<-sr2*S 
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
    GCV <- sum(e^2)/rdf^2
 out <- list(call = C, beta.est = beta.est,coefficients=object.lm$coefficients,
 fitted.values =object.lm$fitted.values, residuals = object.lm$residuals,
 H=H,df = df,r2=r2, GCV=GCV,sr2 = sr2, Vp=Vp,l = l,lambda=lambda,P=P, fdata.comp=pc,
 lm=object.lm,fdataobj = fdataobj,y = y)
 #pc=pc,pf = pf,Z=Z
 class(out) = "fregre.fd"
 return(out)
}
####################################################################
####################################################################



#################################################################
#################################################################
fregre.ppls=function(fdataobj, y=NULL, l = NULL,lambda=0,P=c(0,0,1),...){
if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    if  (is.null(l)) l<-1:nrow(pc$rotation)
    if  (is.null(y)) y<-pc$y
    else if (all(y!=pc$y)) warning("y is different from that calculated on the pls basis")
   }
else {
 if (is.null(l)) l<- 1:3
# omit<-omit.fdata(fdataobj,y)
# fdataobj<-omit[[1]]
# y<-omit[[2]]
   pc<-fdata2ppls(fdataobj,y,ncomp=max(l),lambda=lambda,P=P,...)
}
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x); np <- ncol(x);lenl = length(l)
    if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    if (is.null(rownames(x)))        rownames(x) <- 1:n
    ycen = y - mean(y)
    if (length(l) == 1)      vs <- pc$rotation$data
    else                     vs <- t(pc$rotation$data)
    Z<-(pc$x[,l])
    xcen<-pc$fdataobj.cen
    cnames<-colnames(pc$x)[l]
    response = "y"
    df<-data.frame(y,Z)
    colnames(df)<-c("y",cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf,"+",cnames[i],sep="")
    object.lm = lm(formula = pf, data =df , x = TRUE,y = TRUE)
    beta.est<-object.lm$coefficients[2:(lenl+1)]*pc$rotation[l]
    beta.est$data<-apply(beta.est$data,2,sum)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
            if  (pc$type=="pls") {
             if (pc$norm)  {
              sd.X <- sqrt(apply(fdataobj$data, 2, var))
              beta.est$data<-  beta.est$data/sd.X
             }      
            }     
#    H<-diag(hat(Z, intercept = TRUE),ncol=n)
 # H2<-lm.influence(object.lm, do.coef = T)$hat# o bien
 #    I <- diag(1/(n*pc$lambdas[l]), ncol = lenl) #1/n
    Z=cbind(rep(1,len=n),Z)
    S=solve(t(Z)%*%Z)
    H<-Z%*%S%*%t(Z)
    e<-object.lm$residuals
#    df = traza(H)
   df<-pc$df[lenl]+1
    rdf<-n-df
    sr2 <- sum(e^2)/rdf
    Vp<-sr2*S 
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
    GCV <- sum(e^2)/rdf^2
 out <- list(call = C, beta.est = beta.est,coefficients=object.lm$coefficients,
 fitted.values =object.lm$fitted.values, residuals = object.lm$residuals,
 H=H,df = df,r2=r2, GCV=GCV,sr2 = sr2, Vp=Vp,l = l,lambda=lambda,P=P, fdata.comp=pc,
 lm=object.lm,fdataobj = fdataobj,y = y)
    class(out) = "fregre.fd"
    return(out)
}
#################################################################
#################################################################


#################################################################
#################################################################
fregre.ppls.cv=function (fdataobj, y, kmax=8,lambda=0,P=c(0,0,1),
 criteria = "SIC",...) {
if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    kmax<-nrow(pc$basis)
   }
else {
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 omit<-omit.fdata(fdataobj,y)
 fdataobj<-omit[[1]]
 y<-omit[[2]]
}
  x<-fdataobj[["data"]]
  tt<-fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  n <- nrow(x);    nc <- ncol(x)
  if (min(lambda)>0) {
   if (is.vector(P)) {    P<-P.penalty(tt,P=P)   }
   if (is.logical(lambda[1]))   {
      normx<-sqrt(sum(abs((t(x)%*%x)^2)))
      normp<-sqrt(sum(abs(P^2)))
      normip<-sqrt(sum(abs(diag(nc)+P)^2))
      pii <-traza(P)
      lambda0<-(-2*pii+sqrt(4*pii^2-4*(nc-normx^2)*normp^2))/(2*normp^2)
      #  lambda1<-normx/normp
      #print(lambda0)
      #print(lambda1)
      lambda<-seq(0,sqrt(lambda0),len=10)
    }
  }
    lenrn<-length(lambda)
    ind =1:kmax
    l = l2 = list()
    ck = 1
    tab = list("AIC", "AICc","SIC", "SICc","HQIC","rho","CV")
    type.i = pmatch(criteria, tab)
    MSC.min<-Inf
    cv.AIC <- matrix(NA,nrow=lenrn,ncol=kmax)
    if (is.na(type.i))     stop("Error: incorrect criteria")
    else {
    if (type.i < 7) {
#        cv.AIC <- rep(NA, kmax)
      for (r in 1:lenrn) {
        pls<-fdata2ppls(fdataobj,y,ncomp=kmax,lambda=lambda[r],P=P,...)
        for (j in 1:kmax) {
            pls2<-pls
            pls2$rotation<-pls$rotation[1:j]
            out = fregre.ppls(pls2,y,lambda=lambda[r],P=P,...)
            ck<-out$df
            s2 <- sum(out$residuals^2)/n  #(n-ck)
            cv.AIC[r,j]<-switch(criteria,
              "AIC"=log(s2) + 2 * (ck)/n,
              "AICc"=log(s2) + 2 * (ck)/(n - ck - 2),
              "SIC"=log(s2) + log(n) * ck/n,
              "SICc"=log(s2) + log(n) * ck/(n-ck-2),
              "HQIC"=log(s2) + 2*log(log(n)) * ck/n,
              "rho"={A<-out$residuals;B<-1-diag(out$H)/n; D1<-(A/B)^2;sum(D1)})
     if ( MSC.min>cv.AIC[r,j]) {
       rn.opt<-r
       pc.opt<-j
       MSC.min= cv.AIC[r,j]
       }
    }
#    min.AIC = min(cv.AIC)
#    pc.opt <- which.min(cv.AIC)
    }
    }
# CV criteria
    else {
#        pc2<-pc
        for (j in 1:kmax) {
         residuals2<-rep(NA,n)
         for (r in 1:lenrn) {
          for (i in 1:n){
            out = fregre.ppls(fdataobj[-i], y[-i],lambda=lambda[r],P=P,...)
            ck<-out$df
            a1<-out$coefficients[1]
            out$beta.est$data<-matrix(out$beta.est$data,nrow=1)
            b1<-inprod.fdata(fdata.cen(fdataobj[i],out$fdata.comp$mean)[[1]],out$beta.est)
            yp<- a1+b1
#            residuals[i] <- y[i] - yp
            residuals2[i] <- ((y[i] - yp)/(n-ck))^2
            }
#         cv.AIC[j] <- mean(residuals^2)/(n-j)^2###
         cv.AIC[r,j] <-sum(residuals2)/n

     if ( MSC.min>cv.AIC[r,j]) {
       rn.opt<-r
       pc.opt<-j
       MSC.min= cv.AIC[r,j]
       }
    }   }
    }   }
    colnames(cv.AIC) = paste("PLS",1:kmax , sep = "")
    rownames(cv.AIC) = paste("lambda=",lambda , sep = "")
#    pc2$basis<-pc$rotation[1:pc.opt]
    fregre=fregre.ppls(fdataobj,y,l=1:pc.opt,lambda=rn.opt,P=P,...)
    MSC.min = cv.AIC[rn.opt,pc.opt]
    return(list("fregre.pls"=fregre,pls.opt = 1:pc.opt,lambda.opt=lambda[rn.opt],
    MSC.min = MSC.min,MSC = cv.AIC))
}
#################################################################
#################################################################


