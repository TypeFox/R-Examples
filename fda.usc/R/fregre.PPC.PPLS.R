################################################################
################################################################
D.penalty<- function(tt) {
 rtt<-diff(tt)
 hh<--(1/mean(1/rtt))/rtt
 Dmat1<-diag(hh)
 p<-length(tt)
 mat<-matrix(rep(0,p-1),ncol=1)
 Dmat12<-cbind(Dmat1,mat)
 Dmat22<-cbind(mat,-Dmat1)
 Dmat<-Dmat12+Dmat22
 return(Dmat)
}
################################################################
################################################################
# P.penalty <-function(tt,order=2) {
#   D.pen2<-D.pen <- D.penalty(tt)
#  p<-ncol(D.pen)
#   if (order > 1) {
#      for (k in 2:order) D.pen <- D.pen2[1:(p-k),1:(p-k+1)] %*% D.pen
#   }
#   D.pen <- t(D.pen) %*% D.pen
# }
################################################################
################################################################
P.penalty <- function(tt,P=c(0,0,1)) {
  lenp<-length(P)
  D.pen2<-D.pen <- D.penalty(tt)
  p<-ncol(D.pen)
  D.pen3<-matrix(0,p,p)
  if  (sum(P)!=0){
   D.pen3<-D.pen3+P[1]*diag(p)
   if (lenp > 1) {
   P<-P[-1]
   for (i in 2:lenp-1){
     D.pen2<-D.pen <- D.penalty(tt)
      if (i > 1) {
       for (k in 2:i) D.pen <- D.pen2[1:(p-k),1:(p-k+1)] %*% D.pen
      }
      D.pen3 <- D.pen3+P[i]*(t(D.pen) %*% D.pen)
    }
  }}
  D.pen3
}
################
################
dA<-function (w, A, dw){
    wa <- sqrt(sum((w * (A %*% w))))
    dummy <- (1/wa) * (diag(length(w)) - w %*% t(w) %*% A/(wa^2)) %*%   dw
    return(dummy)
}
################
################
vvtz<-function (v, z){as.vector(v %*% (t(v) %*% z))}
################
################
dvvtz<-function (v, z, dv, dz) {
    if (is.matrix(v) == FALSE) {
        v <- matrix(v, ncol = 1)
        dv <- array(dv, dim = c(1, nrow(dv), ncol(dv)))
    }
    k = ncol(v)
    p <- nrow(v)
    n <- dim(dv)[3]
    dummy <- matrix(0, dim(dv)[2], dim(dv)[3])
    for (i in 1:k) {
        D <- (v[, i] %*% t(z) + sum(v[, i] * z) * diag(p)) %*%
            dv[i, , ] + v[, i] %*% t(v[, i]) %*% dz
        dummy <- dummy + D
    }
    return(dummy)
}
################
################
fdata2pls<-function(fdataobj,y,ncomp = 2,lambda=0,P=c(0,0,1),norm=TRUE,...) {
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
    outlist = list(call=C,df = DoF, rotation=V2,x=scores,lambda=lambda,P=P,norm=norm,type="pls",
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
fdata2pc<-function (fdataobj,  ncomp = 2,norm = TRUE,lambda=0,P=c(0,0,1),...)
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
    lambda = lambda,P=P, fdataobj.cen = Xcen.fdata,norm=norm,type="pc",
    mean = xmean, fdataobj = fdataobj,l=l,u=u[,1:ncomp,drop=FALSE])
    class(out) = "fdata.comp"
    return(out)
}



#################################################################
#################################################################
fregre.pls=function(fdataobj, y=NULL, l = NULL,lambda=0,P=c(0,0,1),...){
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
   pc<-fdata2pls(fdataobj,y,ncomp=max(l),lambda=lambda,P=P,...)
}
    if (length(l)==1) l<-1:l       
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x); np <- ncol(x);lenl = length(l)
    if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    if (is.null(rownames(x)))        rownames(x) <- 1:n
    ycen = y - mean(y)
    vs <- pc$rotation$data[,,drop=FALSE]
    Z<-pc$x[,l,drop=F]
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
#            if  (pc$type=="pls") {
             if (pc$norm)  {
              sd.X <- sqrt(apply(fdataobj$data, 2, var))
              beta.est$data<-  beta.est$data/sd.X
             }      
#            }     
#    H<-diag(hat(Z, intercept = TRUE),ncol=n)
 # H2<-lm.influence(object.lm, do.coef = T)$hat# o bien
 #    I <- diag(1/(n*pc$lambdas[l]), ncol = lenl) #1/n
    Z=cbind(rep(1,len=n),Z)
     order.deriv<-0 
   if (lambda==0) mat<-0
   else {      
    if (!is.matrix(P)){
      if (is.vector(P)) {
         for (i in 1:length(P))   {      if (P[i]!=0)         order.deriv<-i}
          P<-P.penalty(tt,P)                
         P<-vs%*%P%*%t(vs) 
         P<-P*(diff(rtt)/(np -1))^(order.deriv*2-1)
         }         }
    mat<-diag(lenl+1)
    mat[-1,-1]<-lambda*P
    mat[1,1]<-0   
    }
    S<-t(Z)%*%Z+mat
    S<-Minverse(S)
    H<-Z%*%S%*%t(Z)
    e<-object.lm$residuals
    df = max(traza(H),pc$df[lenl]+1)
    rdf<-n-df
    sr2 <- sum(e^2)/rdf
    Vp<-sr2*S 
    r2 <- 1 - sum(e^2)/sum(ycen^2)
#    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
#    GCV <- sum(e^2)/rdf^2             #GCV=GCV,

    std.error = sqrt(diag(S) *sr2)
    t.value =object.lm$coefficients/std.error
    p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
    coefficients <- cbind(object.lm$coefficients, std.error, t.value, p.value)
    colnames(coefficients) <- c("Estimate", "Std. Error","t value", "Pr(>|t|)")
    
 out <- list(call = C,coefficients=object.lm$coefficients, residuals = object.lm$residuals,
 fitted.values =object.lm$fitted.values, beta.est = beta.est,coefs=coefficients,
 H=H,df = df,r2=r2, sr2 = sr2, Vp=Vp,l = l,lambda=lambda,P=P, fdata.comp=pc,
 lm=object.lm,fdataobj = fdataobj,y = y)
    class(out) = "fregre.fd"
    return(out)
}
#################################################################
#################################################################

 
  
   
   
   
#################################################################
#################################################################
fregre.pls.cv=function (fdataobj, y, kmax=8,lambda=0,P=c(0,0,1),
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
   if (is.vector(P)) {    P2<-P.penalty(tt,P=P)   }
   if (is.logical(lambda[1]))   {
      normx<-sqrt(sum(abs((t(x)%*%x)^2)))
      normp<-sqrt(sum(abs(P2^2)))
      normip<-sqrt(sum(abs(diag(nc)+P2)^2))
      pii <-traza(P2)
      lambda0<-(-2*pii+sqrt(4*pii^2-4*(nc-normx^2)*normp^2))/(2*normp^2)
      #  lambda1<-normx/normp
      lambda<-seq(0,sqrt(lambda0),len=10)
    }
  }
  tol<-sqrt(.Machine$double.eps)
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
        pls<-fdata2pls(fdataobj,y,ncomp=kmax,lambda=lambda[r],P=P,...)
        for (j in 1:kmax) {
            pls2<-pls
            pls2$rotation<-pls$rotation[1:j]
            out = fregre.pls(pls2,y,lambda=lambda[r],P=P,...)
            ck<-out$df
            s2 <- sum(out$residuals^2)/n  #(n-ck)
            cv.AIC[r,j]<-switch(criteria,
              "AIC"=log(s2) + 2 * (ck)/n,
              "AICc"=log(s2) + 2 * (ck)/(n - ck - 2),
              "SIC"=log(s2) + log(n) * ck/n,
              "SICc"=log(s2) + log(n) * ck/(n-ck-2),
              "HQIC"=log(s2) + 2*log(log(n)) * ck/n,
              "rho"={A<-out$residuals;B<-1-diag(out$H)/n; D1<-(A/B)^2;sum(D1)})
     if ( MSC.min>(cv.AIC[r,j]+tol)) {
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
            out = fregre.pls(fdataobj[-i], y[-i],lambda=lambda[r],P=P,...)
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
    rownames(cv.AIC) = paste("lambda=",signif(lambda,4) , sep = "")
#    pc2$basis<-pc$rotation[1:pc.opt]
    fregre=fregre.pls(fdataobj,y,l=1:pc.opt,lambda=rn.opt,P=P,...)
    MSC.min = cv.AIC[rn.opt,pc.opt]
    return(list("fregre.pls"=fregre,pls.opt = 1:pc.opt,lambda.opt=lambda[rn.opt],
    MSC.min = MSC.min,MSC = cv.AIC))
}
#################################################################
#################################################################

#################################################################
#################################################################
fregre.pc.cv=function (fdataobj, y, kmax=8,lambda=0,P=c(1,0,0),criteria = "SIC",weights=rep(1,len=n),...) {
  sequen=FALSE
  if (length(kmax)>1) {
    sequen=TRUE
    l<-kmax
    kmax<-max(kmax)
  }
  if (class(fdataobj)=="fdata.comp") {
    fdataobj<-fdataobj$fdataobj
    if (min(lambda)!=0 | !is.null(P)) warning("The arguments lambda and P are not used")
  }
  else {
    if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
    tt<-fdataobj[["argvals"]]
    x<-fdataobj[["data"]]
    X<-fdata.cen(x)[[1]]$data
    np<-ncol(x)
    if (min(lambda,na.rm=TRUE)>0) {
      if (is.logical(lambda[1]))   {
        if (is.vector(P)) {    P2<-P.penalty(tt,P=P)   }
        normx<-sqrt(sum(abs((t(x)%*%x)^2)))
        normp<-sqrt(sum(abs(P2^2)))
        normip<-sqrt(sum(abs(diag(np)+P2)^2))
        pii <-traza(P2)
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
  rownames(pc.opt2)<-paste("lambda=",signif(lambda,4),sep="")
  colnames(pc.opt2)<-paste("PC(",1:kmax,")",sep="")
  MSC2<-pc.opt2
  MSC.min<-Inf
  min.rn<-lambda[1]
  if (is.na(type.i))     stop("Error: incorrect criteria")
  else {
    if (type.i < 6) {
      for (r in 1:lenrn) {
        pc<-fdata2pc(fdataobj,ncomp=kmax,lambda=lambda[r],P=P)
        #       pc<-fdata2pc(fdataobj,ncomp=kmax,...)
        #       if (!is.matrix(P)) if (is.vector(P)) P<-P.penalty(tt,P)
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
          if (sequen) {
            max.c=1
            c1<-matrix(pc$l[1:k],ncol=1)
          }
          for (j in 1:max.c) {
            pc2$rotation <- pc$rotation#[c1[, j]]
            pc2$l <- pc$l[c1[, j]]
            out = fregre.pc(pc2, y,l=c1[, j],lambda=lambda[r],P=P,weights=weights,...)
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
          #                min.AIC = min(cv.AIC)
          #                pc.opt1 <- c1[, which.min(cv.AIC)]
          #                l[[k]] = pc.opt1[k]
          #                l2[[k]] = min.AIC
          #                use[pc.opt1[k]] = TRUE
          #                l[[k + 1]] = c0[use == FALSE]
          #                c1 = t(expand.grid(l))
          #                ck = nrow(c1) + 1
          #                max.c = ncol(c1)
          if (!sequen){
            min.AIC = min(cv.AIC)
            pc.opt1 <- c1[, which.min(cv.AIC)[1]]
            l[[k]] = pc.opt1[k]
            l2[[k]] = min.AIC
            use[pc.opt1[k]] = TRUE
            l[[k + 1]] = c0[use == FALSE]
            c1 = t(expand.grid(l))
            ck = nrow(c1) + 1
            max.c = ncol(c1)
          }
          else {
            pc.opt1 <- 1:k
            l[[k]] = k
            l2[[k]] =drop(cv.AIC[1])
          }
        }
        mn = which.min(l2)[1]
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
          pcl[[i]]<-fdata2pc(fdataobj[-i,],ncomp=kmax,lambda=lambda[r],P=P)
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
          if (sequen) {   max.c=1;   c1<-matrix(pc$l[1:k],ncol=1)   }
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
              out = fregre.pc(pc2,y[-i],l=c1[,j],weights=weights[-i],...) #####
              ck<-out$df
              residuals2[i] <- ((y[i] - predict(out,fdataobj[i,]))/(n-ck))^2
            }
            cv.AIC[j] <-sum(residuals2)/n
          }
          if (!sequen){
            min.AIC = min(cv.AIC)
            pc.opt1 <- c1[, which.min(cv.AIC)[1]]
            l[[k]] = pc.opt1[k]
            l2[[k]] = min.AIC
            use[pc.opt1[k]] = TRUE
            l[[k + 1]] = c0[use == FALSE]
            c1 = t(expand.grid(l))
            ck = nrow(c1) + 1
            max.c = ncol(c1)
          }
          else { pc.opt1 <- 1:k; l[[k]] = k;l2[[k]] =drop(cv.AIC[1]) }
        }
        mn = which.min(l2)[1]
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
    if (all(is.na(MSC))) stop("System is computationally singular: try with other number of basis elements")
    mn = which.min(l2)[1]
    MSC = as.numeric(l2)
    names(pc.opt3)<-paste("PC", pc.opt3, sep = "")
    rn.opt<-lambda[min.rn]
    fregre=fregre.pc(fdataobj,y,l=drop(pc.opt3),lambda=rn.opt,P=P,weights=weights,...)
    return(list("fregre.pc"=fregre,pc.opt = pc.opt3,lambda.opt=rn.opt,
                PC.order=pc.opt2,MSC.order=MSC2))
  }
}
####################################################################
####################################################################

#################################################################
fregre.pc=function (fdataobj, y, l =NULL,lambda=0,P=c(1,0,0),weights=rep(1,len=n),...){
  if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    if (is.null(l))    {
      l<-pc$l
    }
    else if (length(l)>nrow(pc$rotation)) stop("Incorrect value for  argument l")
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]                                                    
  }
  else {
    if (is.null(l)) l<- 1:3
    if (!is.fdata(fdataobj))    fdataobj = fdata(fdataobj)
    #  omit<-omit.fdata(fdataobj,y)
    #  fdataobj<-omit[[1]]
    #  y<-omit[[2]]
    tt<-fdataobj[["argvals"]]
    x<-fdataobj[["data"]]
    pc<-fdata2pc(fdataobj,ncomp=max(l),lambda=lambda,P=P)
  }  
  rtt <- fdataobj[["rangeval"]]
  names <- fdataobj[["names"]]
  n = nrow(x); np <- ncol(x);lenl = length(l)
  if (is.null(rownames(x)))        rownames(x) <- 1:n
  X <-xcen<-pc$fdataobj.cen
  if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
  C <- match.call()
  ycen = y - mean(y)
  vs<-t(pc$rotation$data[l,,drop=F])
  scores<-Z<-(pc$x[,l,drop=F])
  cnames<-colnames(pc$x)[l]
  df<-lenl+1
  J<-min(np,lenl)
  ymean<-mean(y)
  ycen<- y - ymean
  W<-diag(weights) 
  if (is.logical(lambda)) {
    #   val<-log(.25*(pc$d[1]^2),base=2)
    lambda<-.25*(pc$d[1]^2)#lambda<-c(0,2^seq(0,val,len=10))
  }
  
  if (lambda>0) {
    xmean<-pc$mean
    d<-pc$newd[l]
    D<-diag(d)
    diagJ<-diag(J)
    #    lenrn<-length(rn)
    scores<-cbind(rep(1,n),pc$x[,l])
    order.deriv<-0     
    if (!is.matrix(P)){
      if (is.vector(P)) {
        for (i in 1:length(P))   {      if (P[i]!=0)         order.deriv<-i}
        P<-P.penalty(tt,P)
        
        P<-t(vs)%*%P%*%vs 
        P<-P*(diff(rtt)/(np -1))^(order.deriv*2-1)
      }         }
    mat<-diag(J+1)
    mat[-1,-1]<-lambda*P
    mat[1,1]<-0
    Sb<-t(scores)%*%W%*%scores+mat
    # S<-solve(Sb)       
    #    S=solve(t(Z)%*%W%*%Z)    
    S<-Minverse(Sb) 
    Cinv<-S%*%t(scores)%*%W         
    coefs<-Cinv%*%y
    yp<-drop(scores%*%coefs)
    H<-scores%*%Cinv
    df<-traza(H)
    coefs<-drop(coefs)
    names(coefs)<-c("Intercept",cnames)
    beta.est<-coefs[-1]*pc$rotation[l]
    beta.est$data<-colSums(beta.est$data)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
    e<-y-yp
    rdf<-n-df
    sr2 <- sum(e^2)/ rdf
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
    #    GCV <- sum(e^2)/(n - df)^2
    object.lm = list()
    object.lm$coefficients <- coefs
    object.lm$residuals <- drop(e)
    object.lm$fitted.values <- yp
    object.lm$x<-scores 
    object.lm$y <- y
    object.lm$rank <- df
    object.lm$df.residual <-  rdf
    Z=cbind(rep(1,len=n),Z)
    colnames(Z)[1] = "(Intercept)"
    std.error = sqrt(diag(S) *sr2)
    Vp<-sr2*S 
    t.value = coefs/std.error
    p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
    coefficients <- cbind(coefs, std.error, t.value, p.value)
    colnames(coefficients) <- c("Estimate", "Std. Error",
                                "t value", "Pr(>|t|)")
    class(object.lm) <- "lm"
    out <- list(call = C, beta.est = beta.est,coefficients=coefs,
                fitted.values =yp,residuals = e,H=H,df = df,r2=r2,#GCV=GCV,
                sr2 = sr2,Vp=Vp,l = l,lambda=lambda,fdata.comp=pc,lm=object.lm,
                coefs=coefficients,fdataobj = fdataobj,y = y)
    ##################################
  }
  else {
    #print("no rn")
    response = "y"
    dataf<-data.frame(y,Z,weights)
    colnames(dataf)<-c("y",cnames,"weights")
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf,"+",cnames[i],sep="")
    object.lm = lm(formula = pf,data=data.frame(dataf),weights=weights,x=TRUE, y=TRUE)
    beta.est<-object.lm$coefficients[2:(lenl+1)]*pc$rotation[l]
    beta.est$data<-colSums(beta.est$data)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
    Z=cbind(rep(1,len=n),Z)
    #    S=solve(t(Z)%*%W%*%Z)    
    S<-t(Z)%*%W%*%Z
    S<-Minverse(S) 
    H<-Z%*%S%*%t(Z)
    e<-object.lm$residuals
    df<-traza(df)#n- object.lm$df
    sr2 <- sum(e^2)/(n - df)
    Vp<-sr2*S 
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    #     r2.adj<- 1 - (1 - r2) * ((n -    1)/(n-df))
    #     GCV <- sum(e^2)/(n - df)^2
    out <- list(call = C, coefficients=object.lm$coefficients,residuals = e,
                fitted.values =object.lm$fitted.values,weights=weights,beta.est = beta.est,
                df = df,r2=r2,sr2 = sr2,Vp=Vp,H=H, l = l,lambda=lambda,P=P,fdata.comp=pc,
                lm=object.lm,XX=Z, fdataobj = fdataobj,y = y)
  }
  class(out) = "fregre.fd"
  return(out)
}
#############################################################









