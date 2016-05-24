ncross.rq.fitXB <-
function(y, x, B=NULL, X=NULL, taus, interc=FALSE, monotone=FALSE, adj.middle=FALSE,
    ndx=10, lambda=0, deg=3, dif=3, eps=.0001, plott=0, var.pen=NULL, ...){
#Stima dei non-crossing rq, possibly monotone
#A differenza di ncross.rq.fit1() questa usa semplici B-spline con un linear inequality constraint on the
#   B-spline coefficients
#B: la base di spline, se NULL viene costruita attraverso la variabile x
#x: la variabile rispetto a cui viene costruita la base, ammesso che B sia NULL. Quando B e' fornita
#   questa viene usata per la stima, e x viene usata solo per disegnare, ammesso che plott>0
#plott {0,1,2} se 0 non disegna, se 1 aggiunge se 2 apre un nuovo device. Se x non e' fornita
#   plott viene posto a 0.
#Aggiungere una matrice di esplicative lineari?
#--------------------------------------------------------
Rho <- function(u, tau) u * (tau - (u < 0))
#--------------------------------------------------------
bspline <- function(x, ndx, xlr = NULL, knots=NULL, deg = 3, deriv = 0, outer.ok=FALSE) {
    # x: vettore di dati
    # xlr: il vettore di c(xl,xr)
    # ndx: n.intervalli in cui dividere il range
    # deg: il grado della spline
    # Restituisce ndx+deg basis functions per ndx-1 inner nodi
    #ci sono "ndx+1" nodi interni + "2*deg" nodi esterni
#    require(splines)
  if(is.null(knots)) {
    if (is.null(xlr)) {
        xl <- min(x) - 0.01 * diff(range(x))
        xr <- max(x) + 0.01 * diff(range(x))
    }
    else {
        if (length(xlr) != 2)
            stop("quando fornito, xlr deve avere due componenti")
        xl <- xlr[1]
        xr <- xlr[2]
    }
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      }
      #else {
      #if(length(knots)!=(ndx+1+2*deg)) stop("errore nel numero di nodi fornito")
      #}
    B <- splineDesign(knots, x, ord = deg + 1, derivs = rep(deriv, length(x)), outer.ok=outer.ok)
    B
}
#-------------------------------------
blockdiag <- function(...) {
  args <- list(...)
  nc <- sapply(args,ncol)
  cumnc <- cumsum(nc)
  ##  nr <- sapply(args,nrow)
  ## NR <- sum(nr)
  NC <- sum(nc)
  rowfun <- function(m,zbefore,zafter) {
    cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
          matrix(0,ncol=zafter,nrow=nrow(m)))
  }
  ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
  for (i in 2:length(args)) {
    ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
  }
  ret
}
#-------------------------------------
build.D<-function(var.pen.ok, p.ok, dif.ok, lambda.ok){
#build the penalty matrix to be appended below the design matrix
      if(is.null(var.pen.ok)){
        xx.var.pen <- rep(1,(p.ok-dif.ok))
        } else {
      f.var.pen <- function(k) eval(parse(text = var.pen.ok))
      xx.var.pen <- 1:(p.ok - dif.ok)
      xx.var.pen <- sqrt(f.var.pen(max(xx.var.pen)))
      }
      D.ok<-lambda.ok*xx.var.pen*diff(diag(p.ok), diff=dif.ok) #ridge: diag(p) #diff(diag(p), diff=1)
      D.ok
      }
#-------------------------------------
#      require(quantreg)
      taus<-sort(taus)
      n<-length(y)
      #lambda<-sqrt(lambda)
      if(is.null(B)) {
        #if(missing(x) || !is.matrix(x)) stop("'B' (list) or 'x' (matrix) have to be supplied")
        if(missing(x) || !(is.matrix(x)||is.data.frame(x))) stop("'B' (list) or 'x' (matrix) have to be supplied")
        B<-vector("list", length=ncol(x))      
        for(j in 1:ncol(x)) B[[j]]<-bspline(x[,j], ndx=ndx[j], deg=deg[j])
      }
      if(missing(x)) plott<-0
      #deve diventare una matrice..
      all.p<-sapply(B, ncol)
      pSmooth<- sum(all.p)
      H<-length(B) #no. of smooth terms
      B<-matrix(unlist(B), n, pSmooth)
      XB<-cbind(X,B)
      pLin<-ncol(XB)-pSmooth #pSmooth=p1 #n. termini lineari
      p<-ncol(XB) #p1+p2
      DD<-D1<-vector("list", length=H)
      for(j in 1:H){
          D1[[j]]<-if(monotone[j]!=0) sign(monotone[j])*diff(diag(all.p[j]), diff=1) else matrix(0,all.p[j]-1,all.p[j])
          DD[[j]]<- build.D(var.pen[j], all.p[j], dif[j], lambda[j])
          }      

      R.monot<-if(length(D1)<=1) D1[[1]] else do.call("blockdiag",D1)
      R.monot<-cbind(matrix(0,nrow=nrow(R.monot), ncol=pLin), R.monot)
      r.monot<-rep(0, sum(all.p-1))
      D.pen<-if(length(DD)<=1) DD[[1]] else do.call("blockdiag",DD)
      P<-blockdiag(diag(rep(0,pLin), ncol=pLin),D.pen)

      id.start.tau<-which.min(abs(taus-0.5))
      start.tau<-taus[id.start.tau]
      pos.taus<-taus[(taus-start.tau)>0]
      neg.taus<-taus[(taus-start.tau)<0]
      n.pos.taus<-length(pos.taus)
      n.neg.taus<-length(neg.taus)

      if(any(lambda>0)){
          XB<-rbind(XB, P)
          }
      if(interc || H>1) XB<-rbind(XB, .001*diag(ncol(XB))) #a small ridge penalty
      
      y<-c(y, rep(0,nrow(XB)-n))
      
      o.start<-if(any(monotone!=0)) rq.fit(x=XB,y=y,tau=start.tau,method="fnc",R=R.monot,r=r.monot)
       else rq.fit(x=XB,y=y,tau=start.tau)
      
#      if(monotone!=0) {
#          D1<-sign(monotone)*diff(diag(p1), diff=1)
#          R<-cbind(matrix(0,nrow(D1),p2), D1) #era cbind(rep(0,p2), D1)
#          DD<- xx.var.pen*diff(diag(p1), diff=dif) #ridge: diag(c(0,rep(1,p-1)))
#          P<-blockdiag(diag(rep(0,p2), ncol=p2),lambda*DD)
#          if(interc) P<-rbind(P, .000001*diag(ncol(XB))) #a small ridge penalty
#          XB<-rbind(XB, P)
#          y<-c(y, rep(0,nrow(P)))
#          o.start<-rq.fit(x=XB,y=y,tau=start.tau,method="fnc",R=R,r=rep(0,p1-1))
#          } else {
#          DD<-xx.var.pen*diff(diag(p1), diff=dif) #ridge: diag(p) #diff(diag(p), diff=1)
#          P<-blockdiag(diag(rep(0,p2),ncol=p2),lambda*DD)
#          if(interc) P<-rbind(P, .000001*diag(ncol(XB))) #a small ridge penalty
#          XB<-rbind(XB, P)
#          y<-c(y, rep(0,nrow(P)))
#          o.start<-rq.fit(x=XB,y=y,tau=start.tau)
#          }

    if(length(taus)<=1){ #se length(taus)==1
      all.COEF<-o.start$coef
      #colnames(all.COEF)<-paste(taus)
      all.df<- sum(round(o.start$residuals[1:n],2)==0)
      all.rho<-sum(Rho(o.start$residuals[1:n], start.tau)) 
      r<-list(coefficients=all.COEF,B=XB, df=all.df, rho=all.rho,
              fitted.values=o.start$fitted.values[1:n],residuals=o.start$residuals[1:n])
      } else { #se length(taus)>1
    Ident<-diag(p)
    COEF.POS<-COEF.NEG<-FIT.POS<-FIT.NEG<-RES.POS<-RES.NEG<-NULL
    df.pos.tau<-df.neg.tau<-rho.pos.tau<-rho.neg.tau<-NULL
    if(n.pos.taus>0){
      rho.pos.tau <-df.pos.tau <- vector(length=n.pos.taus)
      COEF.POS<-matrix(,ncol(XB),n.pos.taus)
      colnames(COEF.POS)<-paste(pos.taus)
      b.start<-o.start$coef
      if(any(monotone!=0)){
          RR<-rbind(Ident,R.monot)
          rr<-c(b.start + eps, r.monot)
          } else {
          RR<- Ident
          rr<- b.start + eps
          }
      #se global penalty
#      if(monotone!=0){
#          RR<-rbind(Ident,R)
#          rr<-c(b.start + eps, rep(0,p1-1))
#          } else {
#          RR<- Ident
#          rr<- b.start + eps
#          }
      FIT.POS<-RES.POS<-matrix(,n,n.pos.taus)
      for(i in 1:n.pos.taus){
            o<-rq.fit(x=XB,y=y,tau=pos.taus[i],method="fnc",R=RR,r=rr)
            FIT.POS[,i]<-o$fitted.values[1:n]
            RES.POS[,i]<-o$residuals[1:n]
            #estrai la f. obiettivo
            df.pos.tau[i] <- sum(abs(o$residuals[1:n])<=0.000001)
            rho.pos.tau[i] <- sum(Rho(o$residuals[1:n], pos.taus[i]))
            b.start<-o$coef
            COEF.POS[,i]<-b.start
            rr<- if(any(monotone!=0)) c(b.start+eps, r.monot) else b.start + eps
#            rr<- if(monotone!=0) c(b.start+eps, rep(0,p1-1) ) else b.start+eps
            }#end for
      }#end if(n.pos.taus>0)

    if(n.neg.taus>0){
      rho.neg.tau <-df.neg.tau <- vector(length=n.pos.taus)
      COEF.NEG<-matrix(,ncol(XB),n.neg.taus)
      colnames(COEF.NEG)<-paste(neg.taus)
      b.start<-o.start$coef
      neg.taus<-sort(neg.taus,TRUE)
#      if(monotone!=0){
#          RR<-rbind(-Ident,R)
#          rr<-c(-b.start+eps, rep(0, p1-1) )
#          } else {
#          RR<- -Ident
#          rr<- -b.start+eps
#          }
      if(any(monotone!=0)){
          Ident<-diag(p)
          RR<-rbind(-Ident,R.monot)
          rr<-c(-b.start + eps, r.monot)
          } else {
          RR<- -Ident
          rr<- -b.start + eps
          }

      FIT.NEG<-RES.NEG<-matrix(,n,n.neg.taus)
      for(i in 1:n.neg.taus){
            o<-rq.fit(x=XB,y=y,tau=neg.taus[i],method="fnc",R=RR,r=rr)
            FIT.NEG[,i]<-o$fitted.values[1:n]
            RES.NEG[,i]<-o$residuals[1:n]
            df.neg.tau[i] <- sum(abs(o$residuals[1:n])<=0.000001)
            rho.neg.tau[i] <- sum(Rho(o$residuals[1:n], neg.taus[i]))
            b.start<-o$coef
            COEF.NEG[,i]<-b.start
#            rr<- if(monotone!=0) c(-b.start+eps, rep(0,p1-1) ) else -b.start+eps
            rr<- if(any(monotone!=0)) c(-b.start+eps, r.monot) else -b.start + eps
            }#end for
      }#end if(n.neg.taus>0)
#------------------------------
#      if(adj.middle){
#            if(monotone!=0){
#              RR<-rbind(Ident,-Ident,R)
#              rr<-c(COEF.NEG[,1],-COEF.POS[,1], rep(0, p-1))
#              } else {
#              RR<-rbind(Ident, -Ident)
#              rr<-c(COEF.NEG[,1],-COEF.POS[,1])
#              }
#      o.start<-rq.fit(x=XB,y=y,tau=start.tau,method="fnc",R=RR,r=rr)
#        }
#-------------------------------
      all.COEF<-cbind(COEF.NEG[,n.neg.taus:1], o.start$coef, COEF.POS)
      colnames(all.COEF)<-paste(taus)
      all.FIT<-cbind(FIT.NEG[,n.neg.taus:1,drop=FALSE], o.start$fitted.values[1:n], FIT.POS)
      colnames(all.FIT)<-paste(taus)
      all.RES<-cbind(RES.NEG[,n.neg.taus:1,drop=FALSE], o.start$residuals[1:n], RES.POS)
      colnames(all.RES)<-paste(taus)
      all.df<- c(df.neg.tau[n.neg.taus:1], sum(abs(o.start$residuals[1:n])<=.0000001), df.pos.tau)
      all.rho<-c(rho.neg.tau[n.neg.taus:1], sum(Rho(o.start$residuals[1:n], start.tau)) , rho.pos.tau)
      r<-list(coefficients=all.COEF,B=XB, df=all.df, rho=all.rho, fitted.values=all.FIT, residuals=all.RES)
      }
      if(plott>0){
          p2=NULL
          if(plott==1) {matlines(x, B[1:n,]%*%all.COEF[-(1:p2),] ,lwd=2,...)
              } else {plot(x,y[1:n]); matpoints(x, B[1:n,]%*%all.COEF[-(1:p2),] ,lwd=2, type="l",...)}
          }
      return(r)
      }
