ncross.rq.fitX <-
function(y, X=NULL, taus, lambda=0, adj.middle=FALSE, eps=.0001, ...){
#Stima dei non-crossing rq con X lineari (la X dovrebbe avere una colonna di 1, se richiesta..)
#Se lambda>0 viene considerata una penalita' ridge (bhu?? non so se funziona..)
#--------------------------------------------------------
Rho <- function(u, tau) u * (tau - (u < 0))
#-------------------------------------
#      require(quantreg)
      if(length(taus)<=1){
        o<-rq.fit(x=X, y=y, tau=taus, ...)
        return(o)
        }
      n<-length(y)
      B<-X
      p<-ncol(B)
      Ident<-diag(p)

      taus<-sort(taus)

      id.start.tau<-which.min(abs(taus-0.5))
      start.tau<-taus[id.start.tau]

      pos.taus<-taus[(taus-start.tau)>0]
      neg.taus<-taus[(taus-start.tau)<0]
      n.pos.taus<-length(pos.taus)
      n.neg.taus<-length(neg.taus)

      if(lambda>0){
        DD<- diag(p)
        B<-rbind(B, lambda*DD)
        y<-c(y, rep(0,nrow(DD)))
        }
      o.start<-rq.fit(x=B,y=y,tau=start.tau)
      if(length(taus)<=1){ #se length(taus)==1
            all.COEF<-o.start$coef
            #colnames(all.COEF)<-paste(taus)
            all.df<- sum(round(o.start$residuals[1:n],2)==0)
            all.rho<-sum(Rho(o.start$residuals[1:n], start.tau)) 
            r<-list(coefficients=all.COEF,B=B, df=all.df, rho=all.rho,
                    fitted.values=o.start$fitted.values[1:n],residuals=o.start$residuals[1:n])
            } else { #se length(taus)>1
    COEF.POS<-COEF.NEG<-FIT.POS<-FIT.NEG<-RES.POS<-RES.NEG<-NULL
    df.pos.tau<-df.neg.tau<-rho.pos.tau<-rho.neg.tau<-NULL
    if(n.pos.taus>0){
      rho.pos.tau <-df.pos.tau <- vector(length=n.pos.taus)
      COEF.POS<-matrix(,ncol(B),n.pos.taus)
      colnames(COEF.POS)<-paste(pos.taus)
      b.start<-o.start$coef

      RR<- Ident
      rr<- b.start + eps
      FIT.POS<-RES.POS<-matrix(,n,n.pos.taus)
      for(i in 1:n.pos.taus){
            o<-rq.fit(x=B,y=y,tau=pos.taus[i],method="fnc",R=RR,r=rr)
            FIT.POS[,i]<-o$fitted.values[1:n]
            RES.POS[,i]<-o$residuals[1:n]
            #estrai la f. obiettivo
            df.pos.tau[i] <- sum(abs(o$residuals[1:n])<=.000001) #length(o$coef) 
            rho.pos.tau[i] <- sum(Rho(o$residuals[1:n], pos.taus[i]))
            b.start<-o$coef
            COEF.POS[,i]<-b.start
            rr<- b.start + eps
            }#end for
      }#end if(n.pos.taus>0)

    if(n.neg.taus>0){
      rho.neg.tau <-df.neg.tau <- vector(length=n.pos.taus)
      COEF.NEG<-matrix(,ncol(B),n.neg.taus)
      colnames(COEF.NEG)<-paste(neg.taus)
      b.start<-o.start$coef
      neg.taus<-sort(neg.taus,TRUE)
      RR<- -Ident
      rr<- -b.start + eps
      FIT.NEG<-RES.NEG<-matrix(,n,n.neg.taus)
      for(i in 1:n.neg.taus){
            o<-rq.fit(x=B,y=y,tau=neg.taus[i],method="fnc",R=RR,r=rr)
            FIT.NEG[,i]<-o$fitted.values[1:n]
            RES.NEG[,i]<-o$residuals[1:n]
            df.neg.tau[i] <- sum(abs(o$residuals[1:n])<=.000001) #length(o$coef)
            rho.neg.tau[i] <- sum(Rho(o$residuals[1:n], neg.taus[i]))
            b.start<-o$coef
            COEF.NEG[,i]<-b.start
            rr<- -b.start + eps
            }#end for
      }#end if(n.neg.taus>0)
#------------------------------
    monotone<-FALSE
    R=NULL
      if(adj.middle){
            if(monotone){
              RR<-rbind(Ident,-Ident,R)
              rr<-c(COEF.NEG[,1],-COEF.POS[,1], rep(0, p-1))
              } else {
              RR<-rbind(Ident, -Ident)
              rr<-c(COEF.NEG[,1],-COEF.POS[,1])
              }
      o.start<-rq.fit(x=B,y=y,tau=start.tau,method="fnc",R=RR,r=rr)
        }
#-------------------------------
      all.COEF<-cbind(COEF.NEG[,n.neg.taus:1,drop=FALSE], o.start$coef, COEF.POS)
      colnames(all.COEF)<-paste(taus)
      all.FIT<-cbind(FIT.NEG[,n.neg.taus:1,drop=FALSE], o.start$fitted.values[1:n], FIT.POS)
      colnames(all.FIT)<-paste(taus)
      all.RES<-cbind(RES.NEG[,n.neg.taus:1,drop=FALSE], o.start$residuals[1:n], RES.POS)
      colnames(all.RES)<-paste(taus)
      all.df<- c(df.neg.tau[n.neg.taus:1], sum(abs(o.start$residuals[1:n])<=.000001), df.pos.tau)
      all.rho<-c(rho.neg.tau[n.neg.taus:1], sum(Rho(o.start$residuals[1:n], start.tau)) , rho.pos.tau)
      r<-list(coefficients=all.COEF,B=B, df=all.df, rho=all.rho, fitted.values=all.FIT, residuals=all.RES)
      }
      return(r)
      }
