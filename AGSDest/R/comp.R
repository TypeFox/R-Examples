`comp.sT` <-
function(pT,iD){
             sT=list(K=K<-length(pT$t)-iD$T,
                  ce0=cerr(h=0,pT=pT,iD=iD),
                  a=pT$a[1:K],
                  b=pbounds(h=0,pT,iD),
                  t=pT$t[(iD$T+1):length(pT$t)]*pT$Imax-pT$t[iD$T]*pT$Imax)
             class(sT)<-"GSD"
             return(sT)
           }

`comp.z1` <-
function(pT,iD,sT,sTo){
                 (iD$z*sqrt(pT$t[iD$T]*pT$Imax)+sTo$z*sqrt(sT$t[sTo$T]*sT$Imax))/sqrt(pT$t[iD$T]*pT$Imax+sT$t[sTo$T]*sT$Imax)
          }

`comp.z2` <-
function(pT,pTo,iD){
               K<-length(pT$t);
               if(iD$T>=pTo$T) { print("Error: Cannot compute z2 for pTo$T<=iD$T") ; rep(NaN,K)
               } else {
                 (pTo$z*sqrt(pT$t[pTo$T]*pT$Imax)-ifelse(iD$T==0,0,iD$z*sqrt(pT$t[iD$T]*pT$Imax)))/
                 sqrt(pT$t[pTo$T]*pT$Imax-ifelse(iD$T==0,0,pT$t[iD$T]*pT$Imax))
                }
          }
