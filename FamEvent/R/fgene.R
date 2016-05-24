fgene <-
function(base.dist, affage, affsex, variation="secondgene", parms, vbeta, alpha,
          pg=0, m.carrier=0, dominant.m=TRUE, aq){
# returns 2x3 matrix , first row for the major gene, 
#	second row for the second gene
pAA <- pAa <- paa <- 0

AAq <- Aaq <- aaq <- 0
AAq[1] <- Pgene(1, pg=pg[1], a.freq=aq[1])
Aaq[1] <- Pgene(2, pg=pg[1], a.freq=aq[1])
aaq[1] <- Pgene(3, pg=pg[1], a.freq=aq[1])

if(variation=="secondgene"){
  if(pg[1]==0) pg[2]<-0
  AAq[2] <- Pgene(1, pg=pg[2], a.freq=aq[2])
  Aaq[2] <- Pgene(2, pg=pg[2], a.freq=aq[2])
  aaq[2] <- Pgene(3, pg=pg[2], a.freq=aq[2])
  Ft <- matrix(0, ncol=2, nrow=2)
  for(i in c(0,1)) for(j in c(0,1)) {
	xbeta <- affsex*vbeta[1] + i*vbeta[2] + j*vbeta[3]
  	Ft[(i+1),(j+1)] <- 1- surv.dist(base.dist, affage, parms, xbeta, alpha[2], res=0)
                    }
  if(!dominant.m){
    pAA[1] <- Ft[2,2]*AAq[1]*(AAq[2]+Aaq[2]) + Ft[2, 1]*AAq[1]*aaq[2]
    pAa[1] <- Ft[1,2]*Aaq[1]*(AAq[2]+Aaq[2]) + Ft[1, 1]*Aaq[1]*aaq[2]
    paa[1] <- Ft[1,2]*aaq[1]*(AAq[2]+Aaq[2]) + Ft[1, 1]*aaq[1]*aaq[2]
    pAA[2] <- Ft[2,2]*AAq[2]*AAq[1]+Ft[1,2]*AAq[2]*(Aaq[1]+aaq[1])*(1-m.carrier)
    pAa[2] <- Ft[2,2]*Aaq[2]*AAq[1]+Ft[1,2]*Aaq[2]*(Aaq[1]+aaq[1])*(1-m.carrier)
    paa[2] <- Ft[2,1]*aaq[2]*AAq[1]+Ft[1,1]*aaq[2]*(Aaq[1]+aaq[1])*(1-m.carrier)
    }
 else{ 
    pAA[1] <- Ft[2,2]*AAq[1]*(AAq[2]+Aaq[2]) + Ft[2,1]*AAq[1]*aaq[2]
    pAa[1] <- Ft[2,2]*Aaq[1]*(AAq[2]+Aaq[2]) + Ft[2,1]*Aaq[1]*aaq[2]
    paa[1] <- Ft[1,2]*aaq[1]*(AAq[2]+Aaq[2]) + Ft[1,1]*aaq[1]*aaq[2]
    pAA[2] <- Ft[2,2]*AAq[2]*(AAq[1]+Aaq[1]) + Ft[1,2]*AAq[2]*aaq[1]*(1-m.carrier)
    pAa[2] <- Ft[2,2]*Aaq[2]*(AAq[1]+Aaq[1]) + Ft[1,2]*Aaq[2]*aaq[1]*(1-m.carrier)
    paa[2] <- Ft[2,1]*aaq[2]*(AAq[1]+Aaq[1]) + Ft[1,1]*aaq[2]*aaq[1]*(1-m.carrier)
    }
}# close if variation=="secondgene"
else{
  Ft <- 0
  for(i in c(0,1)){
  	xbeta <- affsex*vbeta[1] + i*vbeta[2]
  	Ft[i+1] <- 1 - surv.dist(base.dist, affage, parms, xbeta, alpha[2], res=0)
  } 

pAA <- Ft[2]*AAq[1]
if(dominant.m) pAa <- Ft[2]*Aaq[1]
else pAa <- Ft[1]*Aaq[1]
paa <- Ft[1]*aaq[1]
}
return( cbind(pAA, pAa, paa))
}
