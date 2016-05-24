booteval.yhat<-function(regrOut,boot.out,bty="bca",level=.95,prec=3){

if((bty != "perc") & (bty != "norm") & (bty != "basic") & (bty != "bca")) stop ("Only norm, basic, perc, and bca interval types are supported")

l<-nrow(regrOut$PredictorMetrics)-1
m<-ncol(regrOut$PredictorMetrics)

n<-nrow(regrOut$APSRelatedMetrics)-1
o<-ncol(regrOut$APSRelatedMetrics)-1

lowerCIpm<-upperCIpm<-regrOut$PredictorMetrics[-(l+1),]
for (i in 1: l){
  for (j in 1:m){
      indpm<-(j-1)*l+i
      CI<-boot.ci(boot.out,index=indpm,type=bty,conf=level)
      CIl<-ci.yhat(bty,CI)
      lowerCIpm[i,j]<-CIl[1]
      upperCIpm[i,j]<-CIl[2]
   }
}
lowerCIpm<-round(lowerCIpm,digits=prec)
upperCIpm<-round(upperCIpm,digits=prec)
combCIpm<-combCI(lowerCIpm,upperCIpm,regrOut$PredictorMetrics[-nrow(regrOut$PredictorMetrics),])

pmc<-combn(1:l,2)
upperCIpmDiff<-lowerCIpmDiff<-pmDiff<-matrix(nrow=ncol(pmc),ncol=m)
colnames(upperCIpmDiff)<-colnames(lowerCIpmDiff)<-colnames(pmDiff)<-colnames(regrOut$PredictorMetrics)
rownames(upperCIpmDiff)<-rownames(lowerCIpmDiff)<-rownames(pmDiff)<-paste(rownames(regrOut$PredictorMetrics)[pmc[1,]], 
                                                                          rownames(regrOut$PredictorMetrics)[pmc[2,]],sep="-")

boot.outl<-boot.out                                                                       
for (i in 1:m){
  for (j in 1:ncol(pmc)){
    pmDiff[j,i]<-boot.outl$t0<-boot.out$t0[l*(i-1)+pmc[1,j]]-boot.out$t0[l*(i-1)+pmc[2,j]]
    boot.outl$t<-as.matrix(boot.out$t[,l*(i-1)+pmc[1,j]]-boot.out$t[,l*(i-1)+pmc[2,j]])
    CI<-boot.ci(boot.outl,index=1,type=bty,conf=level)
    CIl<-ci.yhat(bty,CI)
    lowerCIpmDiff[j,i]<-CIl[1]
    upperCIpmDiff[j,i]<-CIl[2]
  }
}
pmDiff<-round(pmDiff,digits=prec)
lowerCIpmDiff<-round(lowerCIpmDiff,digits=prec)
upperCIpmDiff<-round(upperCIpmDiff,digits=prec)
combCIpmDiff<-combCI(lowerCIpmDiff,upperCIpmDiff,pmDiff,TRUE)

nd<-nrow(regrOut$PairedDominanceMetrics)
domboot<-matrix(nrow=nd*3,ncol=7)
rownames(domboot)<-c(paste("Com_",rownames(regrOut$PairedDominanceMetrics),sep=""),
                     paste("Con_",rownames(regrOut$PairedDominanceMetrics),sep=""),
                     paste("Gen_",rownames(regrOut$PairedDominanceMetrics),sep=""))
colnames(domboot)<-c("Dij","Mean","SE","Pij","Pji","Pijno","Reprod")
domboot<-as.data.frame(domboot)
inds<-indpm+1
inde<-indpm+(l*(l-1)/2)*ncol(regrOut$PairedDominanceMetrics)
ind<-indpm
boots<-nrow(boot.out$t)

domboot$Dij<-boot.out$t0[inds:inde]
domboot$Mean<-colMeans(boot.out$t[,c(inds:inde)])

for (i in 1:(nd*3)){
  domboot$SE[i]<-sd(boot.out$t[,(ind+i)])
  domboot$Pij[i]<-sum(boot.out$t[,(ind+i)] == 1)/boots
  domboot$Pji[i]<-sum(boot.out$t[,(ind+i)] == 0)/boots
  domboot$Pijno[i]<-sum(boot.out$t[,(ind+i)] == .5)/boots
  if (domboot$Dij[i]==1)
    domboot$Reprod[i]<-domboot$Pij[i]
  else if (domboot$Dij[i]==0)
    domboot$Reprod[i]<-domboot$Pji[i]
  else
    domboot$Reprod[i]<-domboot$Pijno[i]
}

lowerCIaps<-upperCIaps<-regrOut$APSRelatedMetrics[-(n+1),-2]
for (i in 1: ncol(lowerCIaps)){
  for (j in 1:nrow(lowerCIaps)){
      ind<-inde+(i-1)*n+j
      if (!is.na(lowerCIaps[j,i])){
        CI<-boot.ci(boot.out,index=ind,type=bty,conf=level)
        CIl<-ci.yhat(bty,CI)
        lowerCIaps[j,i]<-CIl[1]
        upperCIaps[j,i]<-CIl[2]
      }
  }
}


lowerCIaps<-round(lowerCIaps,digits=prec)
upperCIaps<-round(upperCIaps,digits=prec)
combCIaps<-combCI(lowerCIaps,upperCIaps,regrOut$APSRelatedMetrics[-nrow(regrOut$APSRelatedMetrics),-2])

s<-1
e<-l
d<-NULL
for (i in 1:(l-1)){
  d<-rbind(d,t(combn(s:e,2)))
  s<-e+1
  e<-e+ncol(combn(1:l,(i+1)))
}

upperCIapsDiff<-lowerCIapsDiff<-apsDiff<-matrix(nrow=nrow(d),ncol=2)
colnames(upperCIapsDiff)<-colnames(lowerCIapsDiff)<-colnames(apsDiff)<-colnames(regrOut$APSRelatedMetrics)[c(1,3)]
rownames(upperCIapsDiff)<-rownames(lowerCIapsDiff)<-rownames(apsDiff)<-paste(rownames(regrOut$APSRelatedMetrics)[d[,1]], 
                                                                          rownames(regrOut$APSRelatedMetrics)[d[,2]],sep="-")

for (i in 1:2){
  for (j in 1:nrow(d)){
    apsDiff[j,i]<-boot.outl$t0<-boot.out$t0[inde+n*(i-1)+d[j,1]]-boot.out$t0[inde+n*(i-1)+d[j,2]]
    boot.outl$t<-as.matrix(boot.out$t[,inde+n*(i-1)+d[j,1]]-boot.out$t[,inde+n*(i-1)+d[j,2]])
    if (!is.na(boot.outl$t0)){
      CI<-boot.ci(boot.outl,index=1,type=bty,conf=level)
      CIl<-ci.yhat(bty,CI)
      lowerCIapsDiff[j,i]<-CIl[1]
      upperCIapsDiff[j,i]<-CIl[2]
    }
  }
}
apsDiff<-round(apsDiff,digits=prec)
lowerCIapsDiff<-round(lowerCIapsDiff,digits=prec)
upperCIapsDiff<-round(upperCIapsDiff,digits=prec)
combCIapsDiff<-combCI(lowerCIapsDiff,upperCIapsDiff,apsDiff,TRUE)

p<-nrow(regrOut$APSRelatedMetrics)-ncol(combn(1:l,l-1))-2
upperCIincDiff<-lowerCIincDiff<-incDiff<-matrix(nrow=p,ncol=nrow(combCIpmDiff))

indi<-inde+n*2
for (i in 1:ncol(incDiff)){
  for (j in 1:nrow(incDiff)){
    incDiff[j,i]<-boot.outl$t0<-boot.out$t0[indi+n*(pmc[1,i]-1)+j]-boot.out$t0[indi+n*(pmc[2,i]-1)+j]
    boot.outl$t<-as.matrix(boot.out$t[,indi+n*(pmc[1,i]-1)+j]-boot.out$t[,indi+n*(pmc[2,i]-1)+j])
    if (!is.na(boot.outl$t0)){
      CI<-boot.ci(boot.outl,index=1,type=bty,conf=level)
      CIl<-ci.yhat(bty,CI)
      lowerCIincDiff[j,i]<-CIl[1]
      upperCIincDiff[j,i]<-CIl[2]
    }
  }
}

incDiff<-round(rbind(pmDiff[,"CD:0"],incDiff),digits=prec)
lowerCIincDiff<-round(rbind(lowerCIpmDiff[,"CD:0"],lowerCIincDiff),digits=prec)
upperCIincDiff<-round(rbind(upperCIpmDiff[,"CD:0"],upperCIincDiff),digits=prec)
rownames(upperCIincDiff)<-rownames(lowerCIincDiff)<-rownames(incDiff)<-c(".",rownames(regrOut$APSRelatedMetrics)[1:p])
colnames(upperCIincDiff)<-colnames(lowerCIincDiff)<-colnames(incDiff)<-paste(colnames(regrOut$APSRelatedMetrics)[pmc[1,]+3], 
                                                                             colnames(regrOut$APSRelatedMetrics)[pmc[2,]+3],sep="-")
combCIincDiff<-combCI(lowerCIincDiff,upperCIincDiff,incDiff,TRUE)

tauds<-matrix(ncol=2,nrow=m)
colnames(tauds)<-c("Mean","SD")
rownames(tauds)<-colnames(regrOut$PredictorMetrics)

for (i in 1:m){
  tauds[i,"Mean"]<-mean(boot.out$t[,(ind+i)])
  tauds[i,"SD"]<-sd(boot.out$t[,(ind+i)])
}

return(list(combCIpm=combCIpm,lowerCIpm=lowerCIpm,upperCIpm=upperCIpm,
            combCIaps=combCIaps,lowerCIaps=lowerCIaps,upperCIaps=upperCIaps,
            domBoot=round(domboot,digits=prec),tauDS=t(round(tauds,digits=prec)),
            combCIpmDiff=combCIpmDiff,lowerCIpmDiff=lowerCIpmDiff,upperCIpmDiff=upperCIpmDiff,
            combCIapsDiff=combCIapsDiff,lowerCIapsDiff=lowerCIapsDiff,upperCIapsDiff=upperCIapsDiff,
            combCIincDiff=combCIincDiff,lowerCIincDiff=lowerCIincDiff,upperCIincDiff=upperCIincDiff))
}
