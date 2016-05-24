"pedStatSummary"<-function(pedStats,pedStats2=NULL){

  sumData<-array(NA,dim=22)

  names(sumData)<-c("records ","maternities","paternities","full sibs","maternal sibs","maternal half sibs","paternal sibs","paternal half sibs","maternal grandmothers","maternal grandfathers","paternal grandmothers","paternal grandfathers","maximum pedigree depth","founders","mean maternal sibsip size","mean paternal sibsip size","non-zero F","F > 0.125","mean pairwise relatedness","pairwise relatedness>=0.125","pairwise relatedness>=0.25","pairwise relatedness>=0.5")

  sumData[1]<-pedStats$totalSampleSize
  sumData[2]<-pedStats$totalMaternities
  sumData[3]<-pedStats$totalPaternities
  sumData[4]<-pedStats$totalFullSibs
  sumData[5]<-pedStats$totalMaternalSibs
  sumData[6]<-pedStats$totalMaternalSibs-pedStats$totalFullSibs
  sumData[7]<-pedStats$totalPaternalSibs
  sumData[8]<-pedStats$totalPaternalSibs-pedStats$totalFullSibs
  sumData[9]<-pedStats$totalMaternalGrandmothers
  sumData[10]<-pedStats$totalMaternalGrandfathers
  sumData[11]<-pedStats$totalPaternalGrandmothers
  sumData[12]<-pedStats$totalPaternalGrandfathers
  sumData[13]<-max(as.numeric(names(pedStats$pedigreeDepth)))
  sumData[14]<-pedStats$pedigreeDepth[1]
  sumData[15]<-mean(pedStats$maternalSibships[,2])
  sumData[16]<-mean(pedStats$paternalSibships[,2])
  sumData[17]<-sum(pedStats$inbreedingCoefficients!=0+0)
  sumData[18]<-sum(pedStats$inbreedingCoefficients>0.125+0)
  rc<-subset(pedStats$relatednessCategories,is.na(pedStats$relatednessCategories)==FALSE)
  sumData[19]<-weighted.mean(as.numeric(names(rc)),rc,na.rm=TRUE)
  sumData[20]<-sum(subset(rc,as.numeric(names(rc))>=0.125))/sum(rc)
  sumData[21]<-sum(subset(rc,as.numeric(names(rc))>=0.25))/sum(rc)
  sumData[22]<-sum(subset(rc,as.numeric(names(rc))>=0.5))/sum(rc)

  if(is.null(pedStats2)==FALSE) {
      sumData<-as.data.frame(cbind(sumData,NA))
      names(sumData)<-c("Pedigree 1","Pedigree 2")
    sumData[1,2]<-pedStats2$totalSampleSize
    sumData[2,2]<-pedStats2$totalMaternities
    sumData[3,2]<-pedStats2$totalPaternities
    sumData[4,2]<-pedStats2$totalFullSibs
    sumData[5,2]<-pedStats2$totalMaternalSibs
    sumData[6,2]<-pedStats2$totalMaternalSibs-pedStats2$totalFullSibs
    sumData[7,2]<-pedStats2$totalPaternalSibs
    sumData[8,2]<-pedStats2$totalPaternalSibs-pedStats2$totalFullSibs
    sumData[9,2]<-pedStats2$totalMaternalGrandmothers
    sumData[10,2]<-pedStats2$totalMaternalGrandfathers
    sumData[11,2]<-pedStats2$totalPaternalGrandmothers
    sumData[12,2]<-pedStats2$totalPaternalGrandfathers
    sumData[13,2]<-max(as.numeric(names(pedStats2$pedigreeDepth)))
    sumData[14,2]<-pedStats2$pedigreeDepth[1]
    sumData[15,2]<-mean(pedStats2$maternalSibships[,2])
    sumData[16,2]<-mean(pedStats2$paternalSibships[,2])
    sumData[17,2]<-sum(pedStats2$inbreedingCoefficients!=0+0)
    sumData[18,2]<-sum(pedStats2$inbreedingCoefficients>0.125+0)
    rc<-subset(pedStats2$relatednessCategories,is.na(pedStats2$relatednessCategories)==FALSE)
    sumData[19,2]<-weighted.mean(as.numeric(names(rc)),rc,na.rm=TRUE)
    sumData[20,2]<-sum(subset(rc,as.numeric(names(rc))>=0.125))/sum(rc)
    sumData[21,2]<-sum(subset(rc,as.numeric(names(rc))>=0.25))/sum(rc)
    sumData[22,2]<-sum(subset(rc,as.numeric(names(rc))>=0.5))/sum(rc)
  }
  sumData
}

