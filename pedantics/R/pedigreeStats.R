`pedigreeStats` <-
function(Ped,cohorts=NULL,dat=NULL,retain='informative',graphicalReport='y',includeA=TRUE,lowMem=FALSE,grContrast=FALSE) {


  "prune"<-function(pedigree, keep, make.base=FALSE){
     ind.keep<-keep
     nind<-length(ind.keep)+1
     while(length(ind.keep)!=nind){
       nind<-length(ind.keep)
       ind.keep<-union(na.omit(unlist(pedigree[,2:3][match(ind.keep,pedigree[,1]),])), ind.keep)
     }
     pedigree<-pedigree[sort(match(ind.keep, pedigree[,1])),]
     if(make.base){
       if(any(match(pedigree[,2], pedigree[,1])>match(pedigree[,1], pedigree[,1]), na.rm=T)){
                 stop("dams appearing before their offspring: use fixPedigree")}
       if(any(match(pedigree[,3], pedigree[,1])>match(pedigree[,1], pedigree[,1]), na.rm=T)){
                 stop("sires appearing before their offspring: use fixPedigree")}
       phenotyped<-pedigree[,1]%in%keep
       delete<-rep(FALSE, dim(pedigree)[1])
       for(i in 1:dim(pedigree)[1]){
         nlinks<-phenotyped[i]+sum(pedigree[,2]%in%pedigree[,1][i])+sum(pedigree[,3]%in%pedigree[,1][i])+sum(is.na(pedigree[i,][2:3])==FALSE)
         if(nlinks<2 & phenotyped[i]==FALSE){
           pedigree[,2][which(pedigree[,2]==pedigree[,1][i])]<-NA
           pedigree[,3][which(pedigree[,3]==pedigree[,1][i])]<-NA
           delete[i]<-TRUE
         }
       }
       if(any(delete)){
         pedigree<-pedigree[-which(delete),]
       }
     }
     pedigree
  }

  names(Ped)[1]<-"id"
  if(names(Ped)[2]!="dam"|names(Ped)[3]!="sire"){
    if(names(Ped)[3]%in%c("mum","mom","mother","Mum","Mmom","Dam","Mother","MUM","MOM","DAM","MOTHER")){
      cat(paste("'mum' appears to be in third column, reordering to 'id','dam','sire'")); flush.console();
      Ped<-Ped[,c(1,3,2)]
    }
    if(names(Ped)[2]%in%c("mum","mom","mother","Mum","Mom","Dam","Mother","MUM","MOM","DAM","MOTHER")){
      names(Ped)[2]<-"dam"
      names(Ped)[3]<-"sire"
    }
    if(names(Ped)[2]!="dam"|names(Ped)[3]!="sire"){
     stop("Unable to identify column names, expecting 'id','dam','sire', or similar")
    }
  }

  for(x in 1:3) Ped[,x]<-as.character(Ped[,x])
  if(is.null(dat)==FALSE&&is.numeric(dat)==FALSE) dat<-as.character(dat)

  if(is.null(cohorts)==FALSE&&length(Ped[,1])!=length(cohorts)) 
    stop("Pedigree and cohorts differ in length.")

  if(is.null(dat)==FALSE&&is.numeric(dat)&&length(Ped[,1])!=dim(as.data.frame(dat))[1]) 
    stop("Pedigree and available data differ in length.")


  getSibNums <- function(Ped)
    .C("CALC_SIB_NUMBERS",
     as.integer(length(Ped$id)),
     as.integer(Ped$id),
     as.integer(Ped$dam),
     as.integer(Ped$sire),
     as.integer(0),
     as.integer(0),
     as.integer(0),PACKAGE="pedantics")


  if(is.null(dat)==FALSE) {
    if(is.numeric(dat)){
      avail<-rowSums(as.data.frame(dat))
      keep<-subset(Ped$id,avail>0)
    } else {
      keep<-dat
    }

    if(is.null(cohorts)==FALSE) c<-cbind(Ped$id,cohorts)

    if(retain=='informative'){
      Ped<-prune(Ped,keep,make.base=TRUE)
    }else{
      Ped<-prune(Ped,keep,make.base=FALSE)
    }

    if(is.null(cohorts)==FALSE) c<-subset(c,c[,1]%in%Ped[,1])
    if(is.null(cohorts)==FALSE) cohorts<-c[,2]
  }

  # sample size

  totalSampleSize<-length(Ped$id)

  # maternities and paternities

  totalMaternities<-sum(table(Ped$dam))
  totalPaternities<-sum(table(Ped$sire))

  # siblings

  numPed<-makePedigreeNumeric(as.character(Ped$id),as.character(Ped$sire),
                           as.character(Ped$dam),missingVal=-1)

  sibNums<-getSibNums(numPed$numericPedigree)

  totalFullSibs<-sibNums[5]
  totalMaternalSibs<-sibNums[6]
  totalPaternalSibs<-sibNums[7]

  # grandparents

  grandparentData<-Ped
  grandparentData$maternalGM<-grandparentData$dam[match(grandparentData$dam,grandparentData$id)]
  grandparentData$maternalGF<-grandparentData$sire[match(grandparentData$dam,grandparentData$id)]
  grandparentData$paternalGM<-grandparentData$dam[match(grandparentData$sire,grandparentData$id)]
  grandparentData$paternalGF<-grandparentData$sire[match(grandparentData$sire,grandparentData$id)]

  totalMaternalGM<-sum(table(grandparentData$maternalGM))
  totalMaternalGF<-sum(table(grandparentData$maternalGF))
  totalPaternalGM<-sum(table(grandparentData$paternalGM))
  totalPaternalGF<-sum(table(grandparentData$paternalGF))

  # pedigree depth

  pedigreeDepth<-table(kindepth(Ped[,1],Ped[,2],Ped[,3]))

  # inbreeding coefficients

  orderPed<-function(ped){
    reorder<-ped[order(kindepth(ped[,1],ped[,2],ped[,3]), decreasing=FALSE),]
    return(reorder)
  }
  orderedPed<-as.data.frame(orderPed(Ped))
  orderedPed$inbreeding<-inverseA(orderPed(Ped))$inbreeding
  reorderInbreeding<-as.data.frame(Ped$id)
  reorderInbreeding$inbreeding<-orderedPed$inbreeding[match(reorderInbreeding[,1],orderedPed$id)]

  # sibship sizes

  matSibships<-as.data.frame(table(as.character(Ped$dam)))
  patSibships<-as.data.frame(table(as.character(Ped$sire)))

cumulativeRelatedness<-NULL
pairwiseRelatedness<-NULL
relatednessBin<-NULL
if(lowMem==FALSE){

  # relatedness classes
  cutoffs<-seq(-0.0125,0.9875,by=0.025)
  midBins<-seq(0,0.975,by=0.025)
  cumulativeRelatedness<-array(dim=length(midBins))
  names(cumulativeRelatedness)<-midBins
  relatednessBin<-array(dim=length(midBins))
  names(relatednessBin)<-midBins
  A<-kinship(Ped[,1],Ped[,3],Ped[,2])*2
  pairwiseRelatedness<-A*abs(diag(length(A[,1]))-1)
  pairwiseRelatedness[upper.tri(pairwiseRelatedness)]<-0
  pairwiseRelatedness<-c(pairwiseRelatedness)
#  pairwiseRelatedness<-subset(pairwiseRelatedness,pairwiseRelatedness>1e-9)
  for(x in 2:(length(cutoffs)-1)) {
  #  cumulativeRelatedness[x]<-table(pairwiseRelatedness<cutoffs[x+1])["TRUE"]
    relatednessBin[x]<-table(pairwiseRelatedness>cutoffs[x]&pairwiseRelatedness<=cutoffs[x+1])["TRUE"]
  }
  relatednessBin[1]<-((totalSampleSize^2-totalSampleSize)/2)-sum(relatednessBin,na.rm=TRUE)
  relatednessBin[relatednessBin %in% NA] <- 0

  rb<-relatednessBin/sum(relatednessBin)
  for(x in 1:(length(cutoffs)-1)) {
    cumulativeRelatedness[x]<-sum(rb[1:x])
  }  

}

  # MacCluer's pedigree completeness statistics

  missingness<-as.data.frame(cbind(as.character(Ped$id),matrix(0,length(Ped[,1]),max(kindepth(Ped[,1],Ped[,2],Ped[,3])))))
  p<-Ped
  for(x in 1:(length(missingness[1,])-1)){
    pOld<-p
    if(is.null(pOld[,2])==FALSE & is.null(pOld[,3])==FALSE ){
      matsub<-subset(p,is.na(p[,2])==FALSE)[,c(1,2)]
      names(matsub)<-c("id","p")
      patsub<-subset(p,is.na(p[,3])==FALSE)[,c(1,3)]
      names(patsub)<-c("id","p")
      p<-as.data.frame(rbind(matsub,patsub))
    }
    if(is.null(pOld[,2])==FALSE & is.null(pOld[,3])==TRUE){
      p<-as.data.frame(subset(p,is.na(p[,2])==FALSE)[,c(1,2)])
    }
    if(is.null(pOld[,2])==TRUE & is.null(pOld[,3])==FALSE ){
      p<-as.data.frame(subset(p,is.na(p[,3])==FALSE)[,c(1,3)])
    }
    t<-table(p[,1])
    missingness[,1+x]<-t[match(missingness[,1],names(t))]
    p$gm<-Ped[match(p[,2],Ped[,1]),2]
    p$gf<-Ped[match(p[,2],Ped[,1]),3]
    p<-p[,c(1,3,4)]
  }
  missingness[is.na(missingness)]<-0
  for(x in 1:(length(missingness[1,])-1)) missingness[,1+x]<-missingness[,1+x]/(2^x)


  ##### if cohort designations are available

  if(is.null(cohorts)==FALSE) {


    # average relatedness within and between cohorts
    vector.to.design.matrix<-function(f){
      factors<-as.numeric(as.factor(f))
      design.matrix<-matrix(0,length(factors),nlevels(as.factor(f)))
      for(x in 1:length(factors)) design.matrix[x,factors[x]]<-1
      results<-list(design.matrix,levels(as.factor(f)))
      names(results)<-c("DesignMatrix","FactorLevels")
      results
    }
    Y<-vector.to.design.matrix(cohorts)$DesignMatrix
    D<-diag(1/table(cohorts))

meanRelatednessAmongCohorts<-NULL
if(lowMem==FALSE){
    meanRelatednessAmongCohorts<-t(Y%*%D)%*%A%*%(Y%*%D)
    rownames(meanRelatednessAmongCohorts)<-names(table(cohorts))
    colnames(meanRelatednessAmongCohorts)<-names(table(cohorts))
}

    # sample size by cohort

    cohortSampleSizes<-table(cohorts)

    # maternities and paternities by cohort

    cohortMaternities<-array(0,dim=length(table(cohorts)))
    cohortPaternities<-array(0,dim=length(table(cohorts)))
    names(cohortMaternities)<-names(table(cohorts))
    names(cohortPaternities)<-names(table(cohorts))

    cohortFullSibs<-array(0,dim=length(table(cohorts)))
    cohortMaternalSibs<-array(0,dim=length(table(cohorts)))
    cohortPaternalSibs<-array(0,dim=length(table(cohorts)))

    cohortMaternalGM<-array(0,dim=length(table(cohorts)))
    cohortMaternalGF<-array(0,dim=length(table(cohorts)))
    cohortPaternalGM<-array(0,dim=length(table(cohorts)))
    cohortPaternalGF<-array(0,dim=length(table(cohorts)))
    names(cohortMaternalGM)<-names(table(cohorts))
    names(cohortMaternalGF)<-names(table(cohorts))
    names(cohortPaternalGM)<-names(table(cohorts))
    names(cohortPaternalGF)<-names(table(cohorts))

    cohortPedgireeDepth<-matrix(0,length(table(cohorts)),max(as.numeric(names(pedigreeDepth)))+1)
    rownames(cohortPedgireeDepth)<-names(table(cohorts))
    colnames(cohortPedgireeDepth)<-as.character(0:max(as.numeric(names(pedigreeDepth))))

    for(x in 1:length(cohortMaternities)){
      temp<-subset(Ped,as.character(cohorts)==names(table(cohorts))[x])
      cohortMaternities[x]<-sum(table(temp$dam))
      cohortPaternities[x]<-sum(table(temp$sire))

      numPed<-makePedigreeNumeric(as.character(temp$id),as.character(temp$sire),
                           as.character(temp$dam),missingVal=-1)
      sibNums<-getSibNums(numPed$numericPedigree)
      cohortFullSibs[x]<-sibNums[5]
      cohortMaternalSibs[x]<-sibNums[6]
      cohortPaternalSibs[x]<-sibNums[7]

      temp<-subset(Ped,as.numeric(as.character(cohorts))<=as.numeric(names(table(cohorts))[x]))
      pedDepth<-table(kindepth(temp[,1],temp[,2],temp[,3]))
      for(y in 1:length(pedDepth)) cohortPedgireeDepth[x,names(pedDepth[y])]<-pedDepth[y]

      temp<-subset(grandparentData,as.character(cohorts)==names(table(cohorts))[x])
      cohortMaternalGM[x]<-sum(table(temp$maternalGM))
      cohortMaternalGF[x]<-sum(table(temp$maternalGF))
      cohortPaternalGM[x]<-sum(table(temp$paternalGM))
      cohortPaternalGF[x]<-sum(table(temp$paternalGF))
    }
    cohortFullSibs<-as.data.frame(cohortFullSibs)
    cohortMaternalSibs<-as.data.frame(cohortMaternalSibs)
    cohortPaternalSibs<-as.data.frame(cohortPaternalSibs)
    names(cohortFullSibs)<-names(table(cohorts))
    names(cohortMaternalSibs)<-names(table(cohorts))
    names(cohortPaternalSibs)<-names(table(cohorts))


    matSibships$mumBirthYear<-cohorts[match(matSibships[,1],Ped$id)]
    patSibships$mumBirthYear<-cohorts[match(patSibships[,1],Ped$id)]


    
    # MacCluer's pedigree completeness statistics, averaged by cohort

    cohortPedigreeCompleteness<-matrix(NA,length(table(cohorts)),length(missingness[1,])-1)
    for(x in 1:length(table(cohorts))){
      subsetMissingness<-subset(missingness,as.character(cohorts)==names(table(cohorts))[x])
      for(y in 1:length(cohortPedigreeCompleteness[1,])){
        cohortPedigreeCompleteness[x,y]<-mean(subsetMissingness[,y+1])
      }
    }


  }

if(lowMem==FALSE&includeA==TRUE){

  results<-list(totalSampleSize=totalSampleSize,
                totalMaternities=totalMaternities,
                totalPaternities=totalPaternities,
                totalFullSibs=totalFullSibs[[1]],
                totalMaternalSibs=totalMaternalSibs[[1]],
                totalPaternalSibs=totalPaternalSibs[[1]],
                totalMaternalGrandmothers=totalMaternalGM,
                totalMaternalGrandfathers=totalMaternalGF,
                totalPaternalGrandmothers=totalPaternalGM,
                totalPaternalGrandfathers=totalPaternalGF,
                pedigreeDepth=pedigreeDepth,
                inbreedingCoefficients=reorderInbreeding$inbreeding,
                Amatrix=A,
                maternalSibships=matSibships,
                paternalSibships=patSibships,
                cumulativeRelatedness=cumulativeRelatedness,
                relatednessCategories=relatednessBin,
                analyzedPedigree=Ped,
                missingness)

  if(is.null(cohorts)==FALSE) {
    results<-list(totalSampleSize=totalSampleSize,
                totalMaternities=totalMaternities,
                totalPaternities=totalPaternities,
                totalFullSibs=totalFullSibs[[1]],
                totalMaternalSibs=totalMaternalSibs[[1]],
                totalPaternalSibs=totalPaternalSibs[[1]],
                totalMaternalGrandmothers=totalMaternalGM,
                totalMaternalGrandfathers=totalMaternalGF,
                totalPaternalGrandmothers=totalPaternalGM,
                totalPaternalGrandfathers=totalPaternalGF,
                sampleSizesByCohort=cohortSampleSizes,
                maternitiesByCohort=cohortMaternities,
                paternitiesByCohort=cohortPaternities,
                fullSibsByCohort=cohortFullSibs,
                maternalSibsByCohort=cohortMaternalSibs,
                paternalSibsByCohort=cohortPaternalSibs,
                maternalGrandmothersByCohort=cohortMaternalGM,
                maternalGrandfathersByCohort=cohortMaternalGF,
                paternalGrandmothersByCohort=cohortPaternalGM,
                paternalGrandfathersByCohort=cohortPaternalGF,
                pedigreeDepth=pedigreeDepth,
                cumulativePedigreeDepth=cohortPedgireeDepth,
                meanRelatednessAmongCohorts=meanRelatednessAmongCohorts,
                inbreedingCoefficients=reorderInbreeding$inbreeding,
                Amatrix=A,
                maternalSibships=matSibships,
                paternalSibships=patSibships,
                cumulativeRelatedness=cumulativeRelatedness,
                relatednessCategories=relatednessBin,
                cohorts=cohorts,
                analyzedPedigree=Ped)
  }


}else{

  results<-list(totalSampleSize=totalSampleSize,
                totalMaternities=totalMaternities,
                totalPaternities=totalPaternities,
                totalFullSibs=totalFullSibs[[1]],
                totalMaternalSibs=totalMaternalSibs[[1]],
                totalPaternalSibs=totalPaternalSibs[[1]],
                totalMaternalGrandmothers=totalMaternalGM,
                totalMaternalGrandfathers=totalMaternalGF,
                totalPaternalGrandmothers=totalPaternalGM,
                totalPaternalGrandfathers=totalPaternalGF,
                pedigreeDepth=pedigreeDepth,
                inbreedingCoefficients=reorderInbreeding$inbreeding,
                maternalSibships=matSibships,
                paternalSibships=patSibships,
                cumulativeRelatedness=cumulativeRelatedness,
                relatednessCategories=relatednessBin,
                analyzedPedigree=Ped)

  if(is.null(cohorts)==FALSE) {
    results<-list(totalSampleSize=totalSampleSize,
                totalMaternities=totalMaternities,
                totalPaternities=totalPaternities,
                totalFullSibs=totalFullSibs[[1]],
                totalMaternalSibs=totalMaternalSibs[[1]],
                totalPaternalSibs=totalPaternalSibs[[1]],
                totalMaternalGrandmothers=totalMaternalGM,
                totalMaternalGrandfathers=totalMaternalGF,
                totalPaternalGrandmothers=totalPaternalGM,
                totalPaternalGrandfathers=totalPaternalGF,
                sampleSizesByCohort=cohortSampleSizes,
                maternitiesByCohort=cohortMaternities,
                paternitiesByCohort=cohortPaternities,
                fullSibsByCohort=cohortFullSibs,
                maternalSibsByCohort=cohortMaternalSibs,
                paternalSibsByCohort=cohortPaternalSibs,
                maternalGrandmothersByCohort=cohortMaternalGM,
                maternalGrandfathersByCohort=cohortMaternalGF,
                paternalGrandmothersByCohort=cohortPaternalGM,
                paternalGrandfathersByCohort=cohortPaternalGF,
                pedigreeDepth=pedigreeDepth,
                cumulativePedigreeDepth=cohortPedgireeDepth,
                meanRelatednessAmongCohorts=meanRelatednessAmongCohorts,
                inbreedingCoefficients=reorderInbreeding$inbreeding,
                maternalSibships=matSibships,
                paternalSibships=patSibships,
                cumulativeRelatedness=cumulativeRelatedness,
                relatednessCategories=relatednessBin,
                cohorts=cohorts,
                analyzedPedigree=Ped)

  }


}

if(graphicalReport=='y'){

  col1<-'red';  col2<-'blue';
  if(grContrast==TRUE) { col1<-colors()[117]; col2<-colors()[109]; }

  if(is.null(cohorts)==FALSE&lowMem==FALSE) {
    cohortRelatedness<-as.data.frame(results$meanRelatednessAmongCohorts)
    cohortTakeOneRelatedness<-array(dim=length(cohortRelatedness[1,]))
    for(x in 1:(length(cohortTakeOneRelatedness)-1)) cohortTakeOneRelatedness[x+1]<-cohortRelatedness[x,x+1]
    cohortTakeTwoRelatedness<-array(dim=length(cohortRelatedness[1,]))
    for(x in 1:(length(cohortTakeOneRelatedness)-2)) cohortTakeTwoRelatedness[x+2]<-cohortRelatedness[x,x+2]
    cohortTakeThreeRelatedness<-array(dim=length(cohortRelatedness[1,]))
    for(x in 1:(length(cohortTakeOneRelatedness)-3)) cohortTakeThreeRelatedness[x+3]<-cohortRelatedness[x,x+3]
    cohortTakeFourRelatedness<-array(dim=length(cohortRelatedness[1,]))
    for(x in 1:(length(cohortTakeOneRelatedness)-4)) cohortTakeFourRelatedness[x+4]<-cohortRelatedness[x,x+4]
    par(oma=c(5,1,1,1))

    plot(as.numeric(names(cohortRelatedness)),cohortTakeOneRelatedness,type='l',xlab="Cohort",ylab="Pairwise mean cohort relatedness")
    lines(as.numeric(names(cohortRelatedness)),cohortTakeTwoRelatedness,lty='dashed')
    lines(as.numeric(names(cohortRelatedness)),cohortTakeThreeRelatedness,col='gray')
    lines(as.numeric(names(cohortRelatedness)),cohortTakeFourRelatedness,col='gray',lty='dashed')
    mtext("Mean relatedness between individuals born 1",side=1,line=5)
    mtext("(black, solid), 2 (black, dashed), 3 (gray,",side=1,line=6)
    mtext("solid) and 4 (gray, dashed) cohorts apart.",side=1,line=7)
    inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
    if(inp=='s'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      savePlot(paste(s,".jpeg",sep=""),type="jpeg")
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }
  }

#  inbredSubset<-as.data.frame(cbind(results$cohorts,results$inbreedingCoefficients))
  inbredSubset<-subset(results$inbreedingCoefficients, results$inbreedingCoefficients>0)
  proportionInbred<-length(inbredSubset)/results$totalSampleSize
  hist(as.numeric(as.character(inbredSubset)), xlab="Inbreeding coefficient",ylab="Count",main="")
  mtext("Distribution of inbreeding coefficients, as ",side=1,line=5)
  mtext("estimated from the pedigree, among the",side=1,line=6)
  mtext("ndividuals in the pedigree with F>0.",side=1,line=7)
  inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
  if(inp=='s'){
    s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
    savePlot(paste(s,".jpeg",sep=""),type="jpeg")
    readline(prompt = "File saved.  Press <Enter> to continue...")
  }


  if(is.null(cohorts)==FALSE) {
    plot(as.numeric(as.character(results$cohorts)),as.numeric(as.character(results$inbreedingCoefficients)), xlab="Inbreeding coefficient",ylab="Count (note: zeros excluded)")
    mtext("Distribution of inbreeding coefficients, as ",side=1,line=5)
    mtext("estimated from the pedigree, by cohort.",side=1,line=6)

    inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
    if(inp=='s'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      savePlot(paste(s,".jpeg",sep=""),type="jpeg")
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }
  }

  plot(as.numeric(names(results$pedigreeDepth)),as.numeric(as.character(results$pedigreeDepth)),type='l', xlab="Pedigree depth",ylab="Count")
  mtext("Distribution of pedigree depth.",side=1,line=5)
  inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
  if(inp=='s'){
    s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
    savePlot(paste(s,".jpeg",sep=""),type="jpeg")
    readline(prompt = "File saved.  Press <Enter> to continue...")
  }

  if(is.null(cohorts)==FALSE) {
    plot(as.numeric(colnames(results$cumulativePedigreeDepth)),
             as.numeric(as.character(results$cumulativePedigreeDepth[1,])),type='l',
               ylim=c(0,max(results$cumulativePedigreeDepth)), xlab="Pedigree depth",ylab="Count")
    for(x in 2:length(results$cumulativePedigreeDepth[,1])){
      lines(as.numeric(colnames(results$cumulativePedigreeDepth)),
             as.numeric(as.character(results$cumulativePedigreeDepth[x,])))
    }
    mtext("Cumulative distribution of pedigree depth.",side=1,line=5)
    mtext("Earliest cohorts are represented by the lowest lines.",side=1,line=6)
    inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
    if(inp=='s'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      savePlot(paste(s,".jpeg",sep=""),type="jpeg")
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }
  }

  if(is.null(cohorts)==FALSE) {
    thickness<-1
    if(grContrast==TRUE) thickness<-1.7
    sibCohortData<-as.data.frame(cbind(
               as.numeric(as.character(results$fullSibsByCohort)),
               as.numeric(as.character(results$maternalSibsByCohort)),
               as.numeric(as.character(results$paternalSibsByCohort)),
               as.numeric(as.character(results$maternalSibsByCohort))
                            -as.numeric(as.character(results$fullSibsByCohort)),
               as.numeric(as.character(results$paternalSibsByCohort))
                            -as.numeric(as.character(results$fullSibsByCohort)) ))
    plot(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,1],type='l',ylim=c(0,max(sibCohortData)), xlab="Cohort",ylab="Count",lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,2],lty='dashed',col=col1,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,4],col=col1,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,3],lty='dashed',col=col2,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,5],col=col2,lwd=thickness)
    mtext("Known sib pairs by cohort.  Black line: full sibs, red: maternal sibs,",side=1,line=5)
    mtext("blue: paternal sibs, dashed: total maternal and paternal",side=1,line=6)
    mtext("bs, solid coloured: maternal and paternal half sibs.",side=1,line=7)
    inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
    if(inp=='s'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      savePlot(paste(s,".jpeg",sep=""),type="jpeg")
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }

    if(inp=='e'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      postscript(paste(s,".eps",sep=""),width=8,height=8,horizontal=FALSE)
      par(oma=c(5,1,1,1))
    plot(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,1],type='l',ylim=c(0,max(sibCohortData)), xlab="Cohort",ylab="Count",lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,2],lty='dashed',col=col1,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,4],col=col1,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,3],lty='dashed',col=col2,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),sibCohortData[,5],col=col2,lwd=thickness)
    mtext("Known sib pairs by cohort.  Black line: full sibs, red: maternal sibs,",side=1,line=5)
    mtext("blue: paternal sibs, dashed: total maternal and paternal",side=1,line=6)
    mtext("bs, solid coloured: maternal and paternal half sibs.",side=1,line=7)
      dev.off()
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }


  }

  if(is.null(cohorts)==FALSE) {
    thickness<-1
    if(grContrast==TRUE) thickness<-1.7
    ymax<-max(max(results$maternalGrandmothersByCohort),max(results$maternalGrandfathersByCohort),max(results$paternalGrandmothersByCohort),max(results$paternalGrandfathersByCohort))
    plot(as.numeric(names(results$fullSibsByCohort)),as.numeric(as.character(results$maternalGrandmothersByCohort)),col=col1,type='l',ylim=c(0,ymax), xlab="Pedigree depth",ylab="Count",lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),as.numeric(as.character(results$maternalGrandfathersByCohort)),lty='dashed',col=col1,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),as.numeric(as.character(results$paternalGrandmothersByCohort)),col=col2,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),as.numeric(as.character(results$paternalGrandfathersByCohort)),lty='dashed',col=col2,lwd=thickness)
    mtext("Known grandparents by cohort (i.e., cohort of grand-offspring).",side=1,line=5)
    mtext("Red: maternal grandparents, blue: paternal grand",side=1,line=6)
    mtext("parents, solid: grandmothers, dashed: grandfathers.",side=1,line=7)
    inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
    if(inp=='s'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      savePlot(paste(s,".jpeg",sep=""),type="jpeg")
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }

    if(inp=='e'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      postscript(paste(s,".eps",sep=""),width=8,height=8,horizontal=FALSE)
      par(oma=c(5,1,1,1))
    plot(as.numeric(names(results$fullSibsByCohort)),as.numeric(as.character(results$maternalGrandmothersByCohort)),col=col1,type='l',ylim=c(0,ymax), xlab="Pedigree depth",ylab="Count",lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),as.numeric(as.character(results$maternalGrandfathersByCohort)),lty='dashed',col=col1,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),as.numeric(as.character(results$paternalGrandmothersByCohort)),col=col2,lwd=thickness)
    lines(as.numeric(names(results$fullSibsByCohort)),as.numeric(as.character(results$paternalGrandfathersByCohort)),lty='dashed',col=col2,lwd=thickness)
    mtext("Known grandparents by cohort (i.e., cohort of grand-offspring).",side=1,line=5)
    mtext("Red: maternal grandparents, blue: paternal grand",side=1,line=6)
    mtext("parents, solid: grandmothers, dashed: grandfathers.",side=1,line=7)
      dev.off()
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }




  }

if(lowMem==FALSE){

  relatednessInterval<-as.numeric(names(results$relatednessCategories)[3])-as.numeric(names(results$relatednessCategories)[2])
  midBins<-as.numeric(names(results$relatednessCategories))
  binLabels<-paste(midBins-relatednessInterval/2,"-",midBins+relatednessInterval/2,sep="")
  binLabels[1]<-paste("0-",relatednessInterval/2,sep="")
  for(x in seq(2,length(binLabels),by=2)) binLabels[x]<-""
  plotbins<-NULL
  for(x in 1:length(results$relatednessCategories)) { if(is.na(results$relatednessCategories[x])==FALSE) plotbins<-x }
  plotbins<-plotbins+1
  par(mar=c(10, 4, 4, 2))
  barplot(results$relatednessCategories,ylab="Count",xaxt='n')
  axis(1,at=(0:(plotbins-2))*1.2+0.6,labels=binLabels[1:(plotbins-1)],las=3,xlim=c(0,length(binLabels)-1))
  mtext("Non-zero pairwise relatednesses       ",1,7)
  par(mar=c(5, 4, 4, 2)+0.1)
  mtext("Distribution of relatedness across the pedigree.",side=1,line=5)
  inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
  if(inp=='s'){
    s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
    savePlot(paste(s,".jpeg",sep=""),type="jpeg")
    readline(prompt = "File saved.  Press <Enter> to continue...")
  }

}

  hist(results$maternalSibships[,2],xlab="(non-zero) maternal sibship sizes",ylab="count",main="")
  inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
  if(inp=='s'){
    s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
    savePlot(paste(s,".jpeg",sep=""),type="jpeg")
    readline(prompt = "File saved.  Press <Enter> to continue...")
  }

  hist(results$paternalSibships[,2],xlab="(non-zero) paternal sibship sizes",ylab="count",main="")
  inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
  if(inp=='s'){
    s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
    savePlot(paste(s,".jpeg",sep=""),type="jpeg")
    readline(prompt = "File saved.  Press <Enter> to continue...")
  }

  plot(x=1:(length(missingness[1,])-1),y=colMeans(as.matrix(missingness[,2:length(missingness[1,])])),type='l',
                                        xlab="Ancestral generation",ylab="Mean indiviual pedigree completeness")
  mtext("Pedigree-wide average ancestral pedigree completness",side=1,line=5)
  mtext("following MacCluer et al. (1983) J Hered 74:394-399.",side=1,line=6)
  inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
  if(inp=='s'){
    s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
    savePlot(paste(s,".jpeg",sep=""),type="jpeg")
    readline(prompt = "File saved.  Press <Enter> to continue...")
  }

  if(is.null(cohorts)==FALSE) {
    gc<-gray.colors(length(cohortPedigreeCompleteness[,1])-1)
    plot(x=1:(length(missingness[1,])-1),y=cohortPedigreeCompleteness[length(cohortPedigreeCompleteness[,1]),],
                     type='l',xlab="Ancestral generation",ylab="Mean individual pedigree completeness by cohort",
                     ylim=c(0,max(cohortPedigreeCompleteness)))
    for(i in 1:(length(cohortPedigreeCompleteness[,1])-1))
           lines(x=1:(length(missingness[1,])-1),y=cohortPedigreeCompleteness[i,],
                     col=gc[length(cohortPedigreeCompleteness[,1])-i])
    mtext("Average ancestral pedigree completeness by cohort",side=1,line=5)
    mtext("following MacCluer et al. (1983) J Hered 74:394-399.",side=1,line=6)
    mtext("Most recent cohorts are represented by the darkest lines.",side=1,line=7)
    inp<-readline(prompt = "Press <s> to save current plot or press <Enter> to continue...") 
    if(inp=='s'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      savePlot(paste(s,".jpeg",sep=""),type="jpeg")
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }

    if(inp=='e'){
      s<-readline(prompt = "Enter path (including file name but not extension) to which to save image: ") 
      postscript(paste(s,".eps",sep=""),width=8,height=8,horizontal=FALSE)
      par(oma=c(5,1,1,1))
      plot(x=1:(length(missingness[1,])-1),y=cohortPedigreeCompleteness[length(cohortPedigreeCompleteness[,1]),],
                     type='l',xlab="Ancestral generation",ylab="Mean individual pedigree completeness by cohort",
                     ylim=c(0,max(cohortPedigreeCompleteness)))
      for(i in 1:(length(cohortPedigreeCompleteness[,1])-1))
           lines(x=1:(length(missingness[1,])-1),y=cohortPedigreeCompleteness[i,],
                     col=gc[length(cohortPedigreeCompleteness[,1])-i])
      mtext("Average ancestral pedigree completeness by cohort",side=1,line=5)
      mtext("following MacCluer et al. (1983) J Hered 74:394-399.",side=1,line=6)
      mtext("Most recent cohorts are represented by the darkest lines.",side=1,line=7)
      dev.off()
      readline(prompt = "File saved.  Press <Enter> to continue...")
    }

  }



  inp<-readline(prompt = "Do you want to view pedigree images now?  This can take some time for large datasets, and more advanced options are available in drawPedigree().  Type <y> to print pedigree images or press <Enter> otherwise...") 

  if(inp=='y'&is.null(cohorts)==FALSE){
    cat("Generating pedigree image...")
    dev.new()
    drawPedigree(Ped=results$analyzedPedigree,cohorts=results$cohorts,writeCohortLabels='y')
    cat("done.")
    cat('\n')
    readline(prompt = "Pause. Press <Enter> to continue...") 
    cat("Generating image of maternal pedigre links...")
    dev.new()
    drawPedigree(Ped=results$analyzedPedigree,cohorts=results$cohorts,writeCohortLabels='y',links='mums')
    cat("done.")
    cat('\n')
    readline(prompt = "Pause. Press <Enter> to continue...") 
    cat("Generating image of paternal pedigre links...")
    dev.new()
    drawPedigree(Ped=results$analyzedPedigree,cohorts=results$cohorts,writeCohortLabels='y',links='dads')
    cat("done.")
    cat('\n')
  }

  if(inp=='y'&is.null(cohorts)){
    cat("Generating pedigree image...")
    dev.new()
    drawPedigree(Ped=results$analyzedPedigree)
    cat("done.")
    cat('\n')
    readline(prompt = "Pause. Press <Enter> to continue...") 
    cat("Generating image of maternal pedigre links...")
    dev.new()
    drawPedigree(Ped=results$analyzedPedigree,links='mums')
    cat("done.")
    cat('\n')
    readline(prompt = "Pause. Press <Enter> to continue...") 
    cat("Generating image of paternal pedigre links...")
    dev.new()
    drawPedigree(Ped=results$analyzedPedigree,links='dads')
    cat("done.")
    cat('\n')
  }
}##end of the graphicalReport=='y' if statement

  results
}

