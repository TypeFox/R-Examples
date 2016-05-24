`rpederr` <-
function(assumedPedigree,founders=NULL,sex=NULL,samp=NULL,sireE=NULL,damE=NULL,
        sireS=NULL,damS=NULL,cohort=NULL,first=NULL,last=NULL,monoecey=0,modifyAssumedPedigree=0){

  rpederrR <- function(Ped,found,sex,samp,dade,mume,dads,mums,cohort,first,last,mono)
    .C("RPEDERR_R",
     as.integer(Ped$id),
     as.integer(Ped$dam),
     as.integer(Ped$sire),
     as.integer(found),
     as.integer(sex),
     as.integer(samp),
     as.double(dade),
     as.double(mume),
     as.double(dads),
     as.double(mums),
     as.integer(cohort),
     as.integer(first),
     as.integer(last),
     as.integer(length(Ped$id)),
     as.integer(Ped$dam*0),
     as.integer(Ped$sire*0),
     as.integer(Ped$dam*0),
     as.integer(Ped$sire*0),
     as.integer(mono),PACKAGE="pedantics")

    warnings<-NULL

    if(is.null(founders)) {
      founders<-rep(1,length(assumedPedigree[,1]))
      warnings<-rbind(warnings,"All individuals set to founder status")
    }
    if(is.null(sex)) {
      sex<-rep(-1,length(assumedPedigree[,1]))
      warnings<-rbind(warnings,"All individuals sex set to unknown")
    }
    if(is.null(samp)) 
      samp<-rep(1,length(assumedPedigree[,1]))
    if(is.null(sireE)) 
      sireE<-rep(0,length(assumedPedigree[,1]))
    if(is.null(damE))
      damE<-rep(0,length(assumedPedigree[,1]))
    if(is.null(sireS)) 
      sireS<-rep(1,length(assumedPedigree[,1]))
    if(is.null(damS))
      damS<-rep(1,length(assumedPedigree[,1]))
    if(is.null(cohort)) {
      cohort<-rep(1,length(assumedPedigree[,1]))
      warnings<-rbind(warnings,"No cohort designations - impossible pedigree loops could result")
    }
    if(is.null(first)) 
      first<-rep(1,length(assumedPedigree[,1]))
    if(is.null(last)) 
      last<-rep(1,length(assumedPedigree[,1]))
   
    numericPed<-makePedigreeNumeric(as.character(assumedPedigree$id),as.character(assumedPedigree$sire),
                           as.character(assumedPedigree$dam),missingVal=-1)


    rpederrRrun<-rpederrR(Ped=numericPed$numericPedigree,
               found=founders,sex=sex,samp=samp,
               dade=sireE,mume=damE,dads=sireS,mums=damS,
               cohort=cohort,first=first,last=last,mono=0)

    truePed<-makePedigreeFactor(rpederrRrun[[1]],rpederrRrun[[16]],rpederrRrun[[15]],numericPed$idKey)
    truePed<-truePed[,c(1,3,2)]

    supplementalPed<-NULL
    if(modifyAssumedPedigree==1){
    supplementalPed<-assumedPedigree
      for(i in 1:length(assumedPedigree[,1])){
        if(is.na(assumedPedigree$sire[x])==FALSE) supplementalPed$sire[x] <- truePed$sire[x]
        if(is.na(assumedPedigree$dam[x])==FALSE) supplementalPed$dam[x] <- truePed$dam[x]
      }
    }
    if(modifyAssumedPedigree==2){
      supplementalPed<-makePedigreeFactor(rpederrRrun[[1]],rpederrRrun[[18]],rpederrRrun[[17]],numericPed$idKey)[,c(1,3,2)]
    }

    rpederrResults<-list(assumedPedigree=assumedPedigree,truePedigree=truePed)
    if(is.null(supplementalPed)==FALSE) rpederrResults<-list(assumedPedigree=assumedPedigree,truePedigree=truePed,
                                                           supplementalPedigree=supplementalPed)
    if(is.null(warnings)==FALSE) {
      for(x in 1:length(warnings)) cat(paste("Warning: ",warnings[x],'\n'));
    }


    rpederrResults
}

