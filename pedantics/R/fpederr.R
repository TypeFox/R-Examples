`fpederr` <-
function(truePedigree,founders=NULL,sex=NULL,samp=NULL,sireE=NULL,damE=NULL,
        sireA=NULL,damA=NULL,cohort=NULL,first=NULL,last=NULL,monoecey=0,modifyAssumedPedigree=0){

  fpederrR <- function(Ped,found,sex,samp,dade,mume,dadA,mumA,cohort,first,last,mono)
   .C("FPEDERR_R",
     as.integer(Ped$id),
     as.integer(Ped$dam),
     as.integer(Ped$sire),
     as.integer(found),
     as.integer(sex),
     as.integer(samp),
     as.double(dade),
     as.double(mume),
     as.double(dadA),
     as.double(mumA),
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
      founders<-rep(1,length(truePedigree[,1]))
      warnings<-rbind(warnings,"All individuals set to founder status")
    }
    if(is.null(sex)) {
      sex<-rep(-1,length(truePedigree[,1]))
      warnings<-rbind(warnings,"All individuals sex set to unknown")
    }
    if(is.null(samp)) 
      samp<-rep(1,length(truePedigree[,1]))
    if(is.null(sireE)) 
      sireE<-rep(0,length(truePedigree[,1]))
    if(is.null(damE))
      damE<-rep(0,length(truePedigree[,1]))
    if(is.null(sireA)) 
      sireA<-rep(1,length(truePedigree[,1]))
    if(is.null(damA))
      damA<-rep(1,length(truePedigree[,1]))
    if(is.null(cohort)) {
      cohort<-rep(1,length(truePedigree[,1]))
      warnings<-rbind(warnings,"No cohort designations - impossible pedigree loops could result")
    }
    if(is.null(first)) 
      first<-rep(1,length(truePedigree[,1]))
    if(is.null(last)) 
      last<-rep(1,length(truePedigree[,1]))

    if(length(sireE)==1) sireE<-rep(sireE,length(truePedigree[,1]))
    if(length(damE)==1) damE<-rep(damE,length(truePedigree[,1]))
    if(length(sireA)==1) sireA<-rep(sireA,length(truePedigree[,1]))
    if(length(damA)==1) damA<-rep(damA,length(truePedigree[,1]))
   
    numericPed<-makePedigreeNumeric(as.character(truePedigree$id),as.character(truePedigree$sire),
                           as.character(truePedigree$dam),missingVal=-1)

    fpederrRrun<-fpederrR(Ped=numericPed$numericPedigree,
               found=founders,sex=sex,samp=samp,
               dade=sireE,mume=damE,dadA=sireA,mumA=damA,
               cohort=cohort,first=first,last=last,mono=0)

    assumedPed<-makePedigreeFactor(fpederrRrun[[1]],fpederrRrun[[16]],fpederrRrun[[15]],numericPed$idKey)

    supplementalPed<-NULL
    if(modifyAssumedPedigree==1){
    supplementalPed<-assumedPed
      for(i in 1:length(assumedPed[,1])){
        if(is.na(assumedPed$sire[x])==FALSE) supplementalPed$sire[x] <- truePedigree$sire[x]
        if(is.na(assumedPed$dam[x])==FALSE) supplementalPed$dam[x] <- truePedigree$dam[x]
      }
    }
    if(modifyAssumedPedigree==2){
      supplementalPed<-makePedigreeFactor(fpederrRrun[[1]],fpederrRrun[[18]],fpederrRrun[[17]],numericPed$idKey)
    }

    fpederrResults<-list(assumedPedigree=assumedPed,truePedigree=truePedigree)
    if(is.null(supplementalPed)==FALSE) fpederrResults<-list(assumedPedigree=assumedPed,truePedigree=truePedigree,
                                                           supplementalPedigree=supplementalPed)
    if(is.null(warnings)==FALSE) {
      for(x in 1:length(warnings)) cat(paste("Warning: ",warnings[x],'\n'));
    }

    fpederrResults
}

