`rpederrBird` <-
function(assumedPedigree,founders=NULL,sex=NULL,samp=NULL,
        sireS=NULL,damS=NULL,cohort=NULL,first=NULL,last=NULL,broods=NULL,monoecey=0,
        broodData=NULL,EPPsireData=NULL,propEPPbroods=1,propEPPchicksGivenEPPbrood=0.5,
        propEPPbroodsTwoFathers=0,modifyAssumedPedigree=0,EPPlambda=0.001,malePhenotypes=NULL,EPPbeta=0,EPPgamma=0){

  if(is.null(malePhenotypes)) malePhenotypes<-as.data.frame(cbind(assumedPedigree[,1],0))

  rpederrR_bird <- function(Ped,found,sex,samp,dads,mums,cohort,first,last,mono,
                             numericBroods,numericEPPsireIDs,EPPsireData,numericBroodIDs,broodData,EPPlambda,
                             propEPPbroods,propEPPchicksGivenEPPbrood,propEPPbroodsTwoFathers,
                             mPhen,EPPbeta,EPPgamma)
    .C("RPEDERR_R_BIRD",
     as.integer(Ped$id),
     as.integer(Ped$dam),
     as.integer(Ped$sire),
     as.integer(found),
     as.integer(sex),
     as.integer(samp),
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
     as.integer(mono),
     as.integer(numericBroods),

     as.integer(length(EPPsireData[,1])),
     as.integer(numericEPPsireIDs),
     as.integer(as.character(EPPsireData[,2])),
     as.double(as.character(EPPsireData[,3])),
     as.double(as.character(EPPsireData[,4])),
     as.double(as.character(mPhen[,2])),
     as.double(EPPbeta),
     as.double(EPPgamma),
     as.integer(length(broodData[,1])),
     as.integer(as.character(broodData[,2])),
     as.double(as.character(broodData[,3])),
     as.double(as.character(broodData[,4])),
     as.integer(as.character(broodData[,5])),
     as.double(EPPlambda),
     as.double(propEPPbroods),
     as.double(propEPPchicksGivenEPPbrood),
     as.double(propEPPbroodsTwoFathers),PACKAGE="pedantics")

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

    EPPsireIDs<-as.data.frame(as.character(EPPsireData[,1]))
    EPPsireIDs$numericIDs<-numericPed$idKey[match(EPPsireIDs[,1],numericPed$idKey[,2]),1]

    all.broods<-c(as.character(broods),as.character(broodData[,1]))
    all.broods<-as.factor(all.broods)
    numericBroodConversion<-as.data.frame(cbind(seq(1:length(levels(all.broods))),levels(all.broods)))
    all.broods<-as.data.frame(all.broods)
    all.broods$numeric<-numericBroodConversion[match(all.broods[,1],numericBroodConversion[,2]),1]
    numericBroods<-as.numeric(as.character(all.broods$numeric[1:length(broods)]))
    for(x in 1:length(numericBroods)) { if(is.na(numericBroods[x])) numericBroods[x]<- -1 }
    numericBroodDataIDs<-as.numeric(as.character(all.broods$numeric[(length(broods)+1):length(all.broods$numeric)]))
    

    rpederrRrun<-rpederrR_bird(Ped=numericPed$numericPedigree,
               found=founders,sex=sex,samp=samp,
               dads=sireS,mums=damS,
               cohort=cohort,first=first,last=last,mono=0,
               numericBroods=numericBroods,numericEPPsireIDs=EPPsireIDs$numericIDs,
               EPPsireData=EPPsireData,numericBroodIDs=numericBroodDataIDs,broodData=broodData,EPPlambda=EPPlambda,
               propEPPbroods=propEPPbroods,propEPPchicksGivenEPPbrood=propEPPchicksGivenEPPbrood,propEPPbroodsTwoFathers=propEPPbroodsTwoFathers,
               mPhen=malePhenotypes,EPPbeta=EPPbeta,EPPgamma=EPPgamma)

    truePed<-makePedigreeFactor(rpederrRrun[[1]],rpederrRrun[[14]],rpederrRrun[[13]],numericPed$idKey)



   truePed
}

