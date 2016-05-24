#"FAMID","ID","PID","MID","SEX",ZYGOSITY,ADOPTED,INDEP
FSV.frompedi <- function(pedi.dat,phen.dat){
  if(!is.data.frame(pedi.dat)){pedi.dat <- data.frame(pedi.dat)}
  if(!is.data.frame(phen.dat)){phen.dat <- data.frame(phen.dat)}
  extracols <- c("ZYGOSITY","ADOPTED","INDEP") %in% colnames(pedi.dat)
  ftype <- rep(NA,nrow(pedi.dat))
  indiv <- rep(NA,nrow(pedi.dat))
  
  if( !all( c("FAMID","ID","PID","MID","SEX") %in% colnames(pedi.dat) ) | ncol(pedi.dat)<5 ){
    warning("No compatible family structure is identifiable from pedigree data; be sure its column names are correct.")
    phen.dat$FTYPE <- 6
    phen.dat$INDIV <- 1
    return(phen.dat)
  }
  else{
    if( !all(pedi.dat[,c("FAMID","ID")]==pedi.dat[order(pedi.dat$FAMID, pedi.dat$ID),c("FAMID","ID")]) ){
      pedi.dat <- pedi.dat[order(pedi.dat$FAMID, pedi.dat$ID),]
    }
    pedi.dat$PID[is.na(pedi.dat$PID)] <- 0
    pedi.dat$MID[is.na(pedi.dat$MID)] <- 0
  }       
  
  if(all(extracols[1:2]==c(FALSE,FALSE))){ #"FAMID","ID","PID","MID","SEX"
    ftype <- rep(4,nrow(pedi.dat))
    if(extracols[3]){
      ftype[ !is.na(pedi.dat$INDEP) & pedi.dat$INDEP==1 ] <- 6
      indiv[ !is.na(pedi.dat$INDEP) & pedi.dat$INDEP==1 ] <- 1
    }
    knownparents <- which(pedi.dat$ID %in% c(pedi.dat$MID,pedi.dat$PID))
    indiv[pedi.dat$ID %in% pedi.dat$MID] <- 3 #known moms
    indiv[pedi.dat$ID %in% pedi.dat$PID] <- 4 #known dads
    indiv[pedi.dat$MID!=0 | pedi.dat$PID!=0] <- 1 #known kids
    indiv[is.na(indiv)] <- 9 #temporary value
    #indiv[ftype!=6] <- unlist(apply(X=pedi.dat[ftype!=6,c("FAMID","ID","PID","MID","SEX")],MARGIN=1,FUN=function(x){
      #if(x[3]!=0 & x[4]!=0){return(1)} #offspring
      #if(x[3]==0 & x[4]==0 & x[5]==2){return(3)} #mom
      #if(x[3]==0 & x[4]==0 & x[5]==1){return(4)} #dad
    #}))
    indiv[which(duplicated(paste(pedi.dat$"FAMID",indiv,sep=".")) & indiv==1 & ftype!=6)] <- 2 #Turn extra 1s into 2s.
    extrakids <- which(duplicated(paste(pedi.dat$"FAMID",indiv,sep=".")) & indiv==2) #Turn extra 2s into indobs.
    indiv[extrakids] <- 1
    ftype[extrakids] <- 6
    extraparents <- which( #duplicated(paste(pedi.dat$"FAMID",indiv,sep=".")) &
      (paste(pedi.dat$"FAMID",indiv,sep=".") %in% 
        paste(pedi.dat$"FAMID",indiv,sep=".")[duplicated(paste(pedi.dat$"FAMID",indiv,sep="."))]) &
        (indiv %in% c(3,4)) & !(pedi.dat$ID %in% pedi.dat$ID[knownparents]) & ftype!=6 )
    ftype[extraparents] <- 6 #Turn extra parents into indobs.
    indiv[extraparents] <- 1
    ftype[indiv==9] <- 6 #turn temporary 9s into indobs.
    indiv[ftype==6] <- 1 #make sure all indobs are INDIV=1.
    names(ftype) <- pedi.dat$ID
    names(indiv) <- pedi.dat$ID
    if(length(extrakids)>0){
      warning("Some non-founders coerced to family-type 6 (treated as 'independent observations').")
    }
    phen.dat$FTYPE <- ftype[as.character(phen.dat$ID)]
    phen.dat$INDIV <- indiv[as.character(phen.dat$ID)]
    return(phen.dat)
  }
  
  if(all(extracols[1:2]==c(TRUE,FALSE))){#"FAMID","ID","PID","MID","SEX",ZYGOSITY
    ftype <- rep(4,nrow(pedi.dat))
    if(extracols[3]){
      ftype[ !is.na(pedi.dat$INDEP) & pedi.dat$INDEP==1 ] <- 6
      indiv[ !is.na(pedi.dat$INDEP) & pedi.dat$INDEP==1 ] <- 1
    }
    mzfams <- unique(subset(pedi.dat$FAMID, !is.na(pedi.dat$ZYGOSITY) & pedi.dat$ZYGOSITY==1))
    dzfams <- unique(subset(pedi.dat$FAMID, !is.na(pedi.dat$ZYGOSITY) & pedi.dat$ZYGOSITY==2))
    ftype[(pedi.dat$FAMID %in% mzfams) & ftype!=6] <- 1
    ftype[(pedi.dat$FAMID %in% dzfams) & ftype!=6] <- 2
    indiv[!is.na(pedi.dat$ZYGOSITY) & (pedi.dat$ZYGOSITY %in% c(1,2))] <- 1 #known twins
    knownparents <- which(pedi.dat$ID %in% c(pedi.dat$MID,pedi.dat$PID))
    indiv[pedi.dat$ID %in% pedi.dat$MID] <- 3 #known moms
    indiv[pedi.dat$ID %in% pedi.dat$PID] <- 4 #known dads
    indiv[(pedi.dat$MID!=0 | pedi.dat$PID!=0) & is.na(indiv) ] <- 1 #known kids, not already known to be twins
    indiv[is.na(indiv)] <- 9 #temporary value
    indiv[which(duplicated(paste(pedi.dat$"FAMID",indiv,sep=".")) & indiv==1 & ftype!=6)] <- 2 #Turn extra 1s into 2s.
    extrakids <- which(duplicated(paste(pedi.dat$"FAMID",indiv,sep=".")) & indiv==2) #Turn extra 2s into indobs.
    indiv[extrakids] <- 1
    ftype[extrakids] <- 6
    extraparents <- which( #duplicated(paste(pedi.dat$"FAMID",indiv,sep=".")) &
      (paste(pedi.dat$"FAMID",indiv,sep=".") %in% 
         paste(pedi.dat$"FAMID",indiv,sep=".")[duplicated(paste(pedi.dat$"FAMID",indiv,sep="."))]) &
        (indiv %in% c(3,4)) & !(pedi.dat$ID %in% pedi.dat$ID[knownparents]) & ftype!=6 )
    ftype[extraparents] <- 6 #Turn extra parents into indobs.
    indiv[extraparents] <- 1
    ftype[indiv==9] <- 6 #turn temporary 9s into indobs.
    indiv[ftype==6] <- 1 #make sure all indobs are INDIV=1.
    names(ftype) <- pedi.dat$ID
    names(indiv) <- pedi.dat$ID
    if(length(extrakids)>0){
      warning("Some non-founders coerced to family-type 6 (treated as 'independent observations').")
    }
    phen.dat$FTYPE <- ftype[as.character(phen.dat$ID)]
    phen.dat$INDIV <- indiv[as.character(phen.dat$ID)]
    return(phen.dat)
  }
  
  if(extracols[2]==TRUE){
    if(extracols[1]==TRUE){
      mzfams <- unique(subset(pedi.dat$FAMID, !is.na(pedi.dat$ZYGOSITY) & pedi.dat$ZYGOSITY==1))
      dzfams <- unique(subset(pedi.dat$FAMID, !is.na(pedi.dat$ZYGOSITY) & pedi.dat$ZYGOSITY==2))
      ftype[(pedi.dat$FAMID %in% mzfams)] <- 1
      ftype[(pedi.dat$FAMID %in% dzfams)] <- 2
    }
    if(extracols[3]==TRUE){
      ftype[ !is.na(pedi.dat$INDEP) & pedi.dat$INDEP==1 ] <- 6
      indiv[ !is.na(pedi.dat$INDEP) & pedi.dat$INDEP==1 ] <- 1
    }
    uniqfams <- unique(pedi.dat$FAMID)
    ambiguousflag <- FALSE
    
    for(i in 1:length(uniqfams)){
      which.rows.curr <- which(pedi.dat$FAMID==uniqfams[i])
      fam.curr <- pedi.dat[which.rows.curr,]
      ftype.curr <- ftype[which.rows.curr]
      indiv.curr <- indiv[which.rows.curr]
      dadflag <- momflag <- kid1flag <- kid2flag <- adoptflag <- bioflag <- FALSE
      
      for(j in 1:length(which.rows.curr)){
        if(!is.na(ftype.curr[j]) & ftype.curr[j]==6){next} #If already classified as type 6
        if(fam.curr[j,"PID"]!=0 & fam.curr[j,"MID"]!=0){ #If nonfounder
          if(kid1flag==FALSE){indiv.curr[j] <- 1; kid1flag <- bioflag <- TRUE; next}
          if(kid1flag==TRUE & kid2flag==FALSE){indiv.curr[j] <- 2; kid2flag <- bioflag <- TRUE; next }
          if(kid1flag+kid2flag==2){indiv.curr[j] <- 1; ftype.curr[j] <- 6; ambiguousflag <- TRUE; next}
        }
        else{ #If founder
          if( !is.na(fam.curr[j,"ADOPTED"]) & fam.curr[j,"ADOPTED"]==1){ #If adopted
            if(kid2flag==FALSE){indiv.curr[j] <- 2; kid2flag <- adoptflag <- TRUE; next}
            if(kid2flag==TRUE & kid1flag==FALSE){indiv.curr[j] <- 1; kid1flag <- adoptflag <- TRUE; next }
            if(kid1flag+kid2flag==2){indiv.curr[j] <- 1; ftype.curr[j] <- 6; ambiguousflag <- TRUE; next}
          }
          else{ #If not adopted
            if(fam.curr[j,"SEX"]==2){ #If female
              if(momflag==FALSE){indiv.curr[j] <- 3; momflag <- TRUE; next}
              else{indiv.curr[j] <- 1; ftype.curr[j] <- 6; next}
            }
            if(fam.curr[j,"SEX"]==1){ #If male
              if(dadflag==FALSE){indiv.curr[j] <- 4; dadflag <- TRUE; next}
              else{indiv.curr[j] <- 1; ftype.curr[j] <- 6; next}
      }}}}
      if(any(is.na(ftype.curr))){
        if(adoptflag==TRUE){
          if(bioflag==TRUE){ftype.curr[is.na(ftype.curr)] <- 5}
          else{ftype.curr[is.na(ftype.curr)] <- 3}
        }
        else{ftype.curr[is.na(ftype.curr)] <- 4}
      }
      ftype[which.rows.curr] <- ftype.curr
      indiv[which.rows.curr] <- indiv.curr
    }
    if(ambiguousflag==TRUE){warning("Some non-founders coerced to family-type 6 (treated as 'independent observations').")}
    names(ftype) <- pedi.dat$ID
    names(indiv) <- pedi.dat$ID
    phen.dat$FTYPE <- ftype[as.character(phen.dat$ID)]
    phen.dat$INDIV <- indiv[as.character(phen.dat$ID)]
    return(phen.dat)
}}
