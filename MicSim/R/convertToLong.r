####################################################################################
####################################################################################
## FUNCTION TO CONVERT MICROSIMULATION OUTPUT INTO LONG FORMAT                    ##
## SZ, December 2013                                                              ##
####################################################################################
####################################################################################
# Function requires as global variables (from microsimulation definition):  
# `absTranstions': matrix comprising information on all transitions to absorbent states
# `allTransitions': matirix comprising information on all transitions possible between values of state variables 
# `simHorizon': simulation horizon (i.e., start and end date)

convertToLongFormat <- function(pop, migr=FALSE) {
  
  # Take global variables used from global environment
  absTransitions <- mget('absTransitions', envir=globalenv(), ifnotfound=list(0))$absTransitions
  allTransitions <- mget('allTransitions', envir=globalenv(), ifnotfound=list(0))$allTransitions
  stateSpace <- mget('stateSpace', envir=globalenv(), ifnotfound=list(0))$stateSpace
  simHorizon <- mget('simHorizon', envir=globalenv(), ifnotfound=list(0))$simHorizon
  if(migr)
    immigrPop <- mget('immigrPop', envir=globalenv(), ifnotfound=list(0))$immigrPop
  
  if(is.vector(absTransitions))
    absTransitions <- matrix(unlist(absTransitions), ncol=2, nrow=1)  
  absStates <- absTransitions[,1]  
  if(is.null(dim(stateSpace))){    
    stateSpaceTMP <- as.data.frame(matrix(unlist(stateSpace), ncol=1)) 
    colnames(stateSpaceTMP) <- attr(stateSpace,'name')
    stateSpace <- stateSpaceTMP
  }
    
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # Set of auxiliary function used in transformation process 
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # Check whether given year is leap year
  isLeapYear <- function(year){
    return(ifelse(year%%400==0, T, ifelse(year%%4==0 & year%%100!=0, T,F)))}
  # Return `dAte' as year
  exactYear <- function(dAte){   
    if(is.na(dAte))
      return(NA)
    year <- as.numeric(as.character(years(dAte)))
    daysTo <- as.numeric(chron(paste(c('31/12',year-1),collapse='/'), format=c(dates='d/m/Y', times='h:m:s')))
    daysInYear <- ifelse(isLeapYear(year),366,365)
    prop <- c(as.numeric(dAte)-daysTo)/daysInYear
    return(year+prop)
  }
  # Extract transition between values of state variables (using information given in `absTranstions' and `allTransitions')
  getTransition <- function(oSdS){
     oS <- oSdS[1]
     dS <- oSdS[2]
     tr <- ''
     if(is.na(oS) | is.na(dS)){
      tr <- 'cens'
     } else {
       oS <- as.character(unlist(oS))
       dS <- as.character(unlist(dS))
       if(dS %in% absTransitions[,1]){
        tr <- dS
       } else {
        oS <- unlist(strsplit(oS,split='/'))
        dS <- unlist(strsplit(dS,split='/'))
        cid <- which(allTransitions[,1] %in% paste(oS,dS, sep='->'))
        if(length(cid)==1){
          tr <- allTransitions[cid,1]
        } else {   
          cad <- which(oS!=dS)
          tr <- paste(oS[cad],dS[cad],sep='->')
        }
       }     
     }
     return(tr)   
  }
  # According to transition given, create state of destination by replacing value of state variable 
  replaceStateValue <- function(vec){
   tr <- as.character(unlist(vec[length(vec)]))
   tr <- unlist(strsplit(tr,split='->'))
   att <- vec[-length(vec)]
   att[which(att==tr[1])] <- tr[2]
   return(att)
  }
  # Create sequence of numbers
  giveSeq <- function(nu){
     return(1:nu)
  }
 
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # Transform output of microsimulation into data in long format
  # --------------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------------- 
  id <- pop[,'ID', drop=F]
  birthDate <- pop[,'birthDate']
  birthyear <- sapply(pop[,'birthDate'],exactYear)
  trYear <- sapply(pop[,'transitionTime'],exactYear)
  trDate <- pop[,'transitionTime'] 
  Tstop <- pop[,'transitionTime']
  date <- Tstop
  popLong <- data.frame(id,birthDate,birthyear,trYear,trDate,Tstop,stringsAsFactors=FALSE)
  popLong <- popLong[order(as.numeric(popLong[,'ID'])),]
  fromState <- do.call(rbind,strsplit(as.character(pop[,'From']), split='/'))
  initState <- do.call(rbind,strsplit(as.character(pop[,'initState']),split='/'))
  idMissVar <- which(is.na(fromState[,1]))
  fromState[idMissVar,] <- initState[idMissVar,] 
  # Create for each state variable considered a column in the `newly' constructed data set 
  for(i in 1:dim(stateSpace)[2]){
   nam <- names(stateSpace)[i]
   colNam <- c(colnames(popLong), nam)
   popLong <- cbind(popLong,fromState[,i])
   colnames(popLong) <- colNam
  }
  # Describe observation scheme: 
  # statusEntry: 0: left trunction, 1: entry observed
  statusEntry <- rep(1,length(fromState[,1])) 
  # statusExit: 0: right censoring, 1: (relevant) transition observed
  statusExit <- rep(1,length(fromState[,1]))
  statusExit[idMissVar] <- rep(0, length(idMissVar))
  popLong <- data.frame(popLong, statusEntry, statusExit, stringsAsFactors=FALSE)
  # There are newborns who only experienced birth during simulation. 
  # These individuals are (right) censored at simulation end. 
  isCensNB <- which(is.na(popLong[,'Tstop']))
  if(length(isCensNB)>0)
    popLong[isCensNB,'Tstop'] <- simHorizon[2]  
  # Identify from information given in the combined states the transitions between values of the state variables.
  OD <- apply(pop[,c('From','To')],1,getTransition)
  popLong <- data.frame(popLong, OD, stringsAsFactors=FALSE)
  # Beware of transitions into absorbent states and censored transitions. 
  # Furthermore, account for the fact that all individuals who are still alive at simulation end will experience
  # a censoring event then.
  idSet <- popLong[which((popLong[,'OD'] %in% absTransitions[,1]) | (popLong[,'statusExit'] == 0)),'ID']
  idCensSE <- setdiff(unique(popLong[,'ID']),idSet)
  if(length(idCensSE)>0){
    popLongCensSE <- popLong[popLong[,'ID'] %in% idCensSE,]
    popLongCensSE <- popLongCensSE[which(c(diff(as.numeric(popLongCensSE[,'ID'])),2)!=0),]
    odId <- which(names(popLong)=='OD')
    ReplMat <- apply(popLongCensSE[,c(7:(6+dim(stateSpace)[2]),odId)],1, replaceStateValue) 
    popLongCensSE[,7:(6+dim(stateSpace)[2])] <- t(ReplMat)
    popLongCensSE[,'statusExit'] <- 0
    popLongCensSE[,'OD'] <- 'cens'
    popLongCensSE[,'Tstop'] <- simHorizon[2]
    popLong <- rbind(popLong, popLongCensSE)
  }
  popLong <- popLong[order(as.numeric(popLong[,'ID'])),]
  # Extract episode starting time and age at episode starting time.
  popLong$Tstart <- c(chron(NA),popLong$Tstop[-length(popLong[,1])])
  suppressWarnings(popLong$Tstart[!duplicated(popLong[,'ID'])] <- chron(NA))
  # Beware we have not observed the episode start of individuals being part of the starting population (=left truncation). 
  simStartYear <- exactYear(simHorizon[1])  
  bornInSimId <- which(is.na(popLong$Tstart) & popLong$birthyear >= simStartYear)
  bornOutSimId <- which(is.na(popLong$Tstart)& popLong$birthyear < simStartYear)
  popLong$Tstart[bornInSimId] <- popLong$birthDate[bornInSimId] 
  popLong$Tstart[bornOutSimId] <- simHorizon[1]
  popLong$statusEntry[bornOutSimId] <- 0
  # The timing of the first episodes reported for immigrants has to be adjusted to account for the fact that we have not
  # observed when they entered the state which they occupied at immigration.
  if(migr){
    imPopLong <- popLong[popLong[,'ID'] %in% immigrPop[,'ID'],]
    imPopLongFEID <- which((popLong[,'ID'] %in% immigrPop[,'ID']) & !duplicated(popLong[,'ID']))
    imPopLongFE <- popLong[imPopLongFEID,]
    imPopLongFE <- merge(imPopLongFE,immigrPop[,-3], by='ID')
    imPopLongFE$Tstart <- imPopLongFE[,'immigrDate']
    imPopLongFE$statusEntry <- 1
    imId <- which(names(imPopLongFE) %in% c('immigrDate','immigrInitState'))
    imPopLongFE <- imPopLongFE[,-imId]
    popLong <- popLong[-imPopLongFEID,]
    popLong <- rbind(popLong, imPopLongFE)
  }
  popLong <- popLong[order(popLong[,c('Tstart')]),]
  popLong <- popLong[order(as.numeric(popLong[,c('ID')])),]
  popLong <- popLong[,-which(names(popLong)=='trYear')]
  # Count episodes each individual experienced.
  ns <-  data.frame(table(popLong$ID),stringsAsFactors=FALSE)
  ns <- ns[order(as.numeric(as.character(ns[,1]))),]
  colnames(ns) <- c('ID','ns')
  popLong <- merge(popLong, ns, by='ID')
  popLong <- popLong[order(as.numeric(popLong[,c('ID')])),] 
  # For each individual enumerate episodes.
  nsU <- popLong$ns[which(!duplicated(popLong[,'ID']))]
  popLong$Episode <- unlist(sapply(nsU,giveSeq))
  # Select information to return.
  idCovs <- 6:(5+dim(stateSpace)[2])
  idC <- which(colnames(popLong) %in% c('ID','birthDate','Tstart','Tstop','OD','ns','statusEntry','statusExit','Episode'))
  popLong <- popLong[,c(idC,idCovs)]   
  idCovs <- 10:(9+dim(stateSpace)[2])
  popLong <- popLong[,c(1,2,7,3,4,5,6,8,9,idCovs)]
  popLong[,'Tstart'] <- chron(popLong[,'Tstart'], format=c(dates='d/m/y'), out.format=c(dates='d/m/year'))
  return(popLong)
}














