####################################################################################
####################################################################################
## FUNCTION EXECUTING MICROSIMULATION                                             ##
## SZ, November 2013                                                              ##
####################################################################################
####################################################################################
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# I. Execute microsimulation as single thread
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
micSim <- function(initPop, immigrPop=NULL, transitionMatrix, absStates=NULL, initStates=c(), initStatesProb=c(), 
                   maxAge=99, simHorizon, fertTr=c(), dateSchoolEnrol='09/01') {
  
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # A. CHECK INPUT FOR CONSISTENCY
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------                   
  if(is.null(initPop))
    stop('No starting population has been defined.')
  if(!is.null(initPop)){
    if(paste(colnames(initPop),collapse='/')!='ID/birthDate/initState')
      stop('Matrix specifying the starting population has not been defined properly.')
  }
  if(!is.null(immigrPop)){
    if(paste(colnames(immigrPop),collapse='/')!='ID/immigrDate/birthDate/immigrInitState')
      stop('Matrix specifying immigrants has not been defined properly.')
  }
  if(is.null(transitionMatrix))
    stop('Matrix defining transition pattern und functions has not been defined properly.')
  if(maxAge<=0)
    stop('The maximal age until which individual life courses are simulated should exceed zero.')
  if(length(simHorizon)!=2)
    stop('The simulation horizon has not been defined properly.')
  if(class(simHorizon)[1]!='dates' & class(simHorizon)[2]!='dates')
    stop('The simulation horizon has not been defined properly.')   
  if(is.null(absStates))
    absStates <- setdiff(colnames(transitionMatrix),rownames(transitionMatrix))  
  if(length(fertTr)>0){
    if(is.null(initStates) | is.null(initStatesProb))
      stop('For children potentially born during simulation no inital state(s) and/or corresponding occurrence probabilities have been defined.')
  }
  if(length(dateSchoolEnrol)==0){
    dateSchoolEnrol <- '09/01'
  } 
  
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # B. DEFINITION OF GLOBAL PARAMETERS
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------   
  # Event queue
  queue <- matrix(NA,ncol=6,nrow=0)
  colnames(queue) <- c('ID','currTime','currState','currAge','nextState','timeToNextState')
  # Global time
  t.clock <- as.numeric(simHorizon[1])  # counts in days
  # Recording of transitions performedfer
  transitions <- matrix(NA,ncol=5,nrow=0)
  colnames(transitions) <- c('ID', 'From', 'To', 'transitionTime', 'transitionAge') 
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # C. FUNCTIONS REQUIRED FOR SIMULATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Function checks whether given `year' is a leap year.
  isLeapYear <- function(year) {  
    if (((year %% 4 == 0) & (year %% 100 != 0)) || (year %% 400 == 0))
      return(TRUE)    
    return(FALSE) 
  }
  
  # Function computes the correct age in years (accounting for leap years); arguments are the birth date and the current date
  getAgeInYears <- function(days) {    
    date <- chron(days)
    c.y <- as.numeric(as.character(years(date)))
    completeYears <-  c.y - 1970 
    daysInCY <- ifelse(isLeapYear(c.y), 366, 365)
    fracCY <- as.numeric(date - chron(paste(1,'/',1,'/',c.y), format=c(dates='d/m/y')) + 1)/daysInCY    
    return(completeYears + fracCY)
  }
  
  # Function checks whether a transition causes a newborn. (Demands `fert': all values of fertility variable, and 
  # `fertTr': matrix indicating transitions between fertility attributes.)
  isBirthEvent <- function(currState, destState){
    if(length(fertTr)==0)
      return(FALSE)
    fert <- matrix(unlist(strsplit(fertTr,split='->')), ncol=2)
    cS <- unlist(strsplit(currState,'/'))
    if(!("f" %in% cS))
      return(FALSE)      
    dS <- unlist(strsplit(destState,'/'))   
    for(i in 1:nrow(fert)){
      ff <- fert[i,]
      oS <- unlist(strsplit(ff[1],'/'))
      bS <- unlist(strsplit(ff[2],'/'))    
      if(!(F %in% (oS %in% cS)) & !(F %in% (bS %in% dS))){
        return(TRUE)
      } 
    }  
    return(FALSE)
  }
  
  # Funtion adds to simulation population a newborn.
  addNewNewborn <- function(birthTime){    
    birthState <- sample(apply(initStates,1,paste,collapse='/'),size=1,replace=T,prob=initStatesProb)
    if(is.null(immigrPop)){
      id <- as.numeric(max(as.numeric(initPop[,'ID'])))+1
    } else {
      id <- as.numeric(max(c(as.numeric(immigrPop[,'ID']),as.numeric(initPop[,'ID']))))+1
    }
    birthDate <- dates(chron(birthTime,
                             format=c(dates='d/m/Y', times='h:m:s'),out.format=c(dates='d/m/year', times='h:m:s')))    
    newInd <- c(id,as.character(birthDate),birthState)
    #cat('NewBorn: ',newInd,'\n')
    initPop <<- rbind(initPop,newInd)
    nE <- getNextStep(c(id,birthState,0,birthTime))
    #print(nE)
    #cat('\n------------n')
  }
  
  # Function checks whether a transition implies a school enrolment (in the year when child turns seven).
  # (If state 1 comprises value `no' and state 2 comprises value `low', the transition is marked as `school enrolment'.)
  isSchoolEnrolment <- function(currState,destState){
    enrol <- F 
    if(T %in% ('no' %in% unlist(strsplit(currState,'/'))) & T %in% ('low' %in% unlist(strsplit(destState,'/'))))
      enrol <- T
    return(enrol) 
  }
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # D. SIMULATION STEP
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Function to compute the next transition state and time of an individual who at age `currAge' and at time `calTime' 
  # entered its current state `currState'.
  getNextStep <- function(inp, isIMInitEvent=F){     
    # Transform input
    id <- as.numeric(unlist(inp[1]))
    currState <- as.character(unlist(inp[2]))
    currAge <- as.numeric(unlist(inp[3]))         # age in days
    calTime <- as.numeric(unlist(inp[4]))         # calendat time in days
    # first event of an immigrant: he/she enters the population later than sim. starting time
    lagToWaitingTime <- ifelse(isIMInitEvent, (as.numeric(calTime) - as.numeric(simHorizon[1]))/365.25,0) 
     #cat('\n-----\nID: ',id,'\n')
     #print(inp)
    
    ageInYears <- getAgeInYears(currAge)
     #cat('Age: ',ageInYears,' - CalTime: ',years(calTime),'-',months(calTime),'-',days(calTime),'\n')
    # Possible destination states
    possTr <- transitionMatrix[which(rownames(transitionMatrix) %in% currState),]  
    possTr <- possTr[which(possTr !=0)]
    nextEventMatrix <- matrix(0, ncol=2, nrow=length(possTr))   
    # How many years (along age scale) remain until `maxAge'? 
    ranMaxAge <- maxAge-ageInYears       
    # How many years (along cal. time scale) remain until simulation end?        
    ranMaxYear <-  (as.numeric(simHorizon[2]) - calTime)/365.25  
    # Identify the time range that should be considered. 
    ran <- min(ranMaxYear,ranMaxAge)   
    #cat('Ran: ',ran,' - ranMaxAge: ',ranMaxAge,' - ranMaxYear: ',ranMaxYear,'\n')
    ranAge <- c(ageInYears,ageInYears+ranMaxAge)
    ranYear <- chron(c(calTime, as.numeric(calTime)+ran*365.25))
    #cat('RanAge: ',ranAge,' - ranYear: ',ranYear,'\n')
    # Extract transition history of individual until current cal. time.
    historiesInd <- transitions[as.numeric(transitions[,'ID']) %in% id & as.numeric(transitions[,'transitionTime']) <= calTime,,drop=F]
    initPopInd <- initPop[as.numeric(initPop[,'ID']) %in% id,]
    birthTime <-  initPopInd['birthDate']
    initState <-  as.character(unlist(initPopInd['initState']))
    # Extract for each state the duration until transition. 
    # Here, we have to differ between states of which we do not know when they are entered (i.e., `initial states' of members
    # of the starting population and the states of migrants when they entered the country), and 
    # the states we know the `entering date' as well as the `leaving date' (if the state has been left). 
    if(as.numeric(birthTime) < as.numeric(simHorizon[1]) | id %in% as.numeric(immigrPop[,'ID'])) { 
      dur <- rbind(c(initState,NA),cbind(historiesInd[,'To'], historiesInd[,'transitionTime'])) 
      dur <- cbind(dur,c(diff(as.numeric(dur[,2])),0))
      colnames(dur) <- c('TransitionTo','AtTime','durUntil')
      dur[which(is.na(dur[,'AtTime'])),'durUntil'] <- NA
    } else {  # Individual is born during simulation.
      birthTime <-  initPop[as.numeric(initPop[,'ID']) %in% id, 'birthDate']
      dur <- rbind(c(initState,birthTime),cbind(historiesInd[,'To'], historiesInd[,'transitionTime'])) 
      dur <- cbind(dur,c(diff(as.numeric(dur[,2])),0))
      colnames(dur) <- c('TransitionTo','AtTime','durUntil')
    }
    # Compute for each possible destination state a waiting time.  
    for(i in 1:length(possTr)){
      tr <- possTr[i]
      destState <-  names(tr)
      cS <- unlist(strsplit(currState,'/'))
      dS <- unlist(strsplit(destState, '/'))
      # To determine the duration (time elapsed since last transition) that applies for the considered destination state,
      # we have to determine the duration since the last change in the covariate affected. 
      # For example, to specify the time being married, we have to determine the duration since (last) marriage. 
      covToCh <- which((cS==dS)==F)
      durSinceLastCovCh <- Inf  # For the transition to `dead' so far the time elapsed since the last transition does not play any role.  
      if(length(covToCh)==1){
        covHist <- do.call(rbind,sapply(dur[,'TransitionTo'],strsplit,split='/'))[,covToCh]
        idd <- which(covHist==cS[covToCh])
        if(length(idd)>1){
          if(F %in% (diff(idd)==1)){
            y <- rev(idd)[c(-1,diff(rev(idd)))==-1]
            idd <- rev(y)[c(diff(rev(y)),1)==1]
          }
        }
        durSinceLastCovCh <- sum(as.numeric(dur[idd,'durUntil']))       # If I do not how long an individual already is in a state: This gives NA.
        if(is.na(durSinceLastCovCh))
          durSinceLastCovCh <- 0 # Then assume the individual just entered that state. --> This is not a very sophisticated solution.
      }  
      if(length(covToCh)>1 & (!destState %in% absStates)){
        cat('Recognized a possible transition implying a change of two or more covariates.',
            'Concerning the derivation of the time being elapsed since the last transition this feature is not yet implemented.', 
            'Current State: ',currState,' -> Possible transition to ',destState,'\n') 
      }
      indRateFctDET <- function(x){
        res <- eval(do.call(tr, 
                            args=list(age=trunc(ageInYears)+x,calTime=trunc(1970.001+calTime/365.25)+x,duration=trunc(durSinceLastCovCh/365.25)+x)))
        return(res)
      }      
      ranAccuracyInDays <- (0:(trunc(ran*365.25)+1))/365.25
      detE <- indRateFctDET(ranAccuracyInDays)
      daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
      if (Inf %in% detE) {
        timeToNext <- daysToTrInYears
      } else {   
        u <- -log(1-runif(1)) 
        #cat('It: ',i,'--u: ',u,'\n')
        # Extract individual transition rate (depending on age, calendar time, and time elapsed)  
        indRateFct <- function(x){
          ageIn <- ageInYears+x
          calIn <- 1970.001+calTime/365.25+x
          durIn <- durSinceLastCovCh/365.25+x
          res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn,duration=durIn)))   
          if(TRUE %in% (res<0))
            stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
          #cat('x: ',x,' -- res', res,'\n')
          #cat('\n---\n')
          return(res)
        }
        if(sum(indRateFct(0:ran))==0){ # Rate funtion contains only zeros.
          intHaz <- 0
        } else {        
          # Integrated hazard at max. value
          intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
          if(inherits(intHaz, 'try-error')){          
            intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
          }
        }
        # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.      
        if(u<=intHaz){
          invHazFct <- function(x){
            #cat('x: ',x,'\n')
            try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
            #print(try.res)
            if(inherits(try.res, 'try-error')){  
              #cat('Seemingly, divergent intergral for ID ',id,
              # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
              #  ' for random number: ',u,'\n')  
              try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
            } 
            #cat('res: ',try.res,'\n-----\n')
            return(try.res)
          }  
          # Find random waiting time. 
          timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root      
        } else {
          timeToNext <- Inf
        }
      }
      nextEventMatrix[i,1] <- destState    
      nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
    }
    #print(nextEventMatrix)
    nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
    if(dim(nE)[1]>1)
      nE <- nE[1,,drop=F]
    if(nE[1,2]!=Inf){
      # Cal. time of next event of individual. (If there is one.)
      tt <- chron(as.numeric(calTime) + as.numeric(nE[1,2]) - ageInYears%%1, out.format=c(dates='d/m/year', times='h/m/s'))
      #print(tt)
      #cat(nE[1,1],'---',tt,'\n')  
      # Check whether next event implies school enrolment. If yes, adjust transition time to ensure that the individual
      # enters school at Sept. 1 in the year he/she turns seven.
      if(isSchoolEnrolment(currState,nE[1,1])){
        enYear <- years(tt)
        if(as.numeric(months(tt)) <= 9) {
          enDate <- chron(paste(enYear,dateSchoolEnrol,sep='/'), format=c(dates='y/m/d'), out.format=c(dates='d/m/year'))
        } else {
          enYear <- as.numeric(as.character(enYear))#+1
          enDate <- chron(paste(enYear,dateSchoolEnrol,sep='/'), format=c(dates='y/m/d'), out.format=c(dates='d/m/year'))
        }      
        diffToEn <- as.numeric(enDate-tt)
        nE[1,2] <- as.numeric(nE[1,2]) + diffToEn 
      }
      # Enqueue new event (if there is one).
      queue <<- rbind(queue, c(id, t.clock, currState, currAge, nE[1,1], nE[1,2]))
    }  
    #cat('\n----------\n')      
    return(nE)
  }
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # E. INITIALIZATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Compute next events for members of starting population
  cat('Initialization ... \n')
  IN <- data.frame(ID=initPop[,'ID'],currState=initPop[,'initState'],age=(simHorizon[1]-initPop[,'birthDate']),
                   calTime=rep(simHorizon[1],dim(initPop)[1]),stringsAsFactors=FALSE) 
  init <- apply(IN, 1, getNextStep)
  # If immigrants enter the population, compute next events for them.
  if(!is.null(immigrPop)){
    IM <- data.frame(ID=immigrPop[,'ID'], currState=immigrPop[,'immigrInitState'], age=immigrPop[,'immigrDate']-immigrPop[,'birthDate'],
                     calTime=as.numeric(immigrPop[,'immigrDate']),stringsAsFactors=FALSE)
    immigrInitPop <- immigrPop[,c('ID','birthDate','immigrInitState')]
    colnames(immigrInitPop)[3] <- 'initState'
    initPop <- rbind(initPop, immigrInitPop)
    imit <- apply(IM, 1, getNextStep, isIMInitEvent=T)
  } 
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # F. SIMULATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Run simulation either until queue is empty or until simulation horizon has been reached.
  cat('Simulation is running ... \n')
  currYear <- as.numeric(as.character(years(simHorizon[1])))
  while(dim(queue)[1]>0 & t.clock <= as.numeric(simHorizon[2])){
    
    # Sort queue according to soonest event to happen. 
    queue <- queue[order(as.numeric(queue[,'currTime'])+as.numeric(queue[,'timeToNextState'])),,drop=F]
    #print(chron(t.clock, format=c(dates='d/m/Y', times='h:m:s')))
    #print(dim(queue)[1])
    # Enqueue individual who has the soonest event to happen.
    indS <- queue[1,]  
    #print(indS)
    # Remove he/she from queue.
    queue <- queue[-1,,drop=F]
    # Set the global clock.
    t.clock <- as.numeric(indS['currTime']) +  as.numeric(indS['timeToNextState']) # in days
    cY <- as.numeric(as.character(years(t.clock))) 
    # If the global clock exceeds the simulation horizon, stop simulation.
    if(t.clock > as.numeric(simHorizon[2]))
      break 
    if(cY>currYear){
      cat('Year: ',cY,'\n')
      currYear <- cY
    }   
    # Age at current transition  
    age <- as.numeric(indS['currAge']) +  as.numeric(indS['timeToNextState'])
    # Register transition.
    transitions <- rbind(transitions, c(indS[c('ID','currState','nextState')], t.clock, age))
    # If current state is not an absorbent one, compute next event.
    if(!indS['nextState'] %in% absStates){
      # Current transition causes a newborn? If yes, add one to simulation population.
      if(isBirthEvent(indS['currState'],indS['nextState'])){
        addNewNewborn(t.clock)
      }  
      res <- getNextStep(c(indS[c('ID','nextState')], age, t.clock))
      #print(res)
    }
    #cat('\n-----------\n')
  }
  transitions <- transitions[order(as.numeric(transitions[,1])),,drop=F]
  
  if (nrow(transitions) == 0){
    
    transitionsOut <- data.frame(ID=initPop[,'ID'], From= rep(NA,nrow(initPop)), 
                                 To=rep(NA,nrow(initPop)), transitionTime = rep(NA,nrow(initPop)), 
                                 transitionAge = rep(NA,nrow(initPop)), stringsAsFactors = FALSE)
    cat('Simulation has finished.\n')
    cat('Beware that along the simulation horizon the individual/s considered do/es not experience any transition/s.\n')
    cat('------------------\n')
    
  } else {
    
    cat('Simulation has finished.\n------------------\n')
    
    # ----------------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    # G. GENERATE OUTPUT 
    # ----------------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    transitionsOut <- data.frame(ID=transitions[,'ID'], From=transitions[,'From'], To=transitions[,'To'],
                                 transitionTime = dates(chron(as.numeric(transitions[,'transitionTime']),
                                                              out.format=c(dates='d/m/year', times='h:m:s'))), 
                                 transitionAge = round(getAgeInYears(as.numeric(transitions[,'transitionAge'])),2),
                                 stringsAsFactors = FALSE)
  }

  pop <- merge(initPop, transitionsOut, all=T, by='ID')
  pop <- pop[order(as.numeric(pop[,1])),]  
  return(pop)
}

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# II. Execute microsimulation distributed (by executing as many single thread microsimulations in parallel as cores 
#     are available)
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
micSimParallel <- function(initPop, immigrPop=NULL, transitionMatrix, absStates=NULL, initStates=c(), initStatesProb=c(), 
                           maxAge=99, simHorizon, fertTr=c(), dateSchoolEnrol='09/01', cores=1, seeds=1254){
  
  cat('Starting at ');print(Sys.time())
  N <- dim(initPop)[1]
  M <- ifelse(is.null(immigrPop),0,dim(immigrPop)[1])
  # Split starting population and (if available) immigrant population accordinf to available cores     
  if(is.null(cores) | cores==1 | N+M<=20) {
    pop <- micSim(initPop, immigrPop, transitionMatrix, absStates, initStates, initStatesProb, 
                  maxAge, simHorizon, fertTr, dateSchoolEnrol)
  } else { 
    widthV <- max(trunc(N/cores), 10)
    widthW <- max(trunc(M/cores), 10)
    intV <- matrix(NA,ncol=2,nrow=cores)
    intW <- matrix(NA,ncol=2,nrow=cores)
    nI <- trunc(N/widthV)
    nIM <- trunc(M/widthW)
    ni <- 1
    for(i in 1:(nI-1)){
      intV[i,1] <- ni
      intV[i,2] <- ni+widthV-1
      ni <- ni+widthV
    }
    intV[nI,1] <- ni
    intV[nI,2] <- N
    ni <- 1
    if(nIM>1){
      for(i in 1:(nIM-1)){
        intW[i,1] <- ni
        intW[i,2] <- ni+widthW-1
        ni <- ni+widthW
      }
    }
    intW[nIM,1] <- ni
    intW[nIM,2] <- M
    initPopList <- list()
    immigrPopList <- list()  
    for(core in 1:cores){
      if(!is.na(intV[core,1])){
        initPopList[[core]] <- initPop[intV[core,1]:intV[core,2],]
      } else {
        initPopList[[core]] <- NA
      }
      if(!is.na(intW[core,1])){
        immigrPopList[[core]] <- immigrPop[intW[core,1]:intW[core,2],]                   
      } else {
        immigrPopList[[core]] <- NA
      }            
    }
    nL <- cores - sum(unlist((lapply(initPopList, is.na))))
    mL <- cores - sum(unlist((lapply(immigrPopList, is.na))))
    sfInit(parallel=T,cpus=cores,slaveOutfile=NULL)
    sfLibrary(chron)    
    sfExportAll()    
    sfClusterSetupRNGstream(seed=(rep(seeds,35)[1:length(cores)]))
    myPar <- function(itt){ 
      #cat('Starting thread: ',itt,'\n')   
      if(itt<=mL){        
        immigrPopL <- immigrPopList[[itt]]
      } else {
        immigrPopL <- NULL
      }     
      if (itt<=nL) {
        initPopL <- initPopList[[itt]]
      } else {
        initPopL <- NULL
        stop("\nCompared to the number of migrants, the starting population is too small to justify running a distributed simulation on several cores.")
      }
      popIt <- micSim(initPop=initPopL, immigrPop=immigrPopL, transitionMatrix=transitionMatrix, 
                      absStates=absStates, initStates=initStates, initStatesProb=initStatesProb, maxAge=maxAge, 
                      simHorizon=simHorizon, fertTr=fertTr, dateSchoolEnrol=dateSchoolEnrol)
      #cat('Thread: ',itt,' has stopped.\n') 
      return(popIt)      
    }
    pop <- sfLapply(1:max(nL,mL), myPar)   
    # create unique IDs for newborns 
    refID <- 0
    replaceID <- function(rr){
      pop[[i]][which(as.numeric(pop[[i]][,1])==rr[1]),1] <<- rr[2]
      return(NULL)
    }
    for(i in 1:length(pop)){
      allIDs <- c(initPopList[[i]]$ID, immigrPopList[[i]]$ID)
      exIDs <- unique(as.numeric(pop[[i]][,1]))
      repl <- setdiff(exIDs, allIDs)
      if(length(repl)>0) {
        newIDs <- cbind(repl,-(refID+(1:length(repl))))
        idch <- apply(newIDs,1,replaceID)
        refID <- max(abs(newIDs[,2]))
      }   
    }    
    pop <- do.call(rbind,pop)
    pop[as.numeric(pop[,1])<0,1]  <- abs(as.numeric(pop[as.numeric(pop[,1])<0,1]))+N+M
    pop <- pop[order(as.numeric(pop[,1])),]
    sfStop()
  } 
  cat('Stopped at ');print(Sys.time())
  return(pop)
}


