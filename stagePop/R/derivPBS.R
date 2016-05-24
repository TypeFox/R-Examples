#' derivPBS
#'
#' Internal-use function which allows use of the PBSddesolv ODE
#' solver. \code{\link{popModel}}
#'
#' note PBS does not allow lagged variables to be computed at the
#' current timestep (hence the changes from the dede deriv function)
#'
#' @param t Time step
#' @param y state variables
#' @param parms Additional parameters passed in from the original invoker of the solver
#'
#' @return the new model state.
#'
#' @seealso \code{\link{derivDede}}
derivPBS=function(t,y,parms){

#  print('enter solver')
  
  #Recruitment function----------------------------------
  RFunc=function(reqdStage,reqdTime){
#  print('enter RFunc')
    if (reqdTime>0){
      
      if (reqdTime<t){
        xLagged=vecToMatList(pastvalue(reqdTime)[1:lenx],numSpecies,numStages,numStrains,speciesNames,stageNames)
      }else{xLagged=x}

      im=parms$immigrationFunc(reqdStage,xLagged,(reqdTime+startTime),species,strain)

      if (reqdStage==1){

        v=parms$reproFunc(xLagged,(reqdTime+startTime),species,strain)

      }else{
        
        if (timeDependLoss | timeDependDuration){
          if (reqdTime<t){
            PLagMat=matrix(pastvalue(reqdTime)[Pdivs[species,1]:Pdivs[species,2]],nrow=(num.stages-1),ncol=num.strains)
            PLag=PLagMat[(reqdStage-1),strain]
          }else{PLag=P[reqdStage-1]}
        }else{PLag=P[reqdStage-1]}

        
        if (timeDependDuration){

          if (reqdTime<t){
            dTauLagMat=matrix(pastgradient(reqdTime)[Tdivs[species,1]:Tdivs[species,2]],nrow=(num.stages-1),ncol=num.strains)
            devel=1-dTauLagMat[(reqdStage-1),strain]
            TauLag=matrix(pastvalue(reqdTime)[Tdivs[species,1]:Tdivs[species,2]],nrow=(num.stages-1),ncol=num.strains)[(reqdStage-1),strain]
          }else{
            TauLag=Tau[reqdStage-1]
            devel=1-dTau[reqdStage-1]
          }
        }else{
          devel=1
          TauLag=Tau[reqdStage-1]
        }
        
        reqdTime=reqdTime-TauLag
        reqdStage=reqdStage-1
        v=RFunc(reqdStage,reqdTime)*PLag*devel
      }

    }else{v=0;im=0}

    recruit=v+im
    return(recruit)
  }

  #---------------------------------------------------
  startTime=parms$startTime
  modelTime=t+startTime

  lengthTime=parms$lengthTime #number
  numSpecies=parms$definePop$numSpecies #number
  speciesNames=parms$definePop$speciesNames #vector
  numStages=parms$definePop$numStages #vector
  stageNames=parms$definePop$stageNames #list of vectors
  numStrains=parms$definePop$numStrains #vector

  yinit=parms$yinit #vector

  Ndivs=parms$divs[[1]] #vector of 2 indices:start and end of state variables
  Pdivs=parms$divs[[2]] #vector of 2 indices: start and end of survival probabilities
  Tdivs=parms$divs[[3]] #vector of 2 indices: start and end of stage durations
  
  #create list of named matrices (species) with named rows (stages) for putting into the rate functions
  lenx=max(Ndivs)
  x=vecToMatList(y[1:lenx],numSpecies,numStages,numStrains,speciesNames,stageNames)
  xinit=vecToMatList(yinit[1:lenx],numSpecies,numStages,numStrains,speciesNames,stageNames)
  
  alldN={}; alldP={}; alldTau={}

  for (species in seq(1,numSpecies)){
#    print(paste('species',species))
    
    num.stages=numStages[species]
    num.strains=numStrains[species]
    timeDependLoss=parms$definePop$timeDependLoss[species]
    timeDependDuration=parms$definePop$timeDependDuration[species]

    #create matrices for state variables (Nmat), survival probs (Pmat) and stage durations (Tmat)
    Nmat=matrix(y[Ndivs[species,1]:Ndivs[species,2]],nrow=num.stages,ncol=num.strains)
    if (timeDependLoss | timeDependDuration){
      Pmat=matrix(y[Pdivs[species,1]:Pdivs[species,2]],nrow=(num.stages-1),ncol=num.strains)}
    if (timeDependDuration){
      Tmat=matrix(y[Tdivs[species,1]:Tdivs[species,2]],nrow=(num.stages-1),ncol=num.strains)}

    for (strain in seq(1,num.strains)){
#    print(paste('strain',strain))

      N=Nmat[,strain]
      dN=0*N;      R=0*N;      M=0*N
      death=0*N
      for (i in seq(1,num.stages)){
        death[i]=parms$deathFunc(i,x,modelTime,species,strain)+parms$emigrationFunc(i,x,modelTime,species,strain)}

      if (num.stages>1){
        
      #compute Tau and P
        if (timeDependDuration){
          Tau=c(Tmat[,strain],lengthTime)
          dTau=Tau[1:(num.stages-1)]*0
          g=Tau*0; gLagged=Tau*0
        }else{
          Tau=0*N; dTau={}
          for (i in seq(1,num.stages-1)){Tau[i]=parms$durationFunc(i,x,modelTime,species,strain)}
          Tau[num.stages]=lengthTime}
        
        if (timeDependLoss | timeDependDuration){
          P=Pmat[,strain]
          dP=P*0
        }else{
          P=0*seq(1,num.stages-1); dP={} 
          for (i in seq(1,num.stages-1)){P[i]=exp(-(parms$deathFunc(i,x,modelTime,species,strain)+parms$emigrationFunc(i,x,modelTime,species,strain))*Tau[i])}}

        
        if (timeDependLoss | timeDependDuration){
          for (i in seq(1,num.stages-1)){
            tLagged=t-Tau[i]
            if (t==0 | tLagged<0){
              DLagged=parms$deathFunc(i,xinit,startTime,species,strain)+parms$emigrationFunc(i,xinit,startTime,species,strain)
              xLagged=xinit
              tLagged=0
            }else{
              if (tLagged<t){
                xLagged=vecToMatList(pastvalue(tLagged)[1:lenx],numSpecies,numStages,numStrains,speciesNames,stageNames)
              }else{xLagged=x}
              DLagged=parms$deathFunc(i,xLagged,(tLagged+startTime),species,strain)+parms$emigrationFunc(i,xLagged,(tLagged+startTime),species,strain)}
            
            if (timeDependDuration){
              gLagged[i]=parms$develFunc(i,xLagged,(tLagged+startTime),species,strain)
              if (gLagged[i]<=0)stop("please ensure the development rate is always greater than zero")
              g[i]=parms$develFunc(i,x,modelTime,species,strain)
              dTau[i]=1.0-(g[i]/max(gLagged[i],1e-12))
              dP[i]=P[i]*(DLagged*(1-dTau[i])-death[i])
            }else{dP[i]=P[i]*(DLagged-death[i])}}
        }


        for (i in seq(1,num.stages)){
          R[i]=RFunc(i,t)}
        
        for (i in seq(1,num.stages-1)){
          M[i]=R[i+1]-parms$immigrationFunc(i+1,x,modelTime,species,strain)}
        
      }else{
        dTau={}; dP={}
        R[1]=RFunc(1,t)
        M[1]=0
      }
            
      for (i in seq(1,num.stages)){
        dN[i]=R[i]-M[i]-death[i]*N[i]}
      
      alldN=c(alldN,dN)
      alldP=c(alldP,dP)
      alldTau=c(alldTau,dTau)    

      
    }#strain loop
  }#species loop
    
  dy=list(c(alldN,alldP,alldTau),alldN)
    
  return(dy)
}
