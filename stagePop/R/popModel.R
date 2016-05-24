#' popModel
#'
#' Run the core model.
#'
#' The default solver options are:
#'
#'   \code{list('DDEsolver'='deSolve','tol'=1e-7,'hbsize'=1e3,'method'=lsoda,'atol'=1e-7,'dt'=0.1)}
#'
#' but these may be changed by use of the \code{solverOptions}
#' parameter. Please see \code{\link{SolverOptions}} for details of this
#' parameter.
#'
#' @param numSpecies Number of species (integer)
#' @param numStages Number of life-stages in each species (vector)
#' @param numStrains Number of strains in each species (vector). Default is 1.
#' @param timeVec Vector of times the solution should be output out e.g. seq(1,10,0.1)
#' @param rateFunctions A list of rate functions to use. See \code{\link{RateFunctions}}
#' @param timeDependLoss A vector specifying TRUE/FALSE for each species. If a species has any time dependent per capita death or emigration rates, the entry in the vector must be TRUE. The default is TRUE for each species. 
#' @param timeDependDuration A vector specifying TRUE/FALSE for each species. If a species has any time dependent stage durations the entry is TRUE for that species. The default is FALSE for each species.
#' @param ICs is a list of matrices containing the initial conditions for every stage and strain of each species. These must be zero for all stages apart from the reproductive stage (usually the last stage). Each species has a matrix with the number of columns equal to the number of strains in that species and the number of rows equal to the number of stages in that species. E.g. for 2 species, the first with 2 strains and 3 stages, the second with 1 strain and 1 stage, then for zero starting conditions: ICs=list(matrix(0,ncol=2,nrow=3),matrix(0,ncol=1,nrow=1)). Due to the restrictions on initial conditions, it is recommended that more complicated initial conditions are defined through immigration rates in immigrationFunc
#' @param solverOptions Options for the DDE solver. A list containing 'DDEsolver' (can be 'deSolve' or 'PBS'), 'tol' (max error tolerance for DDE solver), 'hbsize' (history buffer size), 'method' (method for DDE solver), 'atol' (absolute tolerance (deSolve only)) and 'dt' (maximum initial timestep (PBS only)). Default is solverOptions=list(DDEsolver='PBS',tol=1e-7,hbsize=1e3,method='lsoda',atol=1e-7,dt=0.1)
#' @param checkForNegs If TRUE the function checkSolution is called and the solution for each variable, x, is checked for negative values that are greater in magnitude than ntol*max(x). If negative values occur then the solution is incorect and either the problem is incorrectly specified or the tolerances in the DDE solver need to be smaller. The default is TRUE.
#' @param ntol This controls the tolerance on the warning given for negative values when checkSolution is called. The default is 0.01 i.e. negative values whose magnitude is less than 1 percent of the max value of the variable are allowed.
#' @param plotFigs If TRUE, results will be automatically plotted during the model run. The default is TRUE.
#' @param speciesNames A vector of strings containing the name for each species. 
#' @param stageNames A list of n vectors (where n is the number of species) containing the names of each stage for each species. 
#' @param saveFig Choose to save the figure (TRUE or FALSE). Default is FALSE. 
#' @param figType Figure format can be 'eps', 'tiff' or 'png'. Default is 'eps' 
#' @param figName filepath to save figure to. Default is 'stagePopFig'
#' @param sumOverStrains If any of the species contain multiple strains then if this is TRUE the output is given as the sum over all the strains in the species. If this is FALSE then the time series for each strain will be in the output. Default is TRUE.
#' @param plotStrainsFig If any of the species contain multiple strains then if this is TRUE these will be plotted. Default is TRUE
#' @param saveStrainsFig If any of the species contain multiple strains then if this is TRUE the figures for the strains will be saveed
#' @param strainsFigType If any of the species contain multiple strains and if saveStrainsFig=TRUE then this is used to choose the type of file the figure is saved as (choose from 'eps', 'pdf', 'png' and 'tiff'). Default is 'eps'.
#' @param strainsFigName If any of the species contain multiple strains and if saveStrainsFig=TRUE then this is used to choose the name of file the figure is saved as. Default is created by paste('strainFig',SpeciesName[i]).

#' @return The model output is a matrix where rows are points in time and the columns are the state variables. These are named according to the species names and stage names supplied in inputs; the prefixes 'prob', 'dur' and 'dot' refer to the probability of survival through the stage, the duration of the stage and the rate of change of the variable. 'prob' type variables only appear if the per capita death (or emigration) rate is variable in time and 'dur' only appears if the stage duration is variable in time.

#' @examples
#' rateFuncs=list(
#'   reproFunc=function(x,time,species,strain){
#'    v=10*x$flies['adults',1]*exp(-x$flies['adults',1]/100)
#'    return(max(v,0))
#'  },
#'  deathFunc=function(stage,x,time,species,strain){
#'    a=c(0.05,0.1,0.1); return(a[stage])
#'  },
#'  durationFunc=function(stage,x,time,species,strain){
#'    a=c(5,10); return(a[stage])
#'  },
#'  immigrationFunc=function(stage,x,time,species,strain){
#'    v=0
#'    if (stage==3 & time<1){v=100}; return(v)},
#'  emigrationFunc=function(stage,x,time,species,strain){return(0)}
#' )
#'
#'modelOutput = popModel(
#'  numSpecies=1,
#'  numStages=3,
#'  ICs=list(matrix(0,nrow=3,ncol=1)),
#'  timeVec=seq(0,100,0.5),
#'  timeDependLoss=FALSE,
#'  timeDependDuration=FALSE,
#'  rateFunctions=rateFuncs,
#'  solverOptions=list(DDEsolver='PBS',tol=1e-4,hbsize=1e4,dt=0.01),
#'  stageNames=list(c('eggs','larvae','adults')),
#'  speciesNames=c('flies')
#'  )

#' @export

popModel=function(
  numSpecies,   
  numStages,
  numStrains=rep(1,numSpecies),
  timeVec,
  speciesNames,
  stageNames,
  rateFunctions=defaultRateFunctions,  
  ICs,
  timeDependLoss=rep(TRUE,numSpecies),
  timeDependDuration=rep(FALSE,numSpecies),
  solverOptions=list(),
  checkForNegs=TRUE,
  ntol=0.01,
  plotFigs=TRUE,
  saveFig=FALSE,
  figType='eps',
  figName='stagePopFig',
  sumOverStrains=TRUE,
  plotStrainsFig=TRUE,
  saveStrainsFig=FALSE,
  strainsFigType='eps',
  strainsFigName='strainFig'
){

  #Check Input Args-------------------------------------------

  if (!is.list(stageNames)){
    if (numSpecies>1){
      stop("stageNames should be a list containing vectors of the stage names for each species. E.g. for 2 species with one stage and three stages respectively this could have the form stageNames=list('adults',c('eggs','juvs','adults')).")
    }else{stageNames=list(stageNames)}}
  
  for (i in seq(1,numSpecies)){
    if (length(stageNames[[i]])!=numStages[i]) stop(paste("The number of stage names supplied for species",i," is not the same as the total number of stages for the species. If there is only one stage in a species simply put 'all' (for example) for this species."))}
  
  if (length(timeDependLoss)!=numSpecies) stop("timeDependLoss must specify whether the per capita death rate is changing with time for any stage in each species (e.g. for 3 species timeDependLoss=c(TRUE,FALSE,FALSE) is correct if just the first species has a time dependent per capita death rate). Note density dependent per capita death rates should be entered as TRUE.")

  if (length(timeDependDuration)!=numSpecies) stop("timeDependDuration must specify whether the stage duration is changing with time for any stage in each species (e.g. for 3 species timeDependDuration=c(TRUE,FALSE,FALSE) is correct if just the first species has a time dependent stage duration).")
  
  if (!is.null(speciesNames) & (length(speciesNames)!=numSpecies)) stop("The number of species names supplied is not the same as the number of species. Please change speciesNames or, since this is an optional argument, if you do not want to use the species naming feature then simply remove this input argument from your call to popModel (it is only used for plotting).")

  checkICs(ICs,numSpecies,numStages,numStrains)
  ICsNamed=addNamesToList(ICs,speciesNames,stageNames,numSpecies)
  
  modelStartTime=min(timeVec)

  if (!rateFuncCheck(ICsNamed,numSpecies,numStages,numStrains,timeDependDuration,rateFunctions,speciesNames,stageNames,modelStartTime)) stop('The model can not run until deathFunc, reproFunc, immigrationFunc, emigrationFunc, durationFunc (and develFunc if stage durations are time dependent) are correctly defined.')
 
#  tryCatch(force(figName),finally=print('Error: figName must be a string'))
#  tryCatch(force(figType),finally=print('Error: figType must be a string'))
#  tryCatch(force(strainsFigName),finally=print('Error: strainFigName must be a string'))
#  tryCatch(force(strainsFigType),finally=print('Error: strainFigType must be a string'))
  
  
  #---------------------------------------------------------------


  
  definePop=list(
    'numSpecies'=numSpecies,
    'numStages'=numStages,
    'numStrains'=numStrains,
    'timeDependLoss'=timeDependLoss,
    'timeDependDuration'=timeDependDuration,
    'speciesNames'=speciesNames,
    'stageNames'=stageNames)

  yinit=initConditions(definePop,ICsNamed,rateFunctions,modelStartTime,speciesNames,stageNames)
  
  
  divides=vecDivideFunc(definePop)
  Ndivs=matrix(divides[,1:2],nrow=numSpecies,ncol=2)
  Pdivs=matrix(divides[,3:4],nrow=numSpecies,ncol=2)
  Tdivs=matrix(divides[,5:6],nrow=numSpecies,ncol=2)
  dividesList=list(Ndivs,Pdivs,Tdivs)
  
  myparams=append(
    list('definePop'=definePop,'divs'=dividesList,'yinit'=yinit,'lengthTime'=max(timeVec),'startTime'=modelStartTime),
    rateFunctions)


  timeVec=timeVec-modelStartTime
  
                                        #DDE solver 
  defaultSolverOptions=list(DDEsolver='PBS',tol=1e-7,hbsize=1e3,method='lsoda',atol=1e-7,dt=0.1)
  solverOptionsFinal=replace(defaultSolverOptions,names(solverOptions),solverOptions)
  
                                        #Run model
  if (solverOptionsFinal$DDEsolver=='deSolve'){
                                        #deSolve dede is based on Adams and BDF formulae
    Out=dede(y = yinit, times = timeVec, func = derivDede, parms=myparams, atol = solverOptionsFinal$atol, rtol=solverOptionsFinal$tol,control=list(mxhist=solverOptionsFinal$hbsize),method=solverOptionsFinal$method)  
  } else if (solverOptionsFinal$DDEsolver=='PBS'){
                                        #PBS dde is based on Runge-Kutta formulae
    Out=dde(y = yinit, times = timeVec, func = derivPBS, parms=myparams, tol = solverOptionsFinal$tol, hbsize=solverOptionsFinal$hbsize,dt=solverOptionsFinal$dt)
    Out=data.matrix(Out)
  } else {
    stop(paste("Unknown solver",solverOptionsFinal$DDEsolver,"invoked. Valid options are 'PBS' or 'deSolve'. Aborting"))
  }

  Out[,1]=Out[,1]+modelStartTime
  allVarNames=namingVariables(numSpecies,numStages,numStrains,speciesNames,stageNames,timeDependLoss,timeDependDuration)
  colnames(Out)=c('time',allVarNames)

  if (plotStrainsFig & max(numStrains)>1){plotStrains(Out,numSpecies,numStages,numStrains,speciesNames,stageNames,strainsFigName,strainsFigType,saveStrainsFig)}


  if (checkForNegs){checkSolution(Out,numSpecies,numStages,numStrains,ntol)}

  summedOut=sumStrains(Out,numSpecies,numStages,numStrains,speciesNames,stageNames,timeDependLoss,timeDependDuration)


  if (plotFigs){genericPlot(summedOut,numSpecies,numStages,speciesNames=speciesNames,stageNames=stageNames,saveFig=saveFig,figType=figType,figName=figName)}

  
  if (sumOverStrains){
    modelOutput=summedOut
  }else{modelOutput=Out}
  
  
  return(modelOutput)
}
