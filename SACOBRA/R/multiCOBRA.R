#######################################################################
# multiCOBRA
#
# Samineh Bagheri, Wolfgang Konen
# Jan-May 2015
# Cologne University of Applied Sciences
# 
########################################################################


######################################################################################
# multiCOBRA
#
#' Perform multiple COBRA runs
#' 
#' Perform multiple COBRA runs. Each run starts with a different seed so that a different
#' start point, a different initial design and different random restarts are choosen.
#' 
#' Side effect: An error plot showing each run and the mean and median of all runs (see 
#' \code{\link{multiRunPlot}}). The results (\code{dfAll} and others) are saved to
#' \code{<fName>.Rdata}.
#' 
#'  @param fn         objective function that is to be minimized, should return a vector of the 
#'                    objective function value and the constraint values
#'  @param lower      lower bound of search space
#'  @param upper      upper bound of search space
#'  @param nrun       [10] number of runs
#'  @param feval      [200] function evaluations per run
#'  @param funcName   ["GXX"] name of the problem
#'  @param fName      file name .Rdata where the results (\code{dfAll} and others) are saved
#'                    (only if saveRdata==TRUE)
#'  @param path       [NULL] optional path 
#'  @param cobra      [NULL] list with COBRA settings. If NULL and for elements not present in
#'                    this list, the defaults from \code{\link{cobraInit}} are used. 
#'  @param optim      [NULL] the true optimum (or best known value) of the problem
#'  @param target     [0.05] a single run meets the target, if the final error is smaller than \code{target}
#'  @param ylim       the y limits
#'  @param plotPDF    [FALSE] if TRUE, plot not only to current graphics device but 
#'                    to \code{<fName>.pdf} as well
#'  @param saveRdata  [FALSE] if TRUE, save results (dfAll,optim,target,fName,funcName) on fName                  
#'  @param startSeed  [41]
#'
#'  @return \code{mres}, a list containing
#'      \item{\code{cobra}}{ the settings and results from \strong{last} run }
#'      \item{\code{dfAll}}{ a data frame with a result summary for \strong{all} runs (see below) }
#'      \item{\code{z}}{ a vector containing for each run the ever-best feasible objective value}
#'      \item{\code{z2}}{ a data frame containing for each run the minimum error 
#'                        (if \code{optim} is available)}
#'
#'   The data frame \code{dfAll} contains one row per iteration with columns (among others)
#'   \describe{
#'      \item{ffc}{ fitness function calls (i.e. the iterations \code{cobra$iter}) }
#'      \item{fitVal}{ true fitness function value }
#'      \item{fitSur}{ surrogate fitness function value }
#'      \item{feas}{  is current iterate feasible on the true constraints?}
#'      \item{feval}{ number of evaluations of the internal optimizer on the surrogate functions
#'                  (\code{NA} if it is a repairInfeasible-step)  }
#'      \item{XI}{  the DRC element used in the current iteration }
#'      \item{everBestFeas}{  the ever-best feasible fitness function value}
#'      \item{run}{ the number of the current run }
#'      \item{X1,X2,...}{ the solution in (original) input space }
#'   }
#' 
#' @examples
#' ## solve G11 problem nrun times  and plot the results of all nrun runs
#' nrun=4
#' feval=30
#' 
#' ## Defining the constraint problem: G11
#' fn <- function(x) {
#'   y<-x[1]*x[1]+((x[2]-1)^2)
#'   y<-as.numeric(y)
#'   
#'   g1 <- as.numeric(+(x[2] - x[1]^2))
#'   
#'   return(c(objective=y, g1=g1))
#' }
#' funcName="G11"
#' d=2
#' nConstraints<-1
#' lower<-c(-1,-1) 
#' upper<-c(+1,+1) 
#' 
#' ## Initializing and running cobra
#' cobra <- cobraInit(xStart=c(0,0), fn=fn, fName=funcName, lower=lower, upper=upper,
#'                    nConstraints=nConstraints,  feval=feval,
#'                    initDesPoints=3*d, DOSAC=1,cobraSeed=1)
#' 
#' mres <- multiCOBRA(fn,lower,upper,nrun=nrun,feval=feval,optim=0.75
#'                   ,cobra=cobra,funcName=funcName
#'                   ,ylim=c(1e-12,1e-0),plotPDF=FALSE,startSeed=42)
#' 
#' #Two true solutions at x1* = c(-sqrt(0.5),0.5) and x2* = c(+sqrt(0.5),0.5)
#' #where the true optimum is f(x1*) = f(x2*) = -0.75
#' print(getXbest(mres$cobra))
#' print(getFbest(mres$cobra))
#' print(mres$z2)
#' 
#' 
#' @seealso   \code{\link{multiRunPlot}}, \code{\link{cobraPhaseII}}
#' @author Wolfgang Konen, Samineh Bagheri, Cologne Univeristy of Applied Sciences
#' @export
#'     
multiCOBRA <- function(fn,lower,upper,nrun=10,feval=200
                      ,funcName="GXX",fName=paste0("mult-",funcName,".Rdata"),path=NULL
                      ,cobra=NULL,optim=NULL,target=0.05, saveRdata=FALSE
                      ,ylim=c(1e-05,1e+04),plotPDF=FALSE, startSeed=41)
{
  dimension <- d<- length(lower)
  nConstraints <- m<- length(fn(lower))-1
  if (is.null(optim)) {
    optimumAvail=FALSE
    optim=0.0
  } else {
    optimumAvail=TRUE
  }
  
  runG <- function(cobraSeed){
    
    set.seed(cobraSeed+10000)
    xStart<-runif(d,min=lower,max=upper)
    cobraRun <- cobraInit(xStart=xStart, 
                       fn=fn, 
                       fName=fName,
                       lower=lower, # lower bound constraints
                       upper=upper, # upper bound constraints
                       nConstraints=nConstraints,  # number of inequality constraints
                       feval=feval,       # maximum number of function evaluations
                       initDesPoints=3*d, 
                       cobraSeed=cobraSeed
                       )
    
    # if cobra is not NULL, set remaining elements of cobraRun from cobra
    cobraRun <- setOpts(cobra,cobraRun)       # setOpts def'd in defaultSAC.R
    
    cobraRun <- startCobra(cobraRun)
    return(cobraRun)
  } # runG
  
  
  ptm <- proc.time();
  dfAll = NULL
  for (run in 1:nrun) {
    cobraSeed = startSeed+run 
    cobraRun <- runG(cobraSeed) 
    df = data.frame( ffc=cobraRun$df$iter
                     ,fitVal=cobraRun$df$y
                     ,fitSur=cobraRun$df$predY
                     ,feas=cobraRun$df$feasible
                     ,feval=cobraRun$df$FEval
                     ,XI=c(rep(NA,cobraRun$initDesPoints),cobraRun$df2$XI)
                     ,everBestFeas=cobraRun$df$Best
                     ,run=rep(run,nrow(cobraRun$df))
                     ,data.frame(cobraRun$A[,],row.names=NULL)  # 'row.names' needed to suppress a warning
    )
    # browser()
    df$everBestFeas[1:cobraRun$initDesPoints]=NA
    #df = tail(df,niter)
    dfAll = rbind(dfAll,df)
  }
  cat(paste("Proc time for all runs on ",funcName,": ",sep=""),(proc.time()-ptm)[1],"\n");
  
  if (saveRdata) {
    filename=fName
    if (!is.null(path)) filename <- paste(path,fName,sep="/")
    save(dfAll,optim,target,fName,funcName, file=filename)
  }
  
  z=multiRunPlot(dfAll,optim, target=target
                 ,fName=fName,legendWhere="bottomleft"
                 ,main=paste(funcName ,"Problem")
                 ,xlim=c(3*d,max(dfAll$ffc))
                 ,ylim=ylim
  )
  
  
  if (plotPDF)
    z=multiRunPlot(dfAll,optim, target=target
                   ,fName=fName, plotPDF=TRUE, legendWhere="bottomleft"
                   ,main=paste(funcName ,"Problem")
                   ,xlim=c(3*d,max(dfAll$ffc))
                   ,ylim=ylim 
    )
  
  z2 = stats::aggregate(dfAll$everBestFeas,list(dfAll$run),min,na.rm=T)
  z2[,2]=z2[,2]-optim
  names(z2) <- c("run","error")
  
  mres <- list(cobra=cobraRun
               ,dfAll=dfAll
               ,z=z
               ,z2=z2
               )
  
  return(mres)
} # multiCOBRA

