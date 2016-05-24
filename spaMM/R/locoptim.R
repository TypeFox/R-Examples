## this function 'optimize' and by default minimizes its function
##    it wraps optim, optimize and nlminb
##    and performs by default a grid search before calling there
## It can be turned back to maximization  (minimization used only in confint 29/08/2014)
## the objective function 'objfn' is by default the objective function HLCor.obj (currently not by name bc do.call(optimize,...name)) fails)
##   the first arg of this fn must be ranefParsVec
## with respect to the variable in the LIST init.optim
## within bounds given by LowUp
## by first evaluating the objfn in selected points before running optimize (1D case) or 
## nlminb/one of the methods in optim depending on the 0ptimizer argument
## anyOptim.args contains arguments common to all optimizers, and arguments for the objective function
## optimizers.args contains arguments specific to the given optimizer
locoptim <- function(init.optim,LowUp,anyObjfnCall.args,trace=list(file=NULL,append=T),Optimizer="L-BFGS-B",optimizers.args,maxcorners=2^11,objfn=HLCor.obj,maximize=FALSE) {
  processedHL1 <- getProcessed(anyObjfnCall.args$processed,"HL[1]",from=1L) ## there's also HLmethod in processed<[[]]>$callargs
  initvec <- unlist(init.optim)
  if ( ! is.null(anyObjfnCall.args$ranefParsVec) ) stop("! is.null(anyObjfnCall.args$ranefParsVec)") ## catch programming errors
  optimizerObjfnCall.args <- notoptimObjfnCall.args <- anyObjfnCall.args
  notoptimObjfnCall.args$ranefParsVec <- initvec ## logscale, inherits names from init.optim
  lower <- unlist(LowUp$lower); upper <- unlist(LowUp$upper)
  if (length(initvec)>1) { ## FR->FR voir si on peut dire >0 car ca evite optimze failement bloque par un max local 
    ## old version before 2014/09/03
    parscale <- (upper-lower) ## unlist because some list elements may have length >1 (eg if rho has several     ###### basic unrefined fit for initvec...
    if (is.character(trace$file)) write("## Call for initial values",file=trace$file,append=T)   
    if (maxcorners) {      
      init.obj <- do.call(objfn,notoptimObjfnCall.args)
      ###### look in the corners 
      ## L-BFGS-B tends to find the local max closest to the initial point. Search for good initial point:
      byvar <- t(rbind(lower,upper)) 
      byvar <- 0.999 * byvar + 0.001 *rowMeans(byvar)
      grillelist <- list()
      gridSteps <- floor(35^(1.05/length(initvec))) ## 6 for 2 pars,  3 for 3 pars, then 2  2  1  1 
      gridSteps <- max(2,gridSteps)
      if (  setequal(sort(names(lower)),c("trNu","trRho")) ) {
        ## special ad hoc case, transect trying to catch the likeliy trRho,trNu ridge
        trRhoSeq <- seq(byvar["trRho",1L],byvar["trRho",2L],length.out=36)
        trNuSeq <- seq(byvar["trNu",2L],byvar["trNu",1L],length.out=36)
        corners <- cbind(trRho=trRhoSeq,trNu=trNuSeq)
      } else { ## general case
        for(name in rownames(byvar)) {grillelist[[name]] <- seq(byvar[name,1L],byvar[name,2L],length.out=gridSteps)}
        if ( (2^length(grillelist)) <= maxcorners ) {
          corners <- expand.grid(grillelist)
        } else { ## otherwise we cannot dream of allocating, say, 2^36 values...
          corners <- replicate(maxcorners,unlist(lapply(grillelist,sample,size=1)))
          corners <- t(corners)
        }
      }
      ## ranefParsVec corresponds to ranPars but as vector, not as list
      ## uses HLCor.obj because of the return value...
      if (is.character(trace$file)) write("## Calls for 'corners'",file=trace$file,append=T)  
      if (FALSE) { ## version for debug can use info attribute
        scorners <- split(corners,rownames(corners))
        ll <- lapply(scorners,function(v) {
          notoptimObjfnCall.args$ranefParsVec <- v
          ## FR->FR il serait mieux de pouvoir ne supprimer que certain warnings
          bla <- suppressWarnings(do.call(objfn,notoptimObjfnCall.args)) 
          c(bla,attr(bla,"info"))
        })
        corners.objWithAttr <- t(sapply(ll,function(v) {c(v,info=attr(v,"info"))}))
        corners.obj <- corners.obj[,1]
      } else {
        corners.obj <- apply(corners,1,function(v) {
          notoptimObjfnCall.args$ranefParsVec <- v
          ## FR->FR il serait mieux de pouvoir ne supprimer que certain warnings
          suppressWarnings(do.call(objfn,notoptimObjfnCall.args)) 
        })
      }
      if (maximize) {
        if (max(corners.obj) > init.obj) {
          initvec <- corners[which.max(corners.obj),] 
        }
      } else {
        if (min(corners.obj) < init.obj) {
          initvec <- corners[which.min(corners.obj),] 
        }
      }
    }
    ## Search for good initial point done. 
    if (is.character(trace$file)) write("## Optimization call",file=trace$file,append=T)   
    if (Optimizer=="nlminb") {
      ## nlminb code
      nlminbArgs <- optimizerObjfnCall.args ## must include $skeleton avec son attr RHOMAX
      if (maximize) {
        nlminbArgs$objective <- function(...) {-objfn(...)}
      } else nlminbArgs$objective <- function(...) {objfn(...)}
      nlminbArgs$start <- initvec    
      nlminbArgs$lower <- lower
      nlminbArgs$upper <- upper
      nlminbArgs$scale <- 1/(upper - lower)
      nlminbArgs$control <- optimizers.args$nlminb$control
      optr <- do.call(nlminb,nlminbArgs) ## optimize HLCor.obj.value (=objective=p_bv by default, except for SEM)
    } else {
      control <- list(parscale=parscale,factr=1e9) ## factr was the stricter 1e8 up to 23/01/13
      if (maximize) control$fnscale <- -1     
      control[names(optimizers.args$optim$control)] <- optimizers.args$optim$control ## ...which may be overwritten 
      ## lower and upper were missing 7/11/2014!
      anyBaseOptim.args <- c(optimizerObjfnCall.args,list(par=initvec,lower=lower,upper=upper,control=control,method=Optimizer)) 
      optr <- do.call("optim",c(anyBaseOptim.args,list(fn=objfn))) ## optimize HLCor.obj.value (=p_bv by default, except for SEM)
    }
    optPars <- relist(optr$par,init.optim)
    attr(optPars,"optr") <- optr
    ## HLCor.args$ranPars[names(optPars)] <- optPars 
  } else if (length(initvec)==1L) { ## one-dimensional optimization
    byvar <- c(lower,upper) ## HLCor expects trLambda...
    byvar <- 0.999 * byvar + 0.001 *mean(byvar)
    corners <- seq(byvar[1],byvar[2],length.out=11)
    init.obj <- do.call(objfn,notoptimObjfnCall.args)
    if (is.character(trace$file)) write("## Calls for 'corners'",file=trace$file,append=T)   
    corners.obj <- sapply(corners,function(v) {
      notoptimObjfnCall.args$ranefParsVec <- v
      do.call(objfn,notoptimObjfnCall.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
    })
    if (maximize) {
      if (max(corners.obj) > init.obj) {
        initvec <- corners[which.max(corners.obj)] 
      }
    } else {
      if (min(corners.obj) < init.obj) {
        initvec <- corners[which.min(corners.obj)] 
      }
    }
    ## optimize uses only an interval, no initial value (!!)
    intervlower <- max(lower,corners[corners<initvec])
    intervupper <- min(upper,corners[corners>initvec])
    anyOptimize.args <- optimizerObjfnCall.args
    if (maximize) anyOptimize.args$maximum <- TRUE
    anyOptimize.args$interval <-c(intervlower,intervupper)
    anyOptimize.args$tol <- optimizers.args$optimize$tol
    if (is.character(trace$file)) write("## Optimization call",file=trace$file,append=T)   
    optr <- do.call(optimize,c(anyOptimize.args,list(f=objfn)))  ## optimize p_bv
    if(maximize) {
      optPars <- relist(optr$maximum,init.optim)
    } else optPars <- relist(optr$minimum,init.optim)
    ## HLCor.args$ranPars[names(optPars)] <- optPars 
  } else { ## nothing to optimize
    ## HLCor.args$ranPars should already contain all it needs
    optPars <- NULL
  }
  return(optPars)
} ## end def locoptim
