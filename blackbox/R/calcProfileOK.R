# optim dans subspace, constraints in full space
# contrast with profileBySubHull which uses constraints in subspace (and may use NLOPT)
# cf comments below...
profileByFullHull <- function(fixedlist=NA, otherlist=NULL,
                          extrap=FALSE,
                          locchull=NULL) {
  INFO <- list(fittedNames=blackbox.getOption("fittedNames"),
               rosglobal=blackbox.getOption("rosglobal"), ## FR->FR dangereux
               ParameterNames=blackbox.getOption("ParameterNames"),
               paramnbr=blackbox.getOption("paramnbr"),
               optimizers=blackbox.getOption("optimizers"),
               Kgtotal=blackbox.getOption("hulls")$Kgtotal)
  ## Nb in metric units for metric plot.
  if ("Nb" %in% c(names(otherlist), names(fixedlist))) {
    message.redef(c("(!) From profile(): 'Nb' in function argument is highly suspect"))
    stop.redef()
  }
  notinKgspace <- names(fixedlist) %w/o% INFO$fittedNames
  checknames <- notinKgspace %w/o% c(INFO$ParameterNames, "latt2Ns2", "Nratio", "NactNfounderratio", "NfounderNancratio")
  if (length(checknames)>0) {
    message.redef(c("Error: incorrect members ", checknames, " of fixedlist argument to profile() function"))
    stop.redef()
  } ## ELSE:
  ## Check that there are values for the 'fixed' variables:
  if (any( ! is.numeric(as.numeric(fixedlist)))) {
    message.redef("Error: no value for some fixed variable in profile() function:")
    print(fixedlist)
    return(list(message="Error: no value for some fixed variable in profile() function."))
  }
  if(is.null(locchull)) locchull <- providefullhull(names(fixedlist))[[1]] ## this works also without any composite variable
  hull_a <- locchull$a
  hull_b <- locchull$b    ##  within-hullness occurs for ui %*% x -ci => 0 ie a %*% x -b <= 0
  lb <- apply(locchull$vertices,2,min)
  ub <- apply(locchull$vertices,2,max)
  #
  xnames <- names(lb)
  subnames <-  setdiff(xnames,names(fixedlist))
  if ( ! is.null(otherlist)) {
    maxpt <- tofullKrigingspace(fittedlist = otherlist,fixedlist = fixedlist)
  } else { ## tries to construct a point that satisfies the constraints from the global max
    maxpt <- fromFONKtoanyspace(FONKinput = INFO$rosglobal$par,outputnames = xnames)
    maxpt <- tofullKrigingspace(fittedlist = maxpt[subnames],fixedlist = fixedlist)
  }
  maxpt <- fromFONKtoanyspace(FONKinput = maxpt,outputnames = xnames)
  #
  if ( ! extrap ) {
    if ( any(hull_a %*% maxpt - hull_b>0) ) {
      subV <-  subHullWrapper(vertices=locchull$vertices,equality=fixedlist,precision="double")$vertices
      if (NROW(subV)>0L) {
        inCurspace <- TRUE
        initpt <- colMeans(subV)
        initpt <- tofullKrigingspace(initpt,fixedlist = fixedlist)
        initpt <- fromFONKtoanyspace(initpt,outputnames = xnames)
        candidates <- generateInitpts(bestProfiledOut=initpt, vertices=locchull$vertices,
                                      ui= - hull_a, ci= - hull_b, hrepr=locchull$Hrep, fixedlist=NA, precision="default")
        ###################
        bestresu <- list(value=- Inf)
        for (candidate in candidates) { ## runs over elements of a list, eahc of which is itself a list with a $par
          resu <- optimWrapper( ##purefn,
            initval=candidate$par, gr=NULL, ## not initpt, 14/10/2011
            chullformats=locchull,
            control=list( ## parscale is provided within optimWrapper
              fnscale=-1, trace=FALSE, maxit=10000))
          if (resu$value>bestresu$value) bestresu <- resu
        } ## end loop over candidate starting points
        maxpt <- bestresu$par
        maxpt <- fromFONKtoanyspace(maxpt, colnames(locchull$vertices))
      } else inCurspace <- FALSE
    } else inCurspace <- TRUE
    if (inCurspace) {
      maxpt <- maxpt[subnames]
      lb <- lb[subnames]
      ub <- ub[subnames]
      # optim dans subspace, constraints in full space
  calcfullx <- function(x) {
    names(x) <- subnames
    fx <- tofullKrigingspace(x, fixedlist)
    return(fromFONKtoanyspace(fx,outputnames = xnames))
  }
  objfn <- function(x,a,b) {
        names(x) <- subnames
        logL <- tofKpredict.nohull(x,fixedlist=fixedlist)
    return( - logL) ## minimiz
      }
  eval_g0 <- function(x,a,b) {
    fullx <- calcfullx(x)
    (a %*% fullx-b) ## avoids max() as in calcBounds1D()
  } ## must be <0
      if ( "NLOPT_LD_MMA" %in% INFO$optimizers) {
        ## HOWEVER, using numerical grad and jacobian is slow => NLOPT_LD_MMA is slow
        eval_jac_g0 <- function( x, a, b ) { ## (jacobian:) gradients of constraints
          # one row for each dim of fullx, one col for each dim of x
          jacfullx <- jacobian(func=calcfullx, x=x)
          a %*% jacfullx
        }
    eval_grad_f0 <- function(x,a,b) {
      locfn <- function(x) {objfn(x,a,b)}
      grad(func=locfn, x=x)
    }
    locOptimizer <- "NLOPT_LD_MMA"
        xtol_rel <- 1e-4 ## -4 is default value; e-12 is never attained;
  } else{
        ## less slow than (NLOPT_LD_MMA with num deriv) but slower than constrOptim
    eval_grad_f0 <- eval_jac_g0 <- NULL
    locOptimizer <- "NLOPT_LN_COBYLA"
        xtol_rel <- 1.0e-8
  }
      #
  optr <- try(nloptr(x0=maxpt, eval_f=objfn, lb=lb, ub=ub,
                     eval_g_ineq=eval_g0, a=hull_a, b=hull_b,
                     eval_grad_f=eval_grad_f0, eval_jac_g_ineq=eval_jac_g0,
                         opts=list(algorithm=locOptimizer, xtol_rel=xtol_rel, maxeval=-1)))
    } else { ## ! inCurspace & ! extrap
      resu <- list(full=rep(NA, INFO$paramnbr),
                   par=rep(NA, INFO$paramnbr),
                   value= NA,
                   message="safe exit from profile() because no subvertices found",
                   edgelevel=0L, ## hum
                   inKgspace=FALSE,
                   edge=list())
      return(resu)
    }
  } else { ## extrap=TRUE... not really tested
    optr <- try(nloptr(x0=maxpt, eval_f=objfn, lb=lb, ub=ub,
                       opts=list(algorithm="NLOPT_LN_BOBYQA", xtol_rel=1e-5, maxeval=-1)))
  }
  ## We reach this if extrap or inCurspace:
  solution <- optr$solution
  names(solution) <- subnames
  vkrig <- tofullKrigingspace(solution,fixedlist)
  canon <- canonize(vkrig)$canonVP ## completion/conversion to canonical
  if (extrap==TRUE || length(notinKgspace)>0L ) {
    inKgspace <- isPointInCHull(vkrig, constraints=INFO$Kgtotal[c("a", "b")])
    # this typically determines whether the relLik appears as "extrapolated" in a plot
  } else inKgspace <- inCurspace
  resu <- list(full=canon,   ## must be a suitable argument for tofullKrigingspace;
              par=solution,  ## in subhull
              value= - optr$objective, ## from minimiz...
              message="",
              edgelevel=0L, ## hum
               inKgspace=inKgspace,
               edge=list())
  return(resu) ## zut$canon is vector in canonical param space while zut$par (not always in return value) is vector of parameters profiled out
} ## end profile(...)
