## **** 1D profile CI bounds by NLOPT ****
calcBounds1D <- function(CIdlr,## <0
                         string1, string2, lowval, hival, maxparval, CIvar) {
  INFO <- list(fittedNames=blackbox.getOption("fittedNames"),
               fittedparamnbr=blackbox.getOption("fittedparamnbr"),
               rosglobal=blackbox.getOption("rosglobal"), ## FR->FR dangereux
               optimizers=blackbox.getOption("optimizers"),
               fitobject=blackbox.getOption("fitobject"),
               scalefactor=blackbox.getOption("scalefactor")
               )
  notinKgspace <- CIvar %w/o% INFO$fittedNames ## e.g. latt2Ns2 %w/o% canonicals is latt2Ns2
  locchull <- providefullhull(notinKgspace)[[1]] ## 'full' dimensional, vertices constraints and affectedConstraints
  #
  if (length(notinKgspace)==0L){ ## no extra composite hull
    if(is.na(maxparval)) message.redef("(!) From calcBounds1D(): maxparval is NA were a value is expected")
    maxpt <- INFO$rosglobal$par
    maxlogL <- INFO$rosglobal$value
  } else {
    ## ML in composite space hull. rosglobal$par may not be in the alternative hull, and a higher max may be found in this hull
    if( ! is.na(maxparval)) message.redef("(!) From calcBounds1D(): maxparval has a value were NA is expected")
    initpt <- INFO$fitobject$x[which.max(INFO$fitobject$fitted.values), ]
    initpt <- fromFONKtoanyspace(initpt, colnames(locchull$vertices))
    ui <- -locchull$a;ci <- -locchull$b    ##  within-hullness occurs for ui %*% x -ci => 0
    candidates <- generateInitpts(bestProfiledOut=initpt, vertices=locchull$vertices, ui=ui, ci=ci, hrepr=locchull$Hrep, fixedlist=NA, precision="default")
    ###################
    bestresu <- list(value=- Inf)
    for (candidate in candidates) {
      resu <- optimWrapper( ##purefn,
        initval=candidate$par, gr=NULL, ## not initpt, 14/10/2011
        chullformats=locchull,
        control=list( ## parscale is provided within optimWrapper
          fnscale=-1/INFO$scalefactor, trace=FALSE, maxit=10000))
      if (resu$value>bestresu$value) bestresu <- resu
    } ## end loop over candidate starting points
    maxpt <- bestresu$par
    maxpt <- fromFONKtoanyspace(maxpt, colnames(locchull$vertices))
    maxparval <- maxpt[CIvar]
    maxlogL <- tofKpredict.nohull(maxpt,fixedlist = NULL)
  }
  #
  hull_a <- locchull$a
  hull_b <- locchull$b    ##  within-hullness occurs for ui %*% x -ci => 0 ie a %*% x -b <= 0
  xnames <- colnames(locchull$vertices)
  parPos <- which(xnames==CIvar)
  seuillogL <- maxlogL+CIdlr
  eval_g0 <- function(x,a,b) {
    names(x) <- xnames
    logL <- tofKpredict.nohull(x,fixedlist=NULL)
    c(seuillogL-logL,  a %*% x-b) ## using max() appears fatal...
  } ## must be <0
  if ( "NLOPT_LD_MMA_for_CI" %in% INFO$optimizers) {
    SignedParPosIndic <- as.numeric(seq(length(xnames))==parPos) ## for eval_grad_f0
    eval_jac_g0 <- function( x, a, b ) { ## (jacobian:) gradients of constraints
      # one row for each constraint, one col for each dim of x
      # first row is -dlogL/dx, other rows are trivial.
      locfn <- function(x) {
        names(x) <- xnames
        return(- tofKpredict.nohull(x,fixedlist=NULL))
      }
      firstrow <- grad(func=locfn, x=x)
      rbind(firstrow,a)
    }
    locOptimizer <- "NLOPT_LD_MMA"
    xtol_rel <- 1e-4 ## -4 is default value; e-12 is never attained;
  } else{
    # NLOPT_LN_COBYLA seems to behave oddly. Even given a good starting point, it may look near the (lower) bound after a few iterations, optimize from there, and terminate at a worse solution than the initial point.
    eval_grad_f0 <- eval_jac_g0 <- NULL
    locOptimizer <- "NLOPT_LN_COBYLA"
    xtol_rel <- 1.0e-12
  }
  lb <- apply(locchull$vertices,2,min)
  ub <- apply(locchull$vertices,2,max)
  CIpoints <- matrix(nrow=0, ncol=INFO$fittedparamnbr)
  CIlo <- NA
  if (abs(lowval-maxparval)<abs(maxparval*1e-08)) { ## not an error ## lowval==maxparval if relative diff is <1e-16 in absolute value
    # CIlo <- NA ## already so
  } else {
    errmesbeg <- paste(string1, " for ", CIvar, " could not be computed", sep="")
    objfn <- function(x,a,b) {x[parPos]} ## pour la borne inf
    if ( locOptimizer == "NLOPT_LD_MMA") {
      eval_grad_f0 <- function(x,a,b) { SignedParPosIndic }
    }
    optr <- try(nloptr(x0=maxpt, eval_f=objfn, lb=lb, ub=ub,
                   eval_g_ineq=eval_g0, a=hull_a, b=hull_b,
                   eval_grad_f=eval_grad_f0, eval_jac_g_ineq=eval_jac_g0,
                   opts=list(algorithm=locOptimizer, xtol_rel=xtol_rel, maxeval=-1)))
    if(inherits(optr,"try-error")) {
      errmsg <- paste(errmesbeg," (maybe out of sampled range)", sep="")
      message.redef(errmsg)
    } else {
      solution <- optr$solution
      names(solution) <- xnames
      check_tol <- xtol_rel*max(abs(hull_a %*% abs(solution))) ## cannot be more stringent ?
      if(any(hull_a %*% solution -hull_b>check_tol)) {
        errmsg <- paste(errmesbeg," (optimization result does not satisfy constraints)", sep="")
        warning(errmsg)
      } else {
      CIpoint <- solution
      if (length(notinKgspace)>0L){ ## extra composite hull
        CIpoint <- tofullKrigingspace(CIpoint) ## will be used by predict.HLfit
      }
        logl <- purefn(CIpoint,testhull = FALSE)
        if( logl<seuillogL-1e-4) {
          errmsg <- paste(errmesbeg," (optimization result has suspiciously low likelihood)", sep="")
          warning(errmsg)
        } else if( logl>seuillogL+1e-4) {
          errmsg <- paste(errmesbeg," (maybe out of sampled range: optimization result has suspiciously high likelihood)", sep="")
          message.redef(errmsg)
        } else {
          CIlo <- solution[parPos]
          CIpoints <- rbind(CIpoints, CIpoint)
        }
      }
    }
  }
  CIup <- NA ## default values...
  if (abs(hival-maxparval)<abs(maxparval*1e-08)) { ## not an error
  } else {
    errmesbeg <- paste(string2, " for ", CIvar, " could not be computed", sep="")
    if ( locOptimizer == "NLOPT_LD_MMA") {
      SignedParPosIndic <- - SignedParPosIndic
    }
    objfn <- function(x,a,b) {-x[parPos]} ## pour la borne sup
    optr <- try(nloptr(x0=maxpt, eval_f=objfn, lb=lb, ub=ub,
                       eval_g_ineq=eval_g0, a=hull_a, b=hull_b,
                       eval_grad_f=eval_grad_f0, eval_jac_g_ineq=eval_jac_g0,
                       opts=list(algorithm=locOptimizer, xtol_rel=xtol_rel, maxeval=-1)))
    if(inherits(optr,"try-error")) {
      errmsg <- paste(errmesbeg," (maybe out of sampled range)", sep="")
      message.redef(errmsg)
    } else {
      solution <- optr$solution
      names(solution) <- xnames
      check_tol <- xtol_rel*max(abs(hull_a %*% abs(solution))) ## cannot be more stringent ?
      if(any(hull_a %*% solution -hull_b>check_tol)) {
        errmsg <- paste(errmesbeg," (optimization result does not satisfy constraints)", sep="")
        warning(errmsg)
      } else {
      CIpoint <- solution
      if (length(notinKgspace)>0L){ ## extra composite hull
        CIpoint <- tofullKrigingspace(CIpoint) ## will be used by predict.HLfit
      }
        logl <- purefn(CIpoint,testhull = FALSE)
        if( logl<seuillogL-1e-4) {
          errmsg <- paste(errmesbeg," (optimization result has suspiciously low likelihood)", sep="")
          warning(errmsg)
        } else if( logl>seuillogL+1e-4) {
          errmsg <- paste(errmesbeg," (maybe out of sampled range: optimization result has suspiciously high likelihood)", sep="")
          message.redef(errmsg)
        } else {
          CIup <- solution[parPos]
      CIpoints <- rbind(CIpoints, CIpoint)
        }
      }
    }
  }
  return(list(CIlo=CIlo, CIup=CIup, CIpoints=CIpoints))
} ## end def bounds1D
