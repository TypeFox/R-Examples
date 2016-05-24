optimWrapper <- function(objectivefn=function(x) {tofKpredict.nohull(x, fixedlist=fixedlist)},
                       initval=NULL, chullformats=NULL,
                       maximum=T, oneDmethod=NA, gr, control, fixedlist=NULL, ...) { ## must return a point in 'full' kriging space
  ## this function performs maximization of 'fn' within the convex hull of 'chullformats''
  ## through a call to constrOptim or constrOptimR[ational].
  ## 'fn' itself does not need to check in-convex-hullness when called through optimWrapper
  ## 'initval' and 'fixedlist' are in 'any' space but together they must define a point in kriging space
  ## first tests whether 1D optimization or not
  ## fixedlist is NULL when the ML is first sought, but also (for LRT on composite variables)
  ## when a ML in composite space must be sought before the profile for the tested value is sought
  if (!is.null(initval)) {
    varnbr <- length(initval)
  } else {message.redef("(!) Error in optimWrapper(): no 'initval' argument")}
  ## objectivefn is the fn 'f' called in constrOptim...
  if (is.null(chullformats)) {
    stop.redef("(!) is.null(chullformats) at call to constrOptimM() in optimWrapper()")
  } ##ELSE
  chull.pts <- chullformats$vertices
  if (is.null(chull.pts)) {
    stop.redef("(!) is.null(chull.pts) at call to constrOptimM() in optimWrapper()")
  } ##ELSE
  ####    objectivefn <- function(x) {tofKpredict.nohull(x, fixedlist=fixedlist)} ## the simple way to pass fixedlist to the inner computations...
  lower <- apply(chull.pts, 2, min)
  upper <- apply(chull.pts, 2, max)
  parscale <- upper-lower
  control <- c(control, list(parscale=parscale/blackbox.getOption("scalefactor")))
  initinkg <- tofullKrigingspace(initval, fixedlist)
  initinfh <- fromFONKtoanyspace(initinkg, names(lower)) ## in 'full' space
  if (varnbr==1) { # uses lower, upper, maximum
    if (is.null(lower) | is.null(upper)) {message.redef("(!) Error in optimWrapper(): no 'lower' or 'upper' arguments")}
    if( ! is.na(oneDmethod) && oneDmethod=="nlminb") {## ideal for local optimization around initval
      ## FR->FR retester l'utilisation de nlminb ? not clear why not always used
      message.redef("(!) messy code in optimWrapper(): rethink.")
      stop.redef("(!) From OptimWrapper: probably called from modprofile, but this has not been chercked recently. Check 'signe' ETC. ")
      resu <- nlminb(initinfh, objectivefn,
                   ##   signe=-1, ## valid for modprofile only !! implies value=-...
                   lower=lower, upper=upper)
      names(resu$par) <- names(lower)
      resuinkg <- tofullKrigingspace(resu$par, fixedlist)
    } else { ## still 1D but no 1D method provided (silly code ?)
      ## an ugly patch on names(x) (optimize loses names but tofullKrigingspace needs them
      objectivefn <- function(x) {names(x) <- names(lower);tofKpredict.nohull(x, fixedlist=fixedlist)} ## the simple way to pass fixedlist to the inner computations...
      resu <- optimize(objectivefn, lower=lower, upper=upper, maximum=maximum)
      names(resu$maximum) <- names(lower)
      resuinkg <- tofullKrigingspace(resu$maximum, fixedlist)
    }
    if ((resuinkg-lower)/(upper-lower)<0.01 ## close to lower
        || (upper-resuinkg)/(upper-lower)<0.01 ## close to upper
    ) {edgelevel <- 1} else {edgelevel <- 0} ## this ad hoc code shows the ad hoc nature of the general method for detecting edgeness...
    return(list(par=resuinkg, value=resu$objective, edgelevel=edgelevel)) ## 2: maximization in a line=> two masking points
  } else { # NOT 1D
    if (ncol(chull.pts)!=(length(initval)+length(fixedlist))) { ## chull is in 'full' dimensional space !
      stop.redef("(!) From optimWrapper(): ncol(chull.pts)!=length(initval)+length(fixedlist)")
    } ##ELSE
    ## first do optim within simple bounds in all cases
    fittedparamnbr <- blackbox.getOption("fittedparamnbr")
    if (TRUE) { ## FR->FR NLOPT_LN_BOBYQA works even if some (upper==lower); 
      ## FR->FR ./. BUT test isPointInCHull will fail if we are trying to optimize on a facet of the hull. 
      parnames <- names(lower)
      nloptr_objectivefn <- function(x) {names(x) <- parnames; -tofKpredict.nohull(x, fixedlist=fixedlist)}
      resu <- nloptr(x0=unlist(initinfh),eval_f=nloptr_objectivefn,lb=lower,ub=upper,
                     opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=-1))
      resu$par <- structure(resu$solution,names=parnames)
      resu$value <- - resu$objective
    } else resu <- optim(initinfh, objectivefn, gr=objectivefn.grad, method="L-BFGS-B", control=control, lower=lower, upper=upper)
    ## resu$par is in 'full' Kriging variable *space* (but not necessary within hull of kriged points) since initinfh is
    ## check the OPTIM result
    if( ! (isPointInCHull(resu$par, constraints=chullformats[c("a", "b")]))) { ## if simple optim result not in convex hull
      ## then we will optimize in the convex hull using constroptim(R). But we need a good starting point
      objectivefn.grad <- function(x) {grad(func=objectivefn, x=x)} ## no need to specify fixedlist here
      notinKgspace <- names(chullformats$vertices) %w/o% blackbox.getOption("fittedNames")
      precision <- switch(blackbox.getOption("redundant.mode"),
                        "noElim"="no.elim",
                        "alwaysRational"="rational",
                        "alwaysDouble"="double",
                        "defaultPrecision"="double")
      class(resu) <- "try-error"
      errorThreshold <- 1e-14
      if (precision=="rational") {
        ## it seems that we still need to consider later errors in floating point operations, which makes rational not so useful here
        ui <- qneg(chullformats$a);ci <- qneg(chullformats$b)
      } else {
        ui <- -chullformats$a;ci <- -chullformats$b
      }    ##  within-hullness occurs for ui %*% x -ci => 0
      while(inherits(resu,"try-error") && errorThreshold<1e-06) {
        errorThreshold <- errorThreshold*10
        check1 <- - isPointInCHull(point=initinfh, constraints=chullformats[c("a", "b")], exactArithmetic=(precision=="rational"), returnValue="vector")
        pbext <- which(check1<0)
        pblim <- (which(check1<errorThreshold) %w/o% pbext) ## ==0 unsurprizingly fails in some cases; 1e-14 not large enough.
        if (length(c(pbext, pblim))>0) { ## initinfh outside or on edge of hull
          ## SAFE procedure to construct a good init point for constrOptim, see C source for more explanations
          ## We first take the minimum number of points (dim+1) necess for a hull of dim fittedparamnbr,
          ## then add points if its actual dimension is lower. Actual dimension is determined from the presence of equality constraints in the H repr
          ##
          hull.liks <- apply(chull.pts, 1, tofKpredict.nohull, fixedlist=NULL)
          orderedhull <- chull.pts[order(hull.liks, decreasing=T), ]
          minimalhull <- orderedhull[seq(fittedparamnbr), ] ## one point will be added in the loop before any test
          minimalVrepr <- cbind(0, 1, minimalhull)
          nr <- length(hull.liks)
          it <- fittedparamnbr
          diagnostic <- 1
          if (precision=="rational") {## because scdd can fail in floating point
            minimalVrepr <- d2q(minimalVrepr)
            while(it<nr && diagnostic>0) {
              it <- it+1
              minimalVrepr <- rbind(minimalVrepr, d2q(c(0, 1, orderedhull[it, ])))
              HH <- scdd(minimalVrepr, representation="V")
              diagnostic <- sum(q2d(HH$output[, 1]))
            }
            minimalVrepr <- q2d(minimalVrepr)
          } else { ## not rational
            while(it<nr && diagnostic>0) {
              it <- it+1
              minimalVrepr <- rbind(minimalVrepr, c(0, 1, orderedhull[it, ]))
              HH <- try(scdd( minimalVrepr, representation="V"), silent=TRUE)
              if (inherits(HH,"try-error")) {
                message.redef("Operation failed. If due to floating-point error, then use 'RArguments=alwaysRational'.")
                stop.redef()
              }
              diagnostic <- sum(HH$output[, 1])
            }
          }
          if (it==nr && diagnostic>0) {
            stop.redef("(!) From optimWrapper: Suspicious it>nr")
          } ## otherwise we have reconstructed a fittedparamnbr-dimensional hull
          ##
          ## we now try to find a point close to initinfh but within the hull in the direction given by the 'mean' of the hull
          bordel <- colMeans(minimalVrepr)[-(1:2)]
          ## we determine the numerical threshold value 'seuil'
          if (length(pbext)>0) { ## initinfh outside of hull
            ## then use appropriate isPointInCHull returnValue to find the right threshold
            check2 <- - isPointInCHull(point=bordel, constraints=chullformats[c("a", "b")], exactArithmetic=(precision=="rational"), returnValue="vector")
            zut <- check2/(check2-check1)
            seuil <- min(zut[pbext]) ## threshold value where all ui %*% (.) -ci become nonnegative
          } else if (length(pblim)>0) { ## initinfh exactly on edge of hull
            seuil <- 1 ## should be enough (modulo floating point errors...)
          }
        } else { ## length(c(pbext, pblim))=0
          seuil <- 1
          bordel <- initinfh
        }
        class(resu) <- "try-error"
        globalshrinkingfactor <- 1
        itershrinkingfactor <- min(1e12, max(c(1e3, abs(bordel-initinfh)))) ## in [1e3, 1e12] as bordel-initinfh increases from 0 to Inf
        itershrinkingfactor <- 1-1/itershrinkingfactor ## from [0.999, 1-1e-12] as bordel-initinfh increases from 0 to Inf
        while( inherits(resu,"try-error")&& (1-globalshrinkingfactor)<10*(1-itershrinkingfactor) ) {
          globalshrinkingfactor <- itershrinkingfactor*globalshrinkingfactor ## slighly decreases
          eps <- 1-globalshrinkingfactor*seuil ## 1 - something that decreases from 1 => eps increases from 0
          reinit <- initinfh+(bordel-initinfh)*eps ## moving a bit inside the hull, EXCEPT if bordel=initinfh
          if (precision=="rational") {
            ## rational here not so useful as constrOptimR()->optim()->(fun)->'R'() not fully rational
            resu <- try(constrOptimR(reinit, objectivefn, grad=objectivefn.grad, ui=ui, ci=ci , mu=1e-08, ## a low mu appear important
                                   method = "BFGS",
                                   control = control)
            )
          } else {
            resu <- try(constrOptim(reinit, objectivefn, grad=objectivefn.grad, ui=ui, ci=ci , mu=1e-08, ## a low mu appear important
                                  method = "BFGS",
                                  control = control)
            )
          }    ##  within-hullness occurs for ui %*% x -ci => 0
          if(inherits(resu,"try-error")) {
            message.redef("(*) From optimWrapper:() constrOptim call possibly found")
            message.redef("    that the initial value was not in the convex hull. Trying another value...")
            message.redef("   Diagnostic information:")
            message.redef(paste("   'reinit': ", paste(reinit)))
            message.redef(paste("   'globalshrinkingfactor': ", globalshrinkingfactor))
            message.redef(paste("   'precision': ", precision))
            message.redef(paste("   'min(check1)': ", min(check1)))
            message.redef("   'fixedlist':")
            message.redef(fixedlist)
          }
        } ## fin inner while (shrinkingfactor)
      } ## fin outer while (errorthreshold)
      if (inherits(resu,"try-error")) {
        stop.redef("(!) optimWrapper() failed to find a good starting point.")
      } ## else we have a valid 'resu'
    } ##else optim result is within hull.
    fullerhull <- matchVertCons(chullformats) ## ideally optimWrapper was called with chullformats containing affectedConstraints
    edgecheck <- handlesPointAtEdge(point=resu$par, fixedlist=fixedlist, fullerhull=fullerhull)
    if (!is.null(edgecheck$edge)) {
      locedge <- t(apply(edgecheck$edge, 1, tofullKrigingspace, fixedlist=fixedlist))
      colnames(locedge) <- blackbox.getOption("fittedNames")  ## FR->FR 10/2015 comprends plus F*ON*K...
      addtoedges(, locedge)
    }
    if (!is.null(edgecheck$resu)) {
      if (edgecheck$resu$value>resu$value) {
        resu <- edgecheck$resu
      } ## else resu near the edge remains better that the best point on the edge => nothing to do
    }
    resuinkg <- tofullKrigingspace(resu$par, fixedlist)
    return(list(par=resuinkg, value=resu$value, convergence=resu$convergence, message=resu$message, edgelevel=edgecheck$edgelevel, edge=list(edgecheck$edge)))
  } ## endif not 1D
} ## end def optimWrapper
