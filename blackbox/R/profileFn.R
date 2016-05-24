profileBySubHull <- function(fixedlist=NA, otherlist=NULL, extrap=F, locchull=NULL,
                  templateinkg=blackbox.getOption("rosglobal")$par, ##but another $par can be passed if this one is not in locchull
                  subHull_V_and_H,
                  max.only=T, ##this is for generateInitpts()
                  usezoom=F
) { ## Nb in metric units for metric plot.
    ## Likelihood profile POINT for fixed values of variables fixed, taking Out all variables not specified in fixed
  ## return value is list(CANON, purefn(canon), ......................)
  ## Values may be specified for fixed vars and for other variables.
  ## For the fixed variables, these values are mandatory (they specify for which value the profile is computed).
  ## For the other variables, these values are optional (they only propose starting values for the optimization algorithms),
  ## starting values WERE taken from maximizing values (canonVPoutput) if no value is explictly given.
  ## NOW trying to use convex hull as much as possible (see call to subchullWrapper() )

  ## Returns NA's for the canonical vector if bas>haut...  which is lethal for
  ## all minimizations algos that somehow use the distance between points

  ## optimization: starting from rosglobal may miss the profile maximum in some cases;
  ## In that case, starting from the maximum for previous Nb helps producing smooth figures.

  fittedNames <- blackbox.getOption("fittedNames")
  ## first some checks of the arguments
  if ("Nb" %in% c(names(otherlist), names(fixedlist))) {
    message.redef(c("(!) From profile(): 'Nb' in function argument is highly suspect"))
    stop.redef()
  }
  notinKgspace <- names(fixedlist) %w/o% fittedNames
  checknames <- notinKgspace %w/o% c(blackbox.getOption("ParameterNames"), "latt2Ns2", "Nratio", "NactNfounderratio", "NfounderNancratio")
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
  ## For extra composite variables
  ## we can test within-hullness in convex hull of Nb-transformed data,
  ## which does not make a convex hull in kriging space
  ## but would have been the hull if we had kriged in 2Ds2-space
  ###### => computing the hull for extra composite space
  ## do this only once in the whole FONKgpointls analysis for a given variable
  if(is.null(locchull)) locchull <- providefullhull(names(fixedlist))[[1]] ## this works also without any composite variable
  ## since it comes from providefullhull it must have the affectedConstraints component
  ## control of test of within subhullness
  if (length(notinKgspace)>0) { ## then default is as given in next line else it is double
    if (  blackbox.getOption("redundant.mode")=="defaultPrecision") {
      precision <- "double"
    } else {precision <- blackbox.getOption("redundant.mode")}
  } else precision <- "double"
  ######
  ## Now we know the variables of locchull ('full' dimensional, e.g. 2Nmu, 2Ds2, g), and this must exactly be those of bestProfiledOut + fixedlist !
  ## we first look whether there is something to profile
  ## It is much more efficient to add constraint to existing constraints $Hrep than using the vertices
  if (missing(subHull_V_and_H)) subHull_V_and_H <- subHullWrapper(vertices=locchull$vertices, equality=fixedlist,
                                                                  include.Hrepr=TRUE,
                                                                  precision=precision ## FR->FR argument added 10/2015
                                                                 )
  subvertices <- subHull_V_and_H$vertices
  if (is.null(subvertices)) {  ## called when no point satisfy the constraints: safe exit from profile() with info
    # suggests that testPoint is out of kriged points...
    paramnbr <- blackbox.getOption("paramnbr")
    return(list(fixedlist=fixedlist,
                full=rep(NA, paramnbr),
                par=rep(NA, paramnbr),
                message="safe exit from profile() because is.null(subvertices)",
                value=NA, ## NA is important else -Inf -> numeric RelLik value -> plotted.
                inKgspace=FALSE,
                subHull_V_and_H=subHull_V_and_H)
    )
  } ## null subvertices exit
  ## else
  subverticesmax <- apply(subvertices, 2, max)
  subverticesmin <- apply(subvertices, 2, min)
  localmaxrange <- subverticesmax-subverticesmin
  if (nrow(subvertices)==1 ##point unique: 'profile' is easily computed
      || max(localmaxrange)<0.000001 ## several rows but very close to each over
      ## so we are at the edge of the hull and numerical problems could occur
      ## + we don't need more precision
      ## FR->FR replaced 0.00001 by 0.000001 03/2011... very ad hoc anyway
  ) { ##point unique: 'profile' is easily computed
    point <- colMeans(subvertices) ## don't lose colnames !
    vkrig <- tofullKrigingspace(point, fixedlist) ## conversion to full kriging param space
    canon <- canonize(vkrig)$canonVP ## completion/conversion to canonical
    if (length(notinKgspace)>0  || extrap==T) { ## the two cases where we don't already know the point to be in Kgspace
      inKgspace <- isPointInCHull(vkrig, constraints=blackbox.getOption("hulls")$Kgtotal[c("a", "b")])
    } else inKgspace <- TRUE
    return(list(message="", full=canon, par=point, value=purefn(vkrig, testhull=FALSE), subHull_V_and_H=subHull_V_and_H, inKgspace=inKgspace))
  } ## ELSE : non trivial profiling to perform
  #######################################################################
  ## currently (10/2011) precision controls
  ## (1) the call to subchullWrapper -> subchull -> redundant (but not the returnType which was always floating point)
  ## (2) whether constrOptim or constrOptimR is called...
  ## note that some code in subchull remains in double precision because there in no rational version of addHeq
  ## within hull if ui.x-ci>0 and Ax-b<0 i.e (-A)x+b>0 ie ui=(-A)= cols[-(1:2)] of the hrepr, ci= - col2
  ## subHullWrapper renvoit $Hrepr= qui inclut en colonnes b=-ci, et ui
  ui <- subHull_V_and_H$Hrepr[, -c(1, 2), drop=FALSE]
  ci <- - subHull_V_and_H$Hrepr[, 2] ##  within-hullness occurs for ui %*% x -ci => 0
  #### FR->FR a perhaps better way would be that subchullWrapper returns vertices as floating point, and constraints as <precision>,
  if (precision=="rational") {ui <- d2q(ui);ci <- d2q(ci)}
  ## note that subvertices was computed in double precision because there in no rational version of addHeq
  ## We need to find a starting point for maximization
  ##see the subvertices as a disk in 3D space
  ## if roglobal is close to the inside of a face
  ## a vector inside the disk deduced heuristically from rosglobal and otherlist can be a good start point
  ## we construct such a point, unsure to be in the convex hull
  ## next we consider the case were the likelihood can be maximized on the edge of the coin
  ## finally we explore the face of the coin
  profiledOutPars <- colnames(locchull$vertices) %w/o% names(fixedlist)
  ### first the unsure vector
  bestProfiledOut <- as.numeric(array(NA, length(profiledOutPars)))
  names(bestProfiledOut) <- profiledOutPars
  ### first the unsure vector
  locnames <- intersect(profiledOutPars, fittedNames)
  bestProfiledOut[locnames] <- templateinkg[locnames] ## rosglobal by default
  ## note that the scale of the hull should already be that of kriging, no additional log transf should be needed ##
  ## otherlist overrides the previous values
  if (!is.null(otherlist)) {
    locnames <- intersect(profiledOutPars, names(otherlist))
    bestProfiledOut[locnames] <- unlist(otherlist[locnames])
    ## note that in current use otherlist values come from a grid generated in Kriging space, so again no additional log transf should be needed ##
  }
  if (any(is.na(bestProfiledOut))) {
    message.redef("(!) any(is.na(bestProfiledOut))")
    message.redef(paste("fixedlist was", fixedlist))
    message.redef(paste("otherlist was", otherlist))
    message.redef(paste("colnames(locchull$vertices) were", colnames(locchull$vertices)))
  }
  ################
  ## designed to construct safe starting point(s) from an unsafe one
  candidates <- generateInitpts(bestProfiledOut=bestProfiledOut, vertices=subvertices, ui=ui, ci=ci, hrepr=subHull_V_and_H$Hrepr,
                              fixedlist=fixedlist,
                              precision=precision, max.only=max.only)
  #    candidates <- generateInitpts(bestProfiledOut=bestProfiledOut, locsubchull=subvertices, ui=ui, ci=ci, fixedlist=fixedlist,
  #                                                 precision=precision, max.only=max.only)
  ################
  if (usezoom) { ## then additionally constructs another starting point (found useful at least once for a local maximum)
    ## a slow but more general methodn that was apparently useful for a LRT on latt2Ns2 when for kriging on condS2 and Nm is better estimated than 2Ns2
    ## usezoom currently true only in LRTfn -> profile
    cP <- zoomProfile(fixedlist=fixedlist, extrap=extrap, locchull=locchull,
                    templateinkg=templateinkg, ## conversion to local hull within the function
                    precision=precision   ##FR->FR NO NEED FOR RATIONAL HERE ??
    ) ## returns in canonical space
    bestProfiledOut <- fromFONKtoanyspace(tofullKrigingspace(cP$par),
                                          colnames(locchull$vertices))[profiledOutPars]
    ################
    ## designed to construct safe starting point(s) from an unsafe one
    candidateszoom <- generateInitpts(bestProfiledOut=bestProfiledOut, vertices=subvertices, ui=ui, ci=ci, hrepr=subHull_V_and_H$Hrepr,
                                    fixedlist=fixedlist,
                                    precision=precision, max.only=max.only)
    #           candidateszoom <- generateInitpts(bestProfiledOut=bestProfiledOut, locsubchull=subvertices, ui=ui, ci=ci, fixedlist=fixedlist,
    #                                                 precision=precision, max.only=max.only)
    ################
    candidates <- c(candidates, candidateszoom)
  }
  ## end of construction of bestProfiledOut

  ### Finally we have (a set of) candidate(s) in the hull of the profile. We maximize
  if ("NLOPT_LN_COBYLA" %in% blackbox.getOption("optimizers")) {
    ## minimization of - logL
    objfn_nloptr <- function(x,ui,ci) { ## all functions should have the same args. ui and ci will be ignored in this fn
      names(x) <- names(bestProfiledOut) ## nloptr tends to lose names
      return(- tofKpredict.nohull(x, fixedlist=fixedlist))
    }
    eval_g_ineq <- function(x,ui,ci) {max(ci - ui %*% x)} ## must be <0
    try_nloptr_wrap <- function(...) {
      resu <- try(nloptr(...))
      if (! inherits(resu,"try-error")) {
        resu$par <-resu$solution
        names(resu$par) <- names(bestProfiledOut)
        resu$value <- - resu$objective # back to + logL
      }
      return(resu)
    }
  } else {
    ## maximization of objfn = + logL     (control$fnscale=-1)
    objfn <- function(x) {tofKpredict.nohull(x, fixedlist=fixedlist)} ## the simple way to pass fixedlist to the inner computations...
    objfn.grad <- function(x) {grad(func=objfn, x=x)} ## no need to specify fixedlist here
  }
  ##  initial value in 'vmin' (or 'vmmin') is not finite: function fun:=R(theta, theta.old, ...) may have returned Inf which means it considers the point is not in the hull represented by ui, ci
  control <- list(fnscale = -1/blackbox.getOption("scalefactor"), trace = FALSE, maxit = 10000)
  parscale <- localmaxrange
  control <- c(control, list(parscale=parscale/blackbox.getOption("scalefactor"))) ## don't forget scalefactor here as it's in fnscale...
  bestresu <- list(value=- Inf)
  for (candidate in candidates) {
    bestProfiledOut <- candidate$par
    if (is.matrix(bestProfiledOut)) stop.redef("(!) is.matrix(bestProfiledOut)")
    if (is.null(names(bestProfiledOut))) stop.redef("(!) is.null(names(bestProfiledOut))")
    if (precision=="rational") {
      resu <- try(constrOptimR(unlist(bestProfiledOut), objfn, grad=objfn.grad, ui=ui, ci=ci , mu=1e-08, ## a low mu appear important
                             method = "BFGS",
                             control = control)
      )
    } else {
      if ("NLOPT_LN_COBYLA" %in% blackbox.getOption("optimizers")) {  ## not efficient
        resu <- try_nloptr_wrap(x0=unlist(bestProfiledOut),eval_f=objfn_nloptr,eval_g_ineq=eval_g_ineq,ui=ui,ci=ci,
                            opts=list(algorithm="NLOPT_LN_COBYLA",xtol_rel=1.0e-5,maxeval=-1)) ## BOBYQA = only BOund constraints
        ## dans constrOptim, le gradient est essentiel pour que l'optimisation contrainte marche
        ## Un constrOptim qui appelle NLOPT devrait donc utiliser une method avec gradient
      } else {
        resu <- try(constrOptim(unlist(bestProfiledOut), objfn, grad=objfn.grad, ui=ui, ci=ci , mu=1e-08, ## a low mu appear important
                                #                                       localmaxrange=localmaxrange, ## arg de tofKpredict.nohull.grad
                                method = "BFGS",
                                #                                       fixedlist = fixedlist,
                                control = control)
        )
      }
    }
    if (! inherits(resu,"try-error") && resu$value>bestresu$value) bestresu <- resu
  } ## end loop over candidate starting points
  if (bestresu$value==-Inf) {
    ## parameter value a edge or outside kriged range ? Seems to work most of the time, but far from always, when at edge
    FONKgNames <- blackbox.getOption("FONKgNames")
    zut <- list(full=rep(NA, length(FONKgNames)), par=rep(NA, length(FONKgNames)), value=NA, message="no valid constrOptim(R) return value within profile().",
              inKgspace=F,
              edgelevel=0, edge=list(), subHull_V_and_H=subHull_V_and_H)
    return(zut)
  }
  # ELSE
  resu <- bestresu
  ## One problem, however, is that such algos are a bit stuck when they meet a face of the hull
  ## so we will find whether this is so, and if so we will optimize on the face found
  ## (it is still possible that we are on a wrong face... => rr=0 or 1...)
  fullerhull <- matchVertCons(locchull) ## In grid calls, time will be saved if locchull already contains the affectedConstraints.
  edgecheck <- handlesPointAtEdge(point=resu$par, fullerhull=fullerhull, fixedlist=fixedlist)
  if (!is.null(edgecheck$edge)) {
    locedge <- t(apply(edgecheck$edge, 1, tofullKrigingspace, fixedlist=fixedlist))
    colnames(locedge) <- blackbox.getOption("FONKgNames")
    addtoedges(, locedge)
  }
  if (!is.null(edgecheck$resu)) {
    if (edgecheck$resu$value>resu$value) {
      resu <- edgecheck$resu
    } ## else resu near the edge remains better that the best point on the edge => nothing to do
  }
  ## resu$par is a vector of parameters profiled out
  vkrig <- tofullKrigingspace(resu$par, fixedlist) ## conversion to full kriging param space
  canon <- canonize(vkrig)$canonVP ## completion/conversion to canonical
  zut <- resu;
  if (length(notinKgspace)>0  || extrap==T) {  ## the two cases where we don't already know the point to be in Kgspace
    inKgspace <- isPointInCHull(vkrig, constraints=blackbox.getOption("hulls")$Kgtotal[c("a", "b")])
  } else inKgspace <- T
  ## 'full' must be a suitable argument for tofullKrigingspace; zut$par is the $par element of the return value of constrOptim
  zut <- c(zut, list(full=canon, message="", edgelevel=edgecheck$edgelevel, inKgspace=inKgspace, edge=list(edgecheck$edge), subHull_V_and_H=subHull_V_and_H))
  return(zut) ## zut$canon is vector in canonical param space while zut$par (not always in return value) is vector of parameters profiled out
} ## end profile(...)
