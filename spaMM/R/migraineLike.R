#require(geometry) ## delaunayn,convhulln


upperPoints <- function(Xpredy,minPtNbr=1#sqrt(nrow(Xpredy)+1)
                        ,D.resp) {
  respName <- attr(Xpredy,"fittedName")
  orderResp <- order(Xpredy[,respName], decreasing = TRUE)
  orderedXy <- Xpredy[orderResp,]
  Dy <- orderedXy[1,respName]-orderedXy[,respName]
  ptNbr <- max(c(which(Dy<D.resp)),minPtNbr) 
  orderResp[seq_len(ptNbr)] ## now returns indices
}

#nH <- upperPoints(predictions)

elim.redundant.V <- function(vertices) { ## removes redundant vertices
  if (nrow(vertices)<=ncol(vertices)) { ## convhulln crashes!
    minimalvertices<-vertices ## stupid case, should never occur
  } else if (ncol(vertices)==1L) {
    minimalvertices<-array(c(min(vertices),max(vertices)),dim=c(2,1))
  } else {
    minimalvertices<-vertices[unique(as.numeric(convhulln(vertices, "Pp"))),] ## removes redundant vertices
  }
  colnames(minimalvertices)<-colnames(vertices)
  minimalvertices
}

#mvH <- elim.redundant.V(nH)

volTriangulation <- function(vertices) { ## 
  if (is.data.frame(vertices)) vertices <- as.matrix(vertices)
  ## probleme with repeated vertices occurs sometimes:
  vertices <- unique(vertices) 
  ## ...otherwise it is possible that a vertex (from $vertices) is later selected, which is not in the $simplicesTable
  tc <- delaunayn(vertices,"Pp") ## triangulation by simplices
  # cf blckbox::volTriangulationWrapper for cases where delaunayn fails
  pmul <- cbind(-1,diag(ncol(vertices)))
  factorialdim <- factorial(ncol(vertices)) ## which is not always try if only a subspace is sampled
  vb <- apply(tc,1,function(v){
    simplex <- vertices[v,,drop=FALSE]
    # pmul %*% simplex is simplex[-1,]-simplex[1,] for each simplex
    # volume = abs(det())/ dim!
    vol <- abs(det(pmul %*% simplex))/factorialdim
    bary <- colMeans(simplex)
    c(vol,bary) ## volume, and barycenter \n\
  }) 
  resu <- list(vol=vb[1,], ## a vector
               bary=t(vb[-1,,drop=FALSE]), ## matrix
               vertices=vertices,simplicesTable=tc)
  class(resu) <- c("volTriangulation",class(resu))
  resu
}
#vT <- volTriangulation(mvH)


## 'simplices' are indices kept, or subtracted if negative
# originally named subset.volTriangulation but subset is a generic with a different usage
subsimplices.volTriangulation <- function(x,simplices,...) {
  vT <- x
  vT$vol <- vT$vol[simplices]
  vT$bary <- vT$bary[simplices,]
  vT$simplicesTable <- vT$simplicesTable[simplices,,drop=FALSE]
  # vT$simplicesTable refers to original point indices => $vertices is unchanged
  vT
}


## to sample uniformly wrt volume within a simplex, not uniformly wrt perimeter
# old complicated code with additional argument 'u' removed 04/2016
#
rsimplex <- function(simplex,
                     expand=NULL, ## *volume* expansion 
                     bary=NULL ## used for expansion
) {   
  d <- NROW(simplex)-1 ## 2 pour triangle
  if(NCOL(simplex)!= d) {stop("(!) From 'rsimplex': simplex should have one more row than columns.")}
  weights <- diff(sort(c(0, runif(n=d), 1))) ## Wilks 1962, p.238 in Rubin 1981
  # http://cs.stackexchange.com/questions/3227/uniform-sampling-from-a-simplex for some more discussion
  if ( ! is.null(expand)) {
    if (is.null(bary)) bary <- colMeans(simplex) ## OK pour bary de d-dim simplex avec densité uniforme 
    simplex <- bary + expand^(1/d) * sweep(simplex,2L,bary,`-`)
  }
  ws <- sweep(simplex,1L,weights,`*`)
  return(colSums(ws))
}

## ideally we would have a rvolume S3 generic with .volTriangulation and .data.frame methods
rvolTriangulation <- function(n=1,volTriangulationObj,replace=TRUE,expand=NULL) { 
  if (! inherits(volTriangulationObj,"volTriangulation")) stop("(!) From 'volTriangulation': input 'volTriangulationObj' is not a volTriangulation object.")
  simplexProbs <- volTriangulationObj$vol
  simplexProbs <- simplexProbs/sum(simplexProbs)
  whichSimplex <- sample(seq_len(length(simplexProbs)),n,prob=simplexProbs,replace=replace)
  vertices <- volTriangulationObj$vertices
  vI <- volTriangulationObj$simplicesTable
  resu <- sapply(whichSimplex,function(idx) {rsimplex(vertices[vI[idx,],,drop=FALSE],expand=expand)})
  if (ncol(vertices)==1L) {
    resu <- as.matrix(resu)
  } else resu <- t(resu)
  rownames(resu) <- NULL ## bc all rownames=parm when ncol=1
  colnames(resu) <- colnames(vertices,do.NULL=FALSE) 
  return(resu)
} ## end def rvolTriangulation


# old.nextPoints <- function(n=1,optr,replace=TRUE) { ## random sampling of volume defined from previous fit
#   uP <- upperPoints(optr$predictions) ## indices
#   uP <- optr$predictions[uP,attr(optr$predictions,"fittedPars")]
#   uP <- rbind(uP,optr$par) ## not sure this is useful for volumetric sampling
#   erV <- elim.redundant.V(uP)  
#   vT <- volTriangulation(erV)
#   rvolTriangulation(n,vT,replace=replace)
# }




spaMM_rhullByEI <- function(n, tryn=100*n ,vT, object, fixed=NULL, outputVars,expand=NULL) {
  ## so that oldXnew matrices will contain less than spaMM.getOption("ff_threshold")~1e7 elements,
  maxn <- floor(spaMM.getOption("ff_threshold")/nrow(object$data))
  if (maxn <= n) {
    locmess <- paste("From spaMM_rhullByEI(): 'maxn': ",maxn,"<=",n," ('n'). 'n' reduced to")
    n <- ceiling(maxn/10)
    message(paste(locmess,n))
  }
  if (tryn > maxn) {
    locmess <- paste("From spaMM_rhullByEI(): 'tryn' reduced from",tryn,"to",maxn)
    message(locmess)
    tryn <- maxn
  }
  trypoints <- data.frame(rvolTriangulation(tryn, vT,expand=expand))
  colnames(trypoints) <- colnames(vT$vertices) ## supposeque non null...
  if (! is.null(fixed)) trypoints <- cbind(trypoints, fixed)
  colnames(trypoints) <- outputVars ## 'apply' feature
  trypred <- predict(object, newdata=as.data.frame(trypoints), variances = list(linPred=TRUE))
  trySE <- attr(trypred, "predVar")
  trySE[trySE<0] <- 0
  trySE <- sqrt(trySE)
  tryQ <- trypred + 1.96*trySE ## improvement function for candidate points
  Qmax <- object$Qmax
  expectedImprovement <- trySE*dnorm((Qmax-tryQ)/trySE)+(tryQ-Qmax)*pnorm((tryQ-Qmax)/trySE) ## 7.5 p. 121
  ## expected improvement (over current maximum) is never negative by def.
  trypoints <- trypoints[order(expectedImprovement, decreasing=TRUE)[seq_len(n)], outputVars, drop=FALSE]
  return(trypoints) 
}

## FR->FR same as in blackbox
locatePointinvT <- function(point, ## numeric (not matrix or data frame: see use of 'point' below)
                            vT,fallback=TRUE) { ## in which simplex ? with fall back if 'numerically outside' vT (but quite distinct from minimal distance)
  pmul <- cbind(-1,diag(ncol(vT$vertices)))
  minw <- apply(vT$simplicesTable[vT$vol>0,,drop=FALSE],1,function(v) {
    simplex <- vT$vertices[v,,drop=FALSE]
    vM <- pmul %*% simplex # is simplex[-1,]-simplex[1,] for each simplex
    vWeights <- try(solve(t(vM),point-simplex[1,]),silent=TRUE) ## as.numeric(point) would be required if point were a (1-row) matrix
    # problem may occur if volume of simplex is nearly zero
    if (inherits(vWeights,"try-error")) {vWeights <- ginv(t(vM)) %*% (point-simplex[1,])}
    vWeights <- c(1-sum(vWeights),vWeights) ## weights for all vertices
    min(vWeights) ## if the point is within/marginally outside/clearly outside the vT, this will return a positive/small negative/large negative value
  })
  resu <- which(minw>0)
  if (length(resu)==0L && fallback) resu <- which.max(minw) ## if numerically outside...
  return(resu)
}

spaMM_bounds1D <- function(object, optr, CIvar, precision) {
  fittedNames <- names(optr$par)
  np <- length(fittedNames)
  profiledNames <- setdiff(fittedNames, CIvar) 
  MLval <- optr$par[CIvar]
  givenmax <- optr$value
  lower <- apply(object$data,2,min)
  upper <- apply(object$data,2,max)
  lowval <- lower[CIvar]
  lowval <- lowval+0.002*(MLval-lowval)
  hival <- upper[CIvar]
  hival <- hival-0.002*(hival-MLval)
  ## def objectivefn
  shift <- givenmax - precision
  if(np==1L) {
    objectivefn1D <- function(x) { ## defined to return 0 at the CI threshold
      return(predict(object,x)+shift)
    }
    objectivefn <- objectivefn1D
  } else {
    argfixedlist <- eval(parse(text=paste("list(", CIvar, "=", NA, ")"))) ## je pourrais creer une list et nommer ensuite... mais cette syntaxe pourra resservir
    parv <- optr$par
    objectivefnmultiD <- function(x, return.optim=FALSE) { ## defined to return 0 at the CI threshold
      parv[CIvar] <- x
      objfn <- function(z) { 
        parv[profiledNames] <- z
        return(predict(object,parv))
      }
      locoptr <- optim(optr$par[profiledNames],fn=objfn,lower=lower[profiledNames],upper=upper[profiledNames],method="L-BFGS-B",
                       control=list(fnscale=-1,parscale=(upper-lower)[profiledNames]))
      if(return.optim) {
        return(locoptr)
      } else return(locoptr$value-shift) ## then returns shifted value for uniroot
    }
    objectivefn <- objectivefnmultiD
  }
  CIlo <- NA
  CIpoints <- matrix(nrow=0, ncol=np)
  fupper <- precision
  if (abs(lowval-MLval)<abs(MLval*1e-08)) { ## not an error ## lowval==MLval if relative diff is <1e-16 in absolute value
    # CIlo <- NA ## already so
  } else {
    flower <- objectivefn(lowval)
    if (is.na(flower)) { ## abnormal case
      stop("From spaMM_bounds1D: 'flower' is NA")
    } else {
      if (flower<0) {CIlo <- try((uniroot(objectivefn, interval=c(lowval, MLval),
                                          f.lower=flower, f.upper=fupper))$root, TRUE)
      }
    }
    if(inherits(CIlo,"try-error")) {
      CIlo <- NA
      errmsg <- paste("Bound", " for ", CIvar, " could not be computed (maybe out of sampled range)", sep="")
      message(errmsg)
    }
    if ( ! is.na(CIlo)) {
      CIpoint <- optr$par
      CIpoint[CIvar] <- CIlo
      if (length(profiledNames)>0) {
        thisfit <- objectivefn(CIlo, return.optim=TRUE)
        CIpoint[profiledNames] <- thisfit$par
      }
      CIpoints <- rbind(CIpoints, CIpoint)
    }
  }
  CIup <- NA ## default values...
  flower <- fupper
  if (abs(hival-MLval)<abs(MLval*1e-08)) { ## not an error
    ## CIup <- NA
  } else {
    fupper <- objectivefn(hival)
    if (is.na(fupper)) { ## abnormal case
      stop("From spaMM_bounds1D: 'fupper' is NA")
    } else {
      if (fupper<0) {CIup <- try((uniroot(objectivefn, c(MLval, hival),
                                          f.lower=flower, f.upper=fupper))$root, TRUE)
      }
    }
    if(inherits(CIup,"try-error")) {
      CIup <- NA
      errmsg <- paste("Bound", " for ", CIvar, " could not be computed (maybe out of sampled range)", sep="")
      message(errmsg)
    }
    if ( ! is.na(CIup)) {
      CIpoint <- optr$par
      CIpoint[CIvar] <- CIup
      if (length(profiledNames)>0) {
        thisfit <- objectivefn(CIup, return.optim=TRUE)
        CIpoint[profiledNames] <- thisfit$par
      }
      CIpoints <- rbind(CIpoints, CIpoint)
    }
  }
  return(CIpoints=CIpoints)
} ## end def bounds1D


# Returns sampled locations
# default minPtNbr affect exploration of (provisionally) suboptimal peaks
# networkExpand=2 allows expansion at least towards a local peak
sampleNextPars <- function(sizes,optr,minPtNbr=1,replace=TRUE,networkExpand=1L,simplexExpand=NULL,D.resp) {
  ## derived from blackbox::sampleByResp, 2015/11
  ###################### slow precomputation for EI
  Krigobj <- optr$Krigobj
  obspred <- predict(Krigobj, variances=list(linPred=TRUE))
  obsSE <- attr(obspred, "predVar")
  obsSE[obsSE<0] <- 0
  Krigobj$Qmax <- max( obspred+1.96 * sqrt(obsSE)) ## best improvement function for already computed points
  ######################
  Xpredy <- optr$predictions
  fittedPars <- attr(Xpredy,"fittedPars")
  vT <- volTriangulation(as.matrix(Xpredy[,fittedPars,drop=FALSE])) 
  if (is.infinite(networkExpand)) {
    upperSimplices <- seq_len(nrow(vT$simplicesTable))
    innerVertexIndices <- seq_len(nrow(Xpredy)) 
  } else {
    ## seek simplices which involve 'interesting' points
    innerVertexIndices <- upperPoints(Xpredy,minPtNbr=minPtNbr,D.resp=D.resp) ## indices; FR->FR maybe not most efficient as Xpredy is already sorted
    ## the following defines a NON CONVEX set
    upperSimplices <- apply(vT$simplicesTable,1,function(v) { any(v %in% innerVertexIndices)}) ## indices !
    for (it in seq_len(networkExpand-1)) {
      innerVertexIndices <- unique(as.vector(vT$simplicesTable[upperSimplices,]))
      upperSimplices <- apply(vT$simplicesTable,1,function(v) { any(v %in% innerVertexIndices)}) ## indices !    
    }
  }
  subvT <- subsimplices.volTriangulation(vT,simplices=upperSimplices)
  goodpoints <- data.frame(rvolTriangulation(n=sizes[1L],subvT,replace=replace,expand=simplexExpand))  ## because (rbind(matrix,..., data.frame)) can generate empty rownames...
  ## some exploratory code; precision is here a shift on response (lik); fails if nu is fixed.
  if (FALSE) {
    nubounds1D <- spaMM_bounds1D(Krigobj, optr, "trNu", precision=0.1) ## FR->FR arbitrary threshold
    goodpoints <- rbind(goodpoints,nubounds1D) 
  }
  ######################
  subvT <- subsimplices.volTriangulation(vT, - upperSimplices)
  candidates <- spaMM_rhullByEI(n=sizes[2L],vT=subvT,object=Krigobj,outputVars = fittedPars)
  goodpoints <- rbind(goodpoints,candidates )
  ###################### => optr$par (6) + EI (others)
  resu <- data.frame(goodpoints)
  #attr(resu,"info") <- vT
  return(resu)
}
# attribute "info" would contain 
#    the original $vertices (indices), 
#    $upperVertexIndices the indices of the best points, 
#    $simplicesTable the subset of simplices involving the best points
#    $vol the volumes of subset of simplices involving the best points


connectedSets <- function(indices) {
  uI <- unique(as.vector(indices)) ## vertex indices
  rownames(indices) <- seq_len(nrow(indices))
  simplicesTable <- list()
  it <- 1
  while (nrow(indices)>1) {
    oldvertexSet <- numeric(0)
    growingvertexSet <- uI[1] # a vertex
    while(length(oldvertexSet)!=length(growingvertexSet)) {
      oldvertexSet <- growingvertexSet
      rowIndices <- which(apply(indices,1,function(v) { any(v %in% growingvertexSet)}))
      connectedRows <- indices[rowIndices,]
      growingvertexSet <- unique(as.vector(connectedRows))
    }
    simplicesTable[[it]] <- growingvertexSet
    uI <- setdiff(uI,growingvertexSet)
    indices <- indices[-rowIndices,,drop=FALSE]
    it <- it+1
  }
  if (nrow(indices)==1L) simplicesTable[[it]] <- as.numeric(indices) ## isolated simplex
  lapply(simplicesTable,as.numeric) 
} ## returns point indices

extrapolhull <-function(Vhull,extrapol=1.4) { ##
  if (nrow(Vhull) <2) return(matrix(nrow=0,ncol=ncol(Vhull)))
  ## ELSE
  vertexbary <- colMeans(Vhull) ## not bary for uniform density, except for simplices
  ## extremely lightweight solution: random extrapol along all directions defined by the vertices and the bary
  extrap <- runif(nrow(Vhull))*(extrapol-1)+1
  t((t(Vhull)-vertexbary)* extrap+ vertexbary)  ## bary+ extrapol*(v-bary)
}  