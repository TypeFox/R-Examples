sampleStep <- function(size, step,
                       locLRTlist=blackbox.getOption("LRTlist"),
                       FONKgNames=blackbox.getOption("FONKgNames"),
                       threshold,
                       lowextrapol=1.1,
                       knots_canon=NULL,
                       knotsInfo=NULL,
                       fitforEI=NULL,
                       previous=NULL,
                       extrapol=if(step=="expand") {NULL} else {1}, ## (NULL will be replaced)
                       verbose=FALSE
                       ) {
  INFO <-blackbox.options()[c("fitobject","rosglobal")]
  resu <- list()
  #
  if (step=="expand") {
    resu$extrapol_dlr <- set_extrapol_dlr_from_LRTs(locLRTlist=locLRTlist,dlr=threshold,
                                               FONKgNames=FONKgNames,lowextrapol=lowextrapol,verbose=verbose)
    threshold <- resu$extrapol_dlr$dlr
    extrapol <- resu$extrapol_dlr$extrapol
  }
  #
  if (verbose) {
    message("<<<")
    message.redef(paste("step=",step,"; extrapol=",extrapol))
    message.redef(paste("Expected improvement will",c(" "," NOT ")[is.null(fitforEI)+1],"be used",sep=""))
  }
  if  (is.null(knotsInfo)) resu$knotsInfo <- knotsInfo <- calcKnotsInfoWrapper(fitobject=INFO$fitobject,
                                            rosglobal=INFO$rosglobal,  ## INFO$rosglobal FR->FR dangereux
                                            dlr=threshold,
                                            locLRTlist=locLRTlist,
                                            FONKgNames=FONKgNames,
                                            knots_canon=knots_canon,
                                            verbose=verbose)
  #
  resu$cumulgoodpoints <- generateNewPoints(size,
                                      extrapol=extrapol,
                                      step=step,
                                      knotsInfo=knotsInfo,
                                      fitforEI=fitforEI,
                                      previous=previous,
                                      verbose=verbose)
  if (verbose) message(">>>")
  return(resu)
}

sampleByResp <- function(size=blackbox.getOption("nextPointNumber"),
                         outfile = NULL,
                         useEI,
                         NextBoundsLevel=0.001,
                         threshold=qchisq(1-NextBoundsLevel, 1)/2,
                         rnd.seed=NULL,
                         verbose=FALSE ## verbosity for development purposes
) {
  if ( ! is.null(outfile)) unlink(outfile) ## this *deletes* *now* the file if it exists, before possible failure of this function.
    
  INFO <- blackbox.options()[c("ycolname","fittedNames","pureRMSE","lambdaEst","FONKgpointls","fitobject",
                               "CovFnParam","CIlevel","CIpointsList","rosglobal","FONKgLow","FONKgNames","ParameterNames")]

  ycolname <- blackbox.getOption("ycolname")
  message.redef("\n*** Generating points for next iteration ***")
  if ( ! is.null(rnd.seed)) set.seed(rnd.seed)
  ###################### slow precomputation for EI => silent code
  form <- as.formula(paste("`", INFO$ycolname, "`~1 + Matern(1|", paste(INFO$fittedNames, collapse="+"), ")", sep=""))
  ## lambda GCV = phi hglm/lambda_HGLM
  locphi <- max(1e-06,INFO$pureRMSE^2)
  loclambda <- locphi/INFO$lambdaEst
  if (verbose) { message.redef("Computing predictor for EI...")  }
  spaMMfit <- corrHLfit(form, data=INFO$FONKgpointls,
                        ranFix=list(rho=1/INFO$CovFnParam[INFO$fittedNames],
                                    #         note '1/...'
                                    nu=INFO$CovFnParam["smoothness"],
                                    lambda=loclambda,phi=locphi))
  obspred <- predict(spaMMfit, variances=list(linPred=TRUE))
  obsSE <- attr(obspred, "predVar")
  obsSE[obsSE<0] <- 0
  ## attention le signe de la prediction est inversé... SPECIFIQUEMENT POUR MIGRAINE
  spaMMfit$Qmax <- max( - obspred+1.96 * sqrt(obsSE)) ## best improvement function for already computed points
  #
  ######################
  ###################### First sampling using expansion
  cumul <- sampleStep(floor(size/5),step="expand",threshold=threshold,verbose=verbose)
  
  goodpoints <- cumul$cumulgoodpoints # (2015/09/21): these goodpoints define an expanded hull of points that extrapolate
  plotcolors <- rep("black",NROW(goodpoints)) ## black= expanded hull
  ## precomputation from global knots(not expanded, but uses lr threshold ~dlr*([dlrtolerance=1.2] + [probErr of roglobal$value])
  globalknotsInfo <- cumul$knotsInfo # knots Info that includes rosglobal$par
  redundvT <- volTriangulation(globalknotsInfo$redundantknots) ## contains rosglobal; would $knotsVH$vertices would provide less stratified sampling?
  redundpred <- predict(INFO$fitobject,
                        x=redundvT$vertices, ## that is, unique(globalknotsInfo$redundantknots)
                        testHull=FALSE)
  ## should consider a narrow upper range here:
  topvertices <- which( max(redundpred)-redundpred < qchisq(1-INFO$CIlevel/2,df=1)/2 ) ## certain to contain rosglobal...
  whichSimplices <- apply(redundvT$simplicesTable,1,function(v) any(v %in% topvertices))
  ###################### => expand (size/5)
  ###################### points from CIs
  CIpointsList <- INFO$CIpointsList
  nCIpts <- sum(unlist(lapply(CIpointsList,nrow)))
  if( nCIpts>0L ) {
    nPerCIpoint <- min(floor(size/(2*nCIpts)),6) ## /2 to let at least half of the points be sampled differently
    fromCI <- generateNextpointsfromCI(n=nPerCIpoint,
                                       CIpointsList = CIpointsList,
                                       posdlr= threshold,
                                       previous=goodpoints,
                                       fitforEI=spaMMfit,
                                       LowUp=globalknotsInfo$LowUp,
                                       verbose=verbose
    )
    goodpoints <- fromCI$cumulgoodpoints
    plotcolors <- c(plotcolors,rep("red",NROW(goodpoints)-length(plotcolors))) ## red= CI
  }
  ###################### => expand (size/5) +CI (<size/2)
  ###################### As it says: Sampling near inferred maximum
  if (verbose) {message.redef("Sampling near inferred maximum...")}
  whichSimplex <- locatePointinvT(INFO$rosglobal$par,redundvT) ## INFO$rosglobal FR->FR dangereux
  n6 <- min(6,(size-nrow(goodpoints))) ## 6 in any good usage
  if (length(whichSimplex)==1L) {
    candidates <- t(replicate(n6,rsimplex(simplex=redundvT$vertices[redundvT$simplicesTable[whichSimplex,],,drop=FALSE])))
  } else {
    subvT <- subsimplices.volTriangulation(redundvT,whichSimplex)
    candidates <- rvolTriangulation(n=n6,subvT)
  }
  candidates <- toCanonical(candidates,FONKgLow = INFO$FONKgLow,
                            othernames=INFO$FONKgNames %w/o% INFO$fittedNames)
  goodpoints <- rbind(goodpoints,candidates )
  plotcolors <- c(plotcolors,rep("orange",NROW(candidates))) ## orange : near inferred maximum
  ###################### => expand (size/5) +CI (<size/2) + rosglobal (6)
  ######################
  if (verbose) {message.redef("Non-convex sampling...")}
  if (any(! whichSimplices)) { ## s'il y a des simplex en dehors du sommet
    subvT <- subsimplices.volTriangulation(redundvT, - whichSimplices)
    candidates <- rhullByEI(n=ceiling(size/10),vT=subvT,object=spaMMfit)
    candidates <- toCanonical(candidates,FONKgLow = INFO$FONKgLow,
                              othernames=INFO$FONKgNames %w/o% INFO$fittedNames)
    goodpoints <- rbind(goodpoints,candidates )
    plotcolors <- c(plotcolors,rep("cyan",NROW(candidates))) ## cyan= EI outside the "top" region
  }
  ###################### => expand (size/5) +CI (<size/2) + rosglobal (6) + EI (size/10)
  ###################### filling "top" simplices:
  subvT <- subsimplices.volTriangulation(redundvT,whichSimplices)
  candidates <- rvolTriangulation(n=ceiling((size-nrow(goodpoints))*3/4),subvT)
  candidates <- toCanonical(candidates,FONKgLow = INFO$FONKgLow,
                            othernames=INFO$FONKgNames %w/o% INFO$fittedNames)
  goodpoints <- rbind(goodpoints,candidates )
  plotcolors <- c(plotcolors,rep("green",NROW(candidates))) ## green = non convex sampling of top (as def'd by CI threshold) vertices
  ###################### => expand (size/5) +CI (<size/2) + rosglobal (6) + EI (size/10) + top (3/4 or rest)
  ###################### Then use inner fill on 'global' hull
  cumul <- sampleStep(size,
                         step="fill",
                         threshold=NULL, ## not used in sampleStep bc knotsInfo available
                         knotsInfo=globalknotsInfo,
                         previous=goodpoints,
                         verbose=verbose)
  goodpoints <- cumul$cumulgoodpoints
  plotcolors <- c(plotcolors,rep("blue",NROW(goodpoints)-length(plotcolors))) ## blue= fill global hull (not expanded)
  ###################### => expand (size/5) +CI (<size/2) + rosglobal (6) + EI (size/10) + top (3/4 or rest) + fill (1/4 of rest)
  ## "top" is in region <= than "fill" is, but the two samples may be the same vertices, 
  # and there are more "top" than "fill" points, so that "fill" ponts appears within "top" points.
  ## diagnostic plot
  provideDevice(bbdefaultPars=TRUE)
  plot(goodpoints,col=plotcolors,pch=20,cex=1.25)
  technicolorTitle(c("expand", "CI", "max","EI","top","fill"), c("black", "red","orange","cyan","green","blue"),line=2.5)
  ordergdpts <- do.call(order, goodpoints)
  ## output to file
  goodpoints <- goodpoints[ ordergdpts , , drop=FALSE]
  if ( ! is.null(outfile)) {
    write(t(goodpoints), file=outfile, ncolumns=ncol(goodpoints)) ## t() still needed here...
    write(paste("# initCovFnParam=",paste(INFO$CovFnParam,collapse=" ")), file=outfile, append=TRUE)
  }
  invisible(goodpoints)
}

# Duncan Murdoch, http://r.789695.n4.nabble.com/title-words-in-different-colors-td878698.html
technicolorTitle <- function(words, colours, cex=1, side=3, line=1) {
  widths <- strwidth(words,cex=cex)
  spaces <- rep(strwidth(" ",cex=cex), length(widths)-1)
  middle <- mean(par("usr")[1:2])
  total <- sum(widths) + sum(spaces)
  start <- c(0,cumsum(widths[-length(widths)] + spaces))
  start <- start + middle - total/2
  mtext(words, side=side, line=line, at=start, adj=0, col=colours,cex=cex)
}

# utility for quick plot of parameter points
plotParPoints <- function(parPoints="nextpoints_1.txt",names=blackbox.getOption("ParameterNames"),minrange=1e-06) {
  if (is.character(parPoints)) {
    parPoints <- read.csv(parPoints,sep=" ",header=FALSE)
    colnames(parPoints) <- names
  }
  if (! inherits(parPoints, "data.frame")) parPoints <- data.frame(parPoints)
  ranges <- apply(parPoints,2,range)
  parPoints <- parPoints[,range[2,]-range[1,]>minrange]
  plot(parPoints)
  invisible(NULL)

}
