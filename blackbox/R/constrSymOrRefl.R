## edge eq y=edge[1, ]+sum_j b_j t(orthog_j)
constrSymOrRefl <- function(point, refpoint=c(), orthog=c(), logLcheck, lower, upper, edge) { ## constrained symmetry or reflexion
  if ( length(refpoint)>0) { ## this seems the obvious thing to do, was surely tried at some stage, then abandoned -- come back 08/12/10
    proj <- refpoint
  } else { ## orthog must have nonzero length! ## old option still directly using the edge but refpoint should be used
    ## proj is projection of (point-edge[1, ]) on the orthogonal vectors -> orthogonal projection of point on plane
    proj <- edge[1, ]+((point-edge[1, ]) %*% orthog) %*% t(orthog)
  }
  refl <- 2*proj-point # unconstr. reflexion
  resu <- refl
  ## previously test logLcheck!=-Inf
  if ( TRUE ) { ## if refl has too low predicted likelihood, we seek a closer point with predicted likelihood logLcheck
    ## first seeks the extreme allowable point within lower, upper
    ## 'solves' point+z*(refl-point)=allowable...
    z1 <- (upper-point)/(refl-point) ## for each element: num always positive, upper bound matters only if element positive
    z2 <- (point-lower)/(point-refl) ## idem
    z <- c(z1, z2)
    ## Inf's are OK here
    ## NaN can occur if point at edge (refl-point contains 0)
    z <- min(1, z[z>=0], na.rm=T)
    resu <- point+z*(refl-point)
    ## following is very crude algorithm but it assumes nothing about the shape of the likelihood
    ## It will likely stop.redef() if the logL (input 'point') is above the logLcheck. Hence it should first be checked that this is not the case before the call to mirror()
    pred <- tofKpredict.nohull(resu, fixedlist=NULL) ## bc resu is in sampling space, not in kriging space
    bon <- ( (!is.na(pred)) & (pred>  logLcheck))
    it <- 0
    while (!bon) {
      z <- z/2
      resu <- point+z*(refl-point) ## doesnot get the names(point)...
      names(resu) <- names(point)
      pred <- tofKpredict.nohull(resu, fixedlist=NULL)
      bon <- ( (!is.na(pred)) & (pred>logLcheck ))
      it <- it+1
      if(it>100) {
        message.redef("(!) Suspicious loop in reflexion(...)")
        shouldbehighenough <- tofKpredict.nohull(resu, fixedlist=NULL)
        ## but boundcheck was undefined in migraine 0.4...
        #if (shouldbehighenough<=boundcheck) {
        #  message.redef("   A wrong input 'point' with too low likelihood was apparently provided")
        #  message.redef(paste("   point: ", point))
        #  message.redef(paste("   predicted likelihood: ", shouldbehighenough))
        #  message.redef(paste("   likelihood threshold: ", boundcheck))
        #} else {
        #  message.redef("    Yet the input 'point' is apparently correct.")
        #}
        stop.redef()
      }
    } ## while !bon
  } ##!=-Inf
  return(resu)
} ## end def reflection within mirror
