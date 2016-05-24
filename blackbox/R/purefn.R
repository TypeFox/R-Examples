purefn <- function(z, fitobj=blackbox.getOption("fitobject"), testhull=T, locus=NULL, constraints=blackbox.getOption("hulls")$Kgtotal[c("a", "b")], ...) {
  ## returns NA or not if out of convex hull depending on testhull=T/F
  ## for minimisation, do set testhull=F and use constrOptimM() !!
  if (is.data.frame(z) || is.matrix(z)) {
    if (nrow(z)>1) { ## useful for plot(function(z) {z <- matrix(z, ncol=1);colnames(z) <- ...;purefn(z)} as plot(fn) calls fn on a list of values
      yvals <- apply(z, 1, purefn, fitobj=fitobj, testhull=testhull, locus=locus, constraints=constraints, ...)
      return(yvals)
    } else { ##single point formatted as data frame or matrix
      colNames <- colnames(z)
    }
  } else if (is.list(z) || is.numeric(z)) {
    colNames <- names(z)
  } else {
    # do not return here NA or -Inf here because this generates an error in the sequence
    #  NA/NaN/Inf in foreign function call (arg 4)
    stop.redef(paste("(!) class of purefn 'z'= ", class(z), "; z: ", z))
  }
  if (any(colNames!= blackbox.getOption("fittedNames"))) {
    llocstring <- "(!) From purefn(): z argument not in kriging space (or col names missing)."
    message.redef(llocstring)
    bla <- paste(blackbox.getOption("fittedNames"))
    llocstring <- paste("   colnames(z): ", paste(colNames), "; krig Vars: ", bla)
    message.redef(llocstring)
    stop.redef()
  }
  if (testhull) {
    ## returning NA is OK for plots; minimisation algos require numeric values, not even 'Inf'
    if (!isPointInCHull(z, constraints=constraints)) return(NA)
  }
  if ("OKrig" %in% class(fitobj)) {
    return(predict(fitobj, rbind(z), locus=locus, ...))
  } else { ## should be a fitobject of class Kriglistplus
    dessous <- which(fitobj$blocmin<z[1]);dessus <- which(fitobj$blocmax>z[1]);blocs <- intersect(dessous, dessus)
    nblocsin <- length(blocs)
    if (nblocsin==1) {## point in only one block
      return(predict(fitobj$Kriglist[[blocs[1]]], rbind(z), bloc=blocs[1], locus=locus, ...))
    } else if (nblocsin==0) {## point in no bloc: use lowest or uppermost block
      if(length(dessous)>0) {
        nblocs <- length(fitobj$Kriglist)
        return(predict(fitobj$Kriglist[[nblocs]], rbind(z), bloc=nblocs, locus=locus, ...))
      } else {
        return(predict(fitobj$Kriglist[[1]], rbind(z), bloc=1, locus=locus, ...))
      }
    } else {## interpolates between two blocks
      lowpred <- predict(fitobj$Kriglist[[blocs[1]]], rbind(z), bloc=blocs[1], locus=locus, ...)
      hipred <- predict(fitobj$Kriglist[[blocs[2]]], rbind(z), bloc=blocs[2], locus=locus, ...)
      w1 <- (z[1]-fitobj$blocmax[blocs[1]]);w1 <- w1**2
      w2 <- (z[1]-fitobj$blocmin[blocs[2]]);w2 <- w2**2
      return((w1*lowpred+w2*hipred)/(w1+w2))
    }
  }
}
