calcGridRelProfile <- function(fixed, profileMethod="profileBySubHull",
                               gridsteps=blackbox.getOption("gridstepsNbr"), xGridloc=NA, yGridloc=NA) {
  ## applies profile() over a grid of values, recomputes ML estimates if needed,
  ## input is 'fixed', a vector of variable names [see profile(fixedlist, ...)]
  ## and gridsteps, the length of vectors of grid values for these two variables
  ## outputs *relative* profile likelihoods
  ## FAILS IF (fixed) are not exactly two (at most) of the boundsvars; with an exception for (Nb in fixed but not in boundsvars)

  INFO <- list(fittedNames=blackbox.getOption("fittedNames"),
               fittedparamnbr=blackbox.getOption("fittedparamnbr"),
               FONKgNames=blackbox.getOption("FONKgNames"),
               profile3passesBool=blackbox.getOption("profile3passesBool"))
  ## Determines x and y grids if not given in the arguments:
  if (is.na(xGridloc)) {xGrid <- gridfn(fixed[1], gridsteps, margefrac=1/(2*(gridsteps-1)), varnameS=fixed)} else {xGrid <- xGridloc}
  if (length(fixed)>1) {
    if (is.na(yGridloc)) {yGrid <- gridfn(fixed[2], gridsteps, margefrac=1/(2*(gridsteps-1)), varnameS=fixed)} else {yGrid <- yGridloc}
  } else {yGrid <- NA} ## and then, length(yGrid)=1... ## FR->FR should perhaps have used NULL here; when I have time...

  # grmf <- list(par=NA, value=NA, inKgspace=F)
  RelLik <- array(NA, c(length(xGrid), length(yGrid)))
  inKrigSpace <- array(TRUE, c(length(xGrid), length(yGrid))) ## an array of T/F, whether max is in hull
  if(profileMethod=="profileBySubHull") subhullinfo <- array(list(NULL), c(length(xGrid), length(yGrid))) ## an array of convex hull informations for each grid point
  testvars <- (fixed %w/o% INFO$fittedNames) ## should be empty... except in case adhocked below
  ## ad hoc fix for Nm bound /Nb profile
  if ("latt2Ns2" %in% fixed && "twoNm" %in% INFO$fittedNames) {testvars <- (testvars %w/o% "latt2Ns2")}
  ## ad hoc fix for Nb bound/Nm profile
  if ("twoNm" %in% fixed && "latt2Ns2" %in% INFO$fittedNames) {testvars <- (testvars %w/o% "twoNm")}
  ## ad hoc fix...
  if ("Nratio" %in% fixed && "twoNmu" %in% INFO$fittedNames) {testvars <- (testvars %w/o% "Nratio")}
  if ("Nancratio" %in% fixed && "twoNmu" %in% INFO$fittedNames) {testvars <- (testvars %w/o% "Nancratio")}
  if ("NactNfounderratio" %in% fixed && "twoNmu" %in% INFO$fittedNames) {testvars <- (testvars %w/o% "NactNfounderratio")}
  if ("NfounderNancratio" %in% fixed && "twoNancmu" %in% INFO$fittedNames) {testvars <- (testvars %w/o% "NfounderNancratio")}
  if(length(testvars)>0L) {
    locmess <- paste("'",fixed,"'",sep="",collapse=", ")
    locmess <- paste("Skipping profile grid computation for",locmess)
    message.redef(locmess)
    return()
  }
  ## Following code creates structure to store values of 'other' variables, to be used as
  ## starting point for optimization at neighboring values of given *grid (fixed)* variables
  othervars <- (INFO$fittedNames %w/o% fixed) ## kriging variables not in fixed
  ## ad hoc fix for Nm bound /Nb profile: only g will remain
  ##... this structure is not used in the first profiling but in later ones.
  ##... When Nb is in fixed, only the 2Nm value of othervals can be used in profiling()
  ##... but current behaviour for case where Nb is in boundvars too is
  ##... (1) not to include 2Nm in othervars, and
  ##... (2) to include g only in othervaras and profpts.
  ##... Given this, a bit more adhockery is needed.
  if ("latt2Ns2" %in% fixed && "twoNm" %in% othervars) {othervars <- (othervars %w/o% "twoNm")}
  ## FR->FR what if 2Nm in fixed and Nb in fittedNames ? does this occur and shouldn't this be handled in a parallel way ?
  if ("Nratio" %in% fixed && "twoNmu" %in% othervars) {othervars <- (othervars %w/o% "twoNmu")}
  if ("Nancratio" %in% fixed && "twoNmu" %in% othervars) {othervars <- (othervars %w/o% "twoNmu")}
  if ("NactNfounderratio" %in% fixed && "twoNmu" %in% othervars) {othervars <- (othervars %w/o% "twoNmu")}
  if ("NfounderNancratio" %in% fixed && "twoNancmu" %in% othervars) {othervars <- (othervars %w/o% "twoNancmu")}
  otherlist <- vector("list", length(othervars))
  names(otherlist) <- othervars
  poslogL <- length(othervars)+1L
  othervarspos <- seq_len(length(othervars))
  ## d'o[`u] : si xGrid est 2Nmu et yGrid est Nb, profpts[x value, y value] est c(g, logL) car Nb a pris la place de Nm
  profpts <- array(NA, c(length(xGrid), length(yGrid), poslogL)) ## profpts[i, j, .] is c(oth1, oth2, ..., logL)
  tempGridRelbestll <- blackbox.getOption("rosglobal")$value
  tempBestfull <- NULL
  ## computing profile log(Likelihood)
  ## first pass: estimate as starting point for maximisation
  if (length(fixed)>2L) {message.redef("(!) Length of 'fixed' is >2 in calcGridRelProfile()");return(NA)}
  fixedlist <- vector("list", length(fixed))
  names(fixedlist) <- fixed
  ntotal <- length(xGrid)*length(yGrid)
  ptit <- 1L
  prevmsglength <- 0L
  nProfileCalls <- 0L
  profilesTime <- 0
  t0 <- proc.time()["user.self"]
  usernames <- sapply(fixed, formatName, format="ASCII")
  prevmsglength <- 0L

  if ( ! INFO$profile3passesBool) {
    passes <- 1L
  } else if (length(fixed)==INFO$fittedparamnbr) {  ## no profile
    passes <- 1L
  } else if (length(fixed)==1L) {
    passes <- 1L ## but more passes could make sense for 1D profile of 2-D surface
  } else passes <- c(1L,2L,3L)

  for (pass in passes) {
    ptit <- 1L
    if (pass==3L) {
      xSeq <- rev(seq_len(length(xGrid)))
      ySeq <- rev(seq_len(length(yGrid)))
    } else {
      xSeq <- seq_len(length(xGrid))
      ySeq <- seq_len(length(yGrid))
    }
    for(i in xSeq) {
      for(j in ySeq) { ## NOTE no meaningful value outside the range covered by kriging
        fixedlist[fixed[1L]] <- xGrid[i]
        if (length(fixed)>1) {fixedlist[fixed[2]] <- yGrid[j]}
        locarglist <- list(fixedlist=fixedlist)
        if (length(fixedlist)==INFO$fittedparamnbr) { ## no profile
          if (INFO$fittedparamnbr==1L) {
            grmf <- list(full=unlist(fixedlist), value=purefn(unlist(fixedlist)), par=NULL, inKgspace=TRUE) ## single param, grid is the Kriged range
          } else { ## 2 params
            Lvalue <- purefn(unlist(fixedlist)) ## FR->FR: alsways use Kgtotal hull : valid ?
            grmf <- list(full=unlist(fixedlist), value=Lvalue, par=NULL, inKgspace= (! is.na(Lvalue)) )
          }
        } else {
          if (pass==1L) {
            toto <- system.time(grmf <- do.call(profileMethod,list(fixedlist=fixedlist))) ## returns canonical vector
            nProfileCalls <- nProfileCalls+1
            profilesTime <- profilesTime+toto["user.self"]
            if(profileMethod=="profileBySubHull") subhullinfo[[i, j]] <- grmf$subHull_V_and_H
          } else {
            ## second pass: S&W neighbors as starting point for maximisation
            ## third pass: N&E neighbors as starting point for maximisation
            # otherlist info
            locpts <- array(NA,c(2L,poslogL))
            if(pass==2L) {
              if (i>1) locpts[1, ] <- profpts[i-1, j, ]
              if (j>1) locpts[2, ] <- profpts[i, j-1, ]
            } else if (pass==3L) {
              if(i<length(xGrid)) locpts[1, ] <- profpts[i+1, j, ]
              if(j<length(yGrid)) locpts[2, ] <- profpts[i, j+1, ]
            }
            zut <- which.max(locpts[, poslogL])
            if (length(zut)==1L) { ## non-NA value in locpts...
              ## otherlist serves to reconstruct full initial vector for maximization (replacing rosglobal used by default)
              otherlist[othervarspos] <- locpts[zut, othervarspos] ## copie dans liste nomm['e]e
              locarglist$otherlist <- otherlist
              if (profileMethod=="profileBySubHull") locarglist$subHull_V_and_H <- subhullinfo[[i, j]]
              grmf <- do.call(profileMethod,locarglist) ## returns canonical vector
              previousvalue <- profpts[i, j, poslogL]
              ######## change 12/2011 (see additional comments in version < 29/01/2016)
              ### new code replaces even in new value is out of space and old one was in:
              if ( ! is.na(grmf$value) && (is.na(previousvalue) || grmf$value>previousvalue)) {
                profpts[i, j, ] <- c(grmf$par[othervars], grmf$value)
                inKrigSpace[i, j] <- grmf$inKgspace
              }
            } else grmf <- NULL
          }
        }
        if ( ! is.null(grmf)) { ## TRUE on 1st pass and if zut was TRUE otherwise
          profpts[i, j, ] <- c(grmf$par[othervars], grmf$value)
          ## = profiled out values which maximize profile | fixed, and profile logL
          if (pass==1L) { if ( ! is.na(grmf$value)) inKrigSpace[i, j] <- grmf$inKgspace}
          if (grmf$inKgspace && grmf$value>tempGridRelbestll) {
            tempGridRelbestll <- grmf$value
            tempBestfull <- grmf$full ## non log
          }
        }
        if (interactive()) {
          msg <- paste(" (step ", pass, ") Already ", ptit, " profile points computed out of ", ntotal, "     ", sep="")
          prevmsglength <- overcat(msg, prevmsglength)
          ptit <- ptit+1L
        }
      } ## y subloop
      if (pass==1L) {
        if ( ! interactive() && INFO$profile3passesBool && length(fixedlist)<INFO$fittedparamnbr) {
          pastHullNProf <- proc.time()["user.self"]-t0
          hullNProf <- pastHullNProf*length(xGrid)/i
          if (INFO$profile3passesBool && length(fixedlist)<INFO$fittedparamnbr) {
            if (length(fixedlist)>1) {estTime <- hullNProf+2*profilesTime} else {estTime <- hullNProf+profilesTime} ## overestim as next profile() calls should take less time
          } ##else (string must be short enough that the end of the line is not met...)
          msg <- paste("Estimated time for ", paste(usernames, collapse=", ") , " profile plot: ", signif(hullNProf, 2), " s. (remaining: ", signif(hullNProf-pastHullNProf, 2), " s.)      ", sep="")
          if (hullNProf>10) prevmsglength <- overcat(msg, prevmsglength)
        }
      }
    } ## x loop
  } ## end iteration on passes
  if (length(fixed)>1L) {
    locfn <- function(fixed,MARGIN,Grid,oGrid) {
      otherMARGIN <- 3L-MARGIN
      locarglist <- list(fixedlist=list(),otherlist=list())
      whichmaxS <- apply(profpts[,,poslogL,drop=FALSE],MARGIN, function(rowOrCol) {
        res <- which.max(rowOrCol)
        if (length(res)==0) res <- NA
        res
      })
      Seq <- seq_len(length(Grid))
      if (MARGIN==1L) {
        rc <- cbind(Seq,whichmaxS)
      } else {rc <- cbind(whichmaxS,Seq)}
      marginKgspace <- inKrigSpace[rc]
      margprof <- profpts[cbind(rc,poslogL)]
      oldmargprof <- .blackbox.data$options$margProfsInfo[[fixed[MARGIN]]]$margprof
      if (is.null(oldmargprof) || length(oldmargprof)!=length(margprof)) {
        oldmargprof <- rep(-Inf,length(margprof))
        oldinKgspace <- NULL
      } else {
        oldinKgspace <- .blackbox.data$options$margProfsInfo[[fixed[MARGIN]]]$inKgSpace
      }
      keepOld <- (margprof <= oldmargprof)
      keepOld[is.na(keepOld)] <- is.na(margprof)[is.na(keepOld)]
      marginKgspace[keepOld] <- oldinKgspace[keepOld] ## not even sure its used...
      margprof[keepOld] <- oldmargprof[keepOld]
      for (i in Seq[ ! keepOld]) {
        locarglist$fixedlist[fixed[MARGIN]] <- Grid[i]
        locarglist$otherlist[fixed[otherMARGIN]] <- oGrid[whichmaxS[i]] ## here otherMARGIN goes into otherlist, not into fixedlist
        locarglist$otherlist[othervars] <- profpts[rc[i,1L],,][rc[i,2L],othervarspos]
        # I cannot provide subHull_V_and_H because the available ones are of dim d-2, not d-1
        grmf <- do.call(profileMethod,locarglist)
        if ( all(fixed %in% INFO$FONKgNames) ## rosglobal in Krig space
             && grmf$value > margprof[i]) { ## should be the case, but...
          margprof[i] <- grmf$value
          marginKgspace[i] <- TRUE
          if (grmf$inKgspace && grmf$value>tempGridRelbestll) {
            tempGridRelbestll <<- grmf$value
            tempBestfull <<- grmf$full ## non log
          }
          #cat("+")
        } #else cat("-") ## and keep margprof[i]
      }
      res <- list(xGrid=Grid,#yGrid=NULL,
                  margprof=margprof,inKgSpace=marginKgspace)
      return(res)
    }
    margresu1 <- locfn(fixed=fixed,MARGIN=1L,Grid=xGrid,oGrid=yGrid)
    margresu2 <- locfn(fixed=fixed,MARGIN=2L,Grid=yGrid,oGrid=xGrid)
  }
  #resets maximum if a better point has been found
  if ( ! is.null(tempBestfull) ) { ## := a better point has been found in Kg space
    ## the arg of tofullKrigingspace must be in .blackbox.data$options$FONKgScale !!
    for(st in names(tempBestfull)) {if (islogscale(st)) {tempBestfull[st] <- log(tempBestfull[st])}}
    rosglobal <- findglobalMLE(initptinfK=tofullKrigingspace(tempBestfull))
    blackbox.options(rosglobal=rosglobal)
    returncode <- rosglobal$convergence
    tmp <- rosglobal$edgelevel
    if (tmp>0) returncode <- returncode+tmp/(10^ceiling(log(tmp, 10))) ## second summand goes in decimal part of returcode
    writeoutput(paste(blackbox.getOption("dataFile"), "(profiling)", sep=""), returncode=returncode, levelSlot=NA, CIloSlot=NA, CIupSlot=NA)
    do.call(blackbox.options, list(CIlo= NA, CIup=NA)) # LRT = NA, ?
  } ## end foundmore
  #computing relative profile Log Likelihoods:
  for(i in xSeq) {
    for(j in ySeq) { ## NOTE no meaningful value outside the range covered by kriging
      RelLik[i, j] <- exp(profpts[i, j, poslogL]-blackbox.getOption("rosglobal")$value)
    }
  }
  if (length(fixed)>1L) {
    margRelLik1 <- exp(margresu1$margprof-blackbox.getOption("rosglobal")$value)
    .blackbox.data$options$margProfsInfo[[fixed[1L]]] <- c(margresu1,list(RelLik=margRelLik1))
    margRelLik2 <- exp(margresu2$margprof-blackbox.getOption("rosglobal")$value)
    .blackbox.data$options$margProfsInfo[[fixed[2L]]] <- c(margresu2,list(RelLik=margRelLik2))
    ## ici l'acces direct aux membres de la liste est utile
  } else if (length(fixed)==1L) {
    .blackbox.data$options$margProfsInfo[[fixed[1L]]] <- list(xGrid=xGrid,#yGrid=NULL,
                                                          margprof=profpts[,, poslogL],RelLik=RelLik,inKgSpace=inKrigSpace[,1L])
  }
  ## noted the following problem 09/2011:
  ## the profile values are computed in a locchull. E.g. if kriging in Nb space and profile for 2Nmu, 2Nm.
  ## the predicted maxi in this space may be (by far) above rosglobal in kriging space. Then the profile plot has values way above 1
  ## We keep track of the info in inKrigSpace
  if (interactive()) {cat("\n")}
  return(list(xGrid=xGrid, yGrid=yGrid, RelLik=RelLik, inKrigSpace=inKrigSpace))
} ## end calcGridRelProfile
