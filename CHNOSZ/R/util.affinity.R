# CHNOSZ/util-affinity.R
# helper functions for affinity()

energy <- function(what,vars,vals,lims,T=get("thermo")$opt$Tr,P="Psat",IS=0,sout=NULL,exceed.Ttr=FALSE,transect=FALSE) {
  # 20090329 extracted from affinity() and made to
  # deal with >2 dimensions (variables)

  # calculate "what" property
  # logK - logK of formation reactions
  # logact.basis - logarithms of activities of basis species in reactions
  # logact.species - logarithms of activities of species in reactions
  # logQ - logQ of formation reactions
  # A - A/2.303RT of formation reactions
  # "vars": vector of variables (T,P,basisnames,IS)
  # "vals": list of values for each variable
  # "lims": used to get the dimensions (resolution for each variable)

  ### some argument checking
  if(length(unique(vars)) != length(vars)) stop("please supply unique variable names")
  ## basis definition / number of basis species 
  mybasis <- basis()
  nbasis <- nrow(mybasis)
  ## species definition / number of species
  myspecies <- get("thermo")$species
  if(is.character(what)) {
    if(is.null(myspecies)) stop('species properties requested, but species have not been defined')
    nspecies <- nrow(myspecies)
    if(!identical(what,"logact.basis")) ispecies <- 1:nspecies
  }
  ## the dimensions of our arrays
  resfun <- function(lim) lim[3]
  mydim <- sapply(lims,resfun)
  if(transect) {
    if(!isTRUE(all.equal(min(mydim),max(mydim)))) stop("variables define a transect but their lengths are not all equal")
    mydim <- mydim[[1]]
  }
  ## the number of dimensions we have
  if(transect) nd <- 1 else nd <- length(vars) # == length(mydim)
  ## basis names / which vars denote basis species
  basisnames <- rownames(mybasis)
  ibasisvar <- match(vars,basisnames)
  varisbasis <- !is.na(ibasisvar)
  ibasisvar <- ibasisvar[!is.na(ibasisvar)]
  ## which vars are in P,T,IS
  varissubcrt <- vars %in% c("P","T","IS")
  if(length(which(varissubcrt)) > 2) stop("sorry, currently only up to 2 of P,T,IS are supported")
  ## categorize the basis species:
  # 0 - not in the vars; 1 - one of the vars
  ibasis <- 1:nbasis
  ibasis0 <- ibasis[!ibasis %in% ibasisvar]
  ibasis1 <- ibasis[ibasis %in% ibasisvar]
  if(identical(what,"logact.basis")) ispecies <- ibasis
  ## what subcrt variable is used to make a 2-D grid?
  if(length(which(varissubcrt)) > 1 & !transect) {
    if("IS" %in% vars) grid <- "IS"
    else grid <- vars[varissubcrt][1]
  } else grid <- NULL
  ### done argument processing

  ### function to index the variables in a permuted order
  # by swapping ivar for the first
  # e.g. for nd=5, ivars(4)=c(4,2,3,1,5)
  ivars <- function(ivar,iv=NULL) {
    if(nd==0) return(1)
    if(is.null(iv)) iv <- 1:nd
    iv.1 <- iv[1]
    iv[1] <- ivar
    iv[ivar] <- iv.1
    return(iv)
  }

  ### functions for logact / logQ
  # a basis species not in var
  logact.basis0.fun <- function(ibasis) {
    logact <- mybasis$logact[ibasis]
    # for the case where this basis species is buffered
    if(!can.be.numeric(logact)) logact <- 0
    else logact <- as.numeric(logact)
    return(array(logact,mydim))
  }
  # a basis species in var
  logact.basis1.fun <- function(ivar) {
    dim.in <- dim(vals[[ivar]])
    if(is.null(dim.in)) dim.in <- 0
    # why doesn't this work?
    # if(identical(dim.in,mydim))
    if(all(dim.in==mydim) | transect) return(vals[[ivar]])
    else return(aperm(array(vals[[ivar]],mydim[ivars(ivar)]),ivars(ivar)))
  }
  # any basis species
  logact.basis.fun <- function(ibasis) {
    if(ibasis %in% ibasis0) return(logact.basis0.fun(ibasis))
    else return(logact.basis1.fun(match(basisnames[ibasis],vars)))
  }
  # all basis species
  logact.basis <- function() lapply(ibasis,logact.basis.fun)
  # logact of a single species
  logact.species.fun <- function(ispecies) array(myspecies$logact[ispecies],mydim)
  ## contributions to logQ
  # from a single basis species 
  logQ.basis.fun <- function(ibasis,coeff) - coeff * logact.basis.fun(ibasis)
  # from all basis species in a single formation reaction
  logQ.basis.species <- function(ispecies) 
    Reduce("+", mapply(logQ.basis.fun,ibasis,myspecies[ispecies,1:nbasis],SIMPLIFY=FALSE))
  # all basis species in all reactions
  logQ.basis <- function() mapply(logQ.basis.species,1:nspecies,SIMPLIFY=FALSE)
  # by a single species
  logQ.species.fun <- function(ispecies,coeff) coeff * logact.species.fun(ispecies)
  # by all species
  logQ.species <- function() 
    mapply(logQ.species.fun,1:nspecies,1,SIMPLIFY=FALSE)
  # total logQ of all reactions
  logQ <- function() lsum(logQ.basis(),logQ.species())

  ### function for calling subcrt
  sout.fun <- function(property="logK") {
    if(!is.null(sout)) return(sout) else {
      ## subcrt arguments
      species <- c(mybasis$ispecies,myspecies$ispecies)
      if("T" %in% vars) T <- vals[[which(vars=="T")]]
      if("P" %in% vars) P <- vals[[which(vars=="P")]]
      if("IS" %in% vars) IS <- vals[[which(vars=="IS")]]
      s.args <- list(species=species,property=property,T=T,P=P,IS=IS,grid=grid,convert=FALSE,exceed.Ttr=exceed.Ttr)
      return(do.call("subcrt",s.args)$out)
    }
  }

  ### functions for logK/subcrt props
  # the logK contribution by any species or basis species
  X.species <- function(ispecies,coeff,X) coeff * sout[[ispecies]][,names(sout[[ispecies]])==X]
  # the logK contribution by all basis species in a reaction
  X.basis <- function(ispecies,X) Reduce("+", mapply(X.species,ibasis,-myspecies[ispecies,ibasis],X,SIMPLIFY=FALSE))
  # the logK of any reaction
  X.reaction <- function(ispecies,X) X.species((ispecies+nbasis),1,X) + X.basis(ispecies,X)
  # to get logK or subcrt props or other values into the dimensions we are using
  dim.fun <- function(x,idim=NULL) {
    if(is.null(idim)) {
      if(transect) idim <- 1
      else if(is.null(grid)) {
        # one of T,P,IS
        ivar <- which(vars %in% c("T","P","IS"))
        if(length(ivar)==0) ivar <- 1
        idim <- ivars(ivar)
      } else {
        # two of T,P,IS
        ivar1 <- which(varissubcrt)[1]
        ivar2 <- which(varissubcrt)[2]
        idim <- ivars(ivar2,ivars(ivar1))
      }
    }
    return(aperm(array(x,mydim[idim]),idim))
  }
  # properties of all species
  X.fun <- function(X) lapply(lapply(ispecies,X.reaction,X),dim.fun)
  logK <- function() lapply(ispecies,X.reaction,"logK")
  # A/2.303RT
  A <- function() {
    out <- lsub(X.fun("logK"),logQ())
    # deal with affinities of protein ionization here 20120527
    if("H+" %in% rownames(mybasis)) {
      # which species are proteins
      isprotein <- grepl("_", myspecies$name)
      if(any(isprotein)) {
        # the rownumbers in thermo$protein
        ip <- iprotein(myspecies$name[isprotein])
        # get the affinity of ionization
        iHplus <- match("H+", rownames(mybasis))
        # as.numeric is needed in case the logact column is character mode
        # due to buffer definition (that means we can't do pH buffers)
        pH <- -as.numeric(mybasis$logact[iHplus])
        A.ionization <- A.ionization(ip, vars, vals, T=T, P=P, pH=pH, transect=transect)
        # add it to the affinities of formation reactions of the non-ionized proteins
        out[isprotein] <- lsum(out[isprotein], A.ionization)
      }
    }
    return(out)
  }

  ### call the necessary functions
  if(!is.character(what)) {
    # expand numeric values into our dimensions
    # (used by energy.args() for calculating pe=f(Eh,T) )
    # TODO: document that sout here denotes the dimension
    # we're expanding into
    return(dim.fun(what,ivars(sout)))
  } else if(what %in% c('G','H','S','Cp','V','E','kT','logK')) {
    # get subcrt properties for reactions
    sout <- sout.fun(what)
    a <- X.fun(what)
  } else if(length(agrep("species",what)) > 0) {
    # get subcrt properties for species
    # e.g. what=G.species, Cp.species etc
    mywhat <- s2c(what,sep=".",keep.sep=FALSE)[1]
    sout <- sout.fun(mywhat)
    a <- lapply(mapply(X.species,ispecies+nbasis,1,mywhat,SIMPLIFY=FALSE),dim.fun)
  } else {
    # get affinities or logact.basis, logact.species or logQ
    if(what=="A") sout <- sout.fun("logK")
    a <- do.call(what,list())
  }

  ### use species indices as the names 
  if(what=="logact.basis") names(a) <- mybasis$ispecies
  else names(a) <- myspecies$ispecies
  ### return the results
  return(list(sout=sout,a=a))
}

energy.args <- function(args) {
  ## extracted from affinity() and modified 20090329 jmd
  # takes a list of arguments which define the dimensions
  # over which to calculate logQ, logK and affinity
  # the names should be T, P, IS and names of basis species
  # (or pH, pe, Eh)
  thermo <- get("thermo")
  ## inputs are like c(T1,T2,res)
  # and outputs are like seq(T1,T2,length.out=res)
  # unless transect: do the variables specify a transect? 20090627
  transect <- any(sapply(args,length) > 3)
  ## convert T, P args and take care of 
  # single values for T, P or IS
  T <- thermo$opt$Tr
  P <- "Psat"
  IS <- 0
  T.is.var <- P.is.var <- IS.is.var <- FALSE
  arg.is.T <- names(args)=="T"
  arg.is.P <- names(args)=="P"
  arg.is.IS <- names(args)=="IS"
  if(any(arg.is.T)) {
    T <- args[[which(arg.is.T)]]
    if(length(T) > 1) T.is.var <- TRUE
    if(transect) args[[which(arg.is.T)]] <- T <- envert(T,'K')
    else args[[which(arg.is.T)]][1:2] <- T[1:2] <- envert(T[1:2],'K')
    T <- T[!is.na(T)]
  }
  if(any(arg.is.P)) {
    P <- args[[which(arg.is.P)]]
    if(length(P) > 1) P.is.var <- TRUE
    if(!identical(P,"Psat")) {
      if(transect) args[[which(arg.is.P)]] <- P <- envert(P,'bar')
      else args[[which(arg.is.P)]][1:2] <- P[1:2] <- envert(P[1:2],'bar')
    }
    P <- P[!is.na(P)]
  }
  if(any(arg.is.IS)) {
    IS <- args[[which(arg.is.IS)]]
    if(length(IS) > 1) IS.is.var <- TRUE
  }
  # report non-variables to user
  if(!T.is.var)
    msgout('energy.args: temperature is ',outvert(T,'K'),' ',T.units(),'\n')
  if(!P.is.var) {
    if(identical(P,"Psat")) msgout("energy.args: pressure is Psat\n")
    else msgout('energy.args: pressure is ',outvert(P,'bar'),' ',P.units(),'\n')
  }
  if(!IS.is.var & !identical(IS,0)) msgout('energy.args: ionic strength is ',IS,'\n')
  # default values for resolution
  res <- 128
  # where we store the output
  what <- "A"
  vars <- character()
  vals <- list(NA)
  # this needs to have 1 as the third component b/c
  # energy() uses it to build an array with the given dimension
  lims <- list(c(NA, NA, 1))
  # clean out non-variables
  if(any(arg.is.T) & !T.is.var) args <- args[names(args)!="T"]
  if(any(arg.is.P) & !P.is.var) args <- args[names(args)!="P"]
  if(any(arg.is.IS) & !IS.is.var) args <- args[names(args)!="IS"]
  # the property we're interested in
  if("what" %in% names(args)) {
    what <- args[[names(args)=="what"]]
    args <- args[names(args)!="what"]
  }
  # assemble the variables
  if(length(args) > 0) {
    for(i in 1:length(args)) {
      nametxt <- names(args)[i]
      if(transect) lims.orig <- c(min(args[[i]]),max(args[[i]]))
      else lims.orig <- args[[i]][1:2]
      if(names(args)[i]=="pH") {
        names(args)[i] <- "H+"
        if(transect) args[[i]] <- -args[[i]]
        else args[[i]][1:2] <- -args[[i]][1:2]
        if(!'H+' %in% rownames(thermo$basis)) 
          msgout('energy.args: pH requested, but no H+ in the basis\n')
      } 
      if(names(args)[i]=="pe") {
        names(args)[i] <- "e-"
        if(!'e-' %in% rownames(thermo$basis)) 
          msgout('energy.args: pe requested, but no e- in the basis\n')
        if(transect) args[[i]] <- -args[[i]]
        else args[[i]][1:2] <- -args[[i]][1:2]
      }
      if(length(args[[i]]) < 3 & !transect) args[[i]] <- c(args[[i]],res)
      vars[length(vars)+1] <- names(args)[i]
      if(transect) {
        vals[[length(vars)]] <- args[[i]]
        lims[[length(vars)]] <- c(lims.orig,length(vals[[i]]))
      } else {
        vals[[length(vars)]] <- seq(args[[i]][1],args[[i]][2],length.out=args[[i]][3])
        lims[[length(vars)]] <- args[[i]]
      }
      names(lims)[length(vars)] <- names(args)[i]
      # say something about the identities, ranges, and units of the variables
      unittxt <- ""
      # number of values
      if(transect) n <- length(args[[i]]) else n <- args[[i]][3]
      # physical state
      ibasis <- match(nametxt, rownames(thermo$basis))
      if(isTRUE(as.logical(ibasis))) {
        if(thermo$basis$state[ibasis]=="gas") nametxt <- paste("log_f(", nametxt, ")", sep="") 
        else nametxt <- paste("log_a(", nametxt, ")", sep="") 
      }
      # temperature and pressure and Eh
      if(nametxt=="T") unittxt <- " K"
      if(nametxt=="P") unittxt <- " bar"
      if(nametxt=="Eh") unittxt <- " V"
      msgout("energy.args: variable ", length(vars), " is ", nametxt, 
        " at ", n, " values from ", lims.orig[1], " to ", lims.orig[2], unittxt, "\n")
    }
  }
  args <- list(what=what,vars=vars,vals=vals,lims=lims,T=T,P=P,IS=IS,transect=transect)

  # convert Eh to pe
  if("Eh" %in% args$vars) {
    # get Eh into our dimensions
    Eh.args <- args
    # what variable is Eh
    Eh.var <- which(args$vars=="Eh")
    Eh.args$what <- args$vals[[Eh.var]]
    Eh.args$sout <- Eh.var
    Eh <- do.call("energy",Eh.args)
    # get temperature into our dimensions
    T.args <- args  
    if("T" %in% args$vars) {
      T.var <- which(args$vars=="T")
      T.args$what <- args$vals[[T.var]]
    } else {
      T.var <- 1
      T.args$what <- T
    }
    T.args$sout <- T.var
    T <- do.call("energy",T.args)
    # do the conversion on vectors
    mydim <- dim(Eh)
    Eh <- as.vector(Eh)
    T <- as.vector(T)
    pe <- convert(Eh,"pe",T=T)
    dim(pe) <- mydim
    # update the arguments list
    args$vars[Eh.var] <- "e-"
    args$vals[[Eh.var]] <- -pe
  }
  return(args)
}

slice.affinity <- function(affinity,d=1,i=1) {
  # take a slice of affinity along one dimension
  a <- affinity
  for(j in 1:length(a$values)) {
    # preserve the dimensions (especially: names(mydim))
    # - fix for change in behavior of aperm in R-devel 2015-11-17
    mydim <- dim(a$values[[j]])
    a$values[[j]] <- as.array(slice(a$values[[j]],d=d,i=i))
    # the dimension from which we take the slice vanishes
    dim(a$values[[j]]) <- mydim[-d]
  }
  return(a)
}

A.ionization <- function(iprotein, vars, vals, T=get("thermo")$opt$Tr, P="Psat", pH=7, transect=FALSE) {
  # a function to build a list of values of A/2.303RT of protein ionization
  # that can be used by energy(); 20120527 jmd
  # some of the variables might not affect the values (e.g. logfO2)
  # what are the variables that affect the values
  T <- convert(T, "C")
  if(!is.na(iT <- match("T", vars))) T <- convert(vals[[iT]], "C")
  if(!is.na(iP <- match("P", vars))) P <- vals[[iP]]
  if(!is.na(iHplus <- match("H+", vars))) pH <- -vals[[iHplus]]
  # is it a transect (single points) or a grid?
  if(transect) TPpH <- list(T=T, P=P, pH=pH)
  else {
    # make a grid of all combinations
    # put T, P, pH in same order as they show up vars
    egargs <- list(T=T, P=P, pH=pH)
    # order(c(NA, NA, NA)) might segfault in some versions of R (seen in R 2.15.3 on Linux)
    if(is.na(iT) & is.na(iP) & is.na(iHplus)) TPpHorder <- c(1, 2, 3)
    else TPpHorder <- order(c(iT, iP, iHplus))
    egargs <- c(egargs[TPpHorder], list(stringsAsFactors=FALSE))
    TPpH <- do.call(expand.grid, egargs)
    # figure out the dimensions of T-P-pH (making sure to drop any that aren't in vars)
    TPpHdim <- numeric(3)
    TPpHdim[iT] <- length(T)
    TPpHdim[iP] <- length(P)
    TPpHdim[iHplus] <- length(pH)
    TPpHdim <- TPpHdim[!TPpHdim==0]
    # figure out the dimensions of the other vars
    othervars <- vars[!vars %in% c("T", "P", "H+")]
    iother <- match(othervars, vars)
    otherdim <- sapply(vals, length)[iother]
    # if Eh was given to energy.args, values of pe were calculated in all dimensions
    # figure out the original length of the Eh variable
    ieminus <- match("e-", vars[iother])
    if(!is.na(ieminus)) {
      otherdim[ieminus] <- NA
      edim <- dim(vals[[iother[ieminus]]])
      # loop TPpHdim and otherdim, taking out each one from edim
      for(dim in c(TPpHdim, otherdim)) {
        id <- match(dim, edim)
        edim[id] <- NA
      }
      otherdim[ieminus] <- edim[!is.na(edim)]
    }
    # the permutation vector
    if(length(iother) > 0) allvars <- c(vars[-iother], vars[iother])
    else allvars <- vars
    perm <- match(vars, allvars)
  }
  # initialize output list
  out <- vector("list", length(iprotein))
  # get aa from iprotein
  aa <- ip2aa(iprotein)
  # calculate the values of A/2.303RT as a function of T-P-pH
  A <- ionize.aa(aa=aa, property="A", T=TPpH$T, P=TPpH$P, pH=TPpH$pH)
  if(transect) {
    # just turn the values into a list
    for(i in 1:length(iprotein)) out[[i]] <- A[, i]
  } else for(i in 1:length(iprotein)) {
    # what are the values of ionization affinity for this protein
    thisA <- A[, i]
    # apply the dimensions of T-P-pH
    tpphdim <- TPpHdim
    if(length(tpphdim)==0) tpphdim <- 1
    thisA <- array(thisA, tpphdim)
    # grow into the dimensions of all vars
    alldim <- c(TPpHdim, otherdim)
    if(length(alldim)==0) alldim <- 1
    thisA <- array(thisA, alldim)
    # now permute the array to put dimensions in same order as the variables
    thisA <- aperm(thisA, perm)
    # store in output list
    out[[i]] <- thisA
  }
  return(out)
}
