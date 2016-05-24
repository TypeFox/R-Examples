powcomp.fast <- function(law.indices,stat.indices,vectn=c(20,50,100),M=10^3,levels=c(0.05,0.1),critval=NULL,alter=create.alter(stat.indices),parlaws=NULL,parstats=NULL,nbclus=1,model=NULL,null.law.index=2,null.law.pars=NULL,Rlaws=NULL,Rstats=NULL,center=FALSE,scale=FALSE) {

    compquant <- 0L # This argument should only be used when one wants to compute the critical values (quantiles) of some test statistics,
                    # through the function many.crit(). The values of the statistics will be returned through the critvalL argument (of size M*vectnlen*statslen*lawlen)
    
  if (any(stat.indices == 0) & is.null(Rstats)) stop("'Rstats' should be a list whose components are R functions.")
  if (any(stat.indices == 0)) {
    if (!is.list(Rstats)) stop("'Rstats' should be a list whose components are R functions.")
    for (i in 1:length(stat.indices)) if ((stat.indices[i] == 0) & !is.function(Rstats[[i]])) stop(paste("The ",i,"th component of 'Rstats' should be an R function",sep=""))
  }
  if(is.null(Rstats)) Rstats <- list(NULL)

  
    Rcpp <- any(law.indices == 0)

    if (Rcpp) {
        if (length(Rlaws) != length(law.indices)) stop("When some law indices in 'law.indices' are equal to 0, this means that you will be using some R random generators. In that case, you should provide the names of the random generation functions in the corresponding components of 'Rlaws' list, the other components should be set to NULL.")
      tmp <- gsub(" ","",paste(text=match.call()$Rlaws))[-1]
        tmpnames <- tmp
      for (i in 1:length(law.indices)) {
          if (!is.null(Rlaws[[i]])) {
              if (class(Rlaws[[i]]) != "function") stop("Each non-null compoment of the list 'Rlaws' should be a (random generation) R function.")
              if (law.indices[i] != 0) stop(paste("law.indices[",i,"] should be set to 0.",sep=""))
          } else {
              Rlaws[[i]] <- function(){}
          }
      }
  }

  
  if(getRversion() < "3.1.0") dontCheck <- identity

  if (nbclus>1) {
#    suppressWarnings(parallel.pkg.present <- require(parallel))
    parallel.pkg.present <- "package:parallel" %in% search()
    Rmpi.pkg.present <- "package:Rmpi" %in% search()
    if (all(!c(parallel.pkg.present,Rmpi.pkg.present))) stop("Either package parallel or Rmpi should be installed!")
#    suppressWarnings(rsprng.pkg.present <- require(rsprng))
#    if (!rsprng.pkg.present) stop("Package rsprng is not installed!")
    cluster.type <- if (parallel.pkg.present) "PSOCK" else "MPI" # We prefer to use "PSOCK" (i.e. parallel) because it's easier.
  }
  
  vectn.len <- length(vectn)
  stats.len <- length(stat.indices)
  laws.len  <- length(law.indices)
  nblevel   <- length(levels)  
  
# Management of critval and creation of critvalL, critvalR  and usecrit
# If we provide a single value in critval$statj then it is critvalR
# If we provide two values in critval$statj then it is c(critvalL,critvalR)  ... IN THAT ORDER!!
  critvalL <- critvalR <- rep(0,vectn.len*stats.len*nblevel)
  usecrit <- rep(0,vectn.len*stats.len)
  if (is.null(critval)) {
      warning(paste("'critval' has been computed internally using function many.crit() with the value of 'law.index'=",null.law.index," (i.e. ",law.cstr(null.law.index)$name,")",sep=""))
      if (is.null(null.law.pars)) {
          tmp2 <- law.cstr(null.law.index)$law.pars
          null.law.pars <- c(tmp2,rep(0.0,4-length(tmp2)))
      }
    critval <- many.crit(law.index=null.law.index,stat.indices,M,vectn,levels,alter,null.law.pars,parstats,center=center,scale=scale)
  }
  if (!is.list(critval)) stop("'critval' should be a list")
  if (is.null(names(critval))) stop("'critval' should be a named list")
  if (any(is.na(names(critval)))) stop("'critval' names should all be defined")
  if (length(critval) != length(stat.indices)) stop("'critval' and 'stat.indices' should have the same length")
  for (s in 1:stats.len) {
    crittmp <- as.matrix(critval[[s]][,-3])
    if (ncol(crittmp) == 1) crittmp <- t(crittmp)
    for (l in 1:nblevel) {
      vals <- crittmp[crittmp[,"level"]==levels[l],]
      for (k in 1:vectn.len) {
        if (vectn.len == 1) {
          vals2 <- vals[-(1:2)]
        } else {
          vals2 <- vals[vals[,"n"]==vectn[k],-(1:2)]
        }
        if (!is.null(vals2)) {
          usecrit[k+vectn.len*(s-1)] <- 1 # I CHANGED HERE!! IF THERE IS A BUG, check here ...
          critvalL[k+vectn.len*(l-1)+nblevel*vectn.len*(s-1)] <- vals2[1]
          critvalR[k+vectn.len*(l-1)+nblevel*vectn.len*(s-1)] <- vals2[2]
        }
      }      
    }
  }
  
# Management of alter
  if (!is.null(alter)) {
    if (!is.list(alter)) stop("'alter' should be a list")
    if (is.null(names(alter))) stop("'alter' should be a named list")
    if (any(is.na(names(alter)))) stop("'alter' names should all be defined")
    if (length(alter) != length(stat.indices)) stop("'alter' and 'stat.indices' should have the same length")
    for (s in 1:stats.len) {
      if (stat.indices[s] != 0) {
        if (names(alter)[s] != paste("stat",stat.indices[s],sep="")) stop(paste("Name of 'alter'[[",s,"]] should be equal to 'stat",stat.indices[s],sep=""))
        if (!(alter[[s]] %in% 0:4)) stop(paste("'alter'[[",s,"]] should be in  {0,1,2,3,4}.",sep=""))
        Cstat.name <- paste("stat",as.character(stat.indices[s]),sep="")
        alter.true <- .C(dontCheck(Cstat.name),as.double(0.0),1L,0.05,1L,rep(" ",50),1L,0.0,0L,0.0,0.0,0.0,0L,alter=as.integer(alter[s]),0L,rep(0.0,4),0L,PACKAGE="PoweR")$alter
        if (alter[[s]] != alter.true) {
          warning(paste("'alter'[[",s,"]] should be set to ",alter.true,". We have done this for you!"),sep="")
          alter[[s]] <- alter.true
        }
      }
    }
    alter <- unlist(alter)
  } else { # alter is NULL
    alter <- rep(0,stats.len)
    names(alter) <- paste("stat",stat.indices,sep="")
  }
 

# Management of parlaws
  parlawtmp <- nbparlaws <- c()
  if (!is.null(parlaws)) {
    if (!is.list(parlaws)) stop("'parlaws' should be a list")
    if (is.null(names(parlaws))) stop("'parlaws' should be a named list")
    if (any(is.na(names(parlaws)))) stop("'parlaws' names should all be defined")
    if (length(parlaws) != length(law.indices)) stop("'parlaws' and 'law.indices' should have the same length")
    for (s in 1:laws.len) {
        if (law.indices[s] != 0) {
            if (names(parlaws)[s] != paste("law",law.indices[s],sep="")) stop(paste("Name of 'parlaws'[[",s,"]] should be equal to 'law",law.indices[s],sep=""))
            if (length(parlaws[[s]]) > 4) stop(paste("Length of 'parlaws'[[",s,"]] should not exceed 4.",sep=""))
            if ((length(parlaws[[s]]) > 1) && any(is.na(parlaws[[s]]))) stop(paste("'parlaws'[[",s,"]] cannot contain NA values unless its length is 1.",sep=""))
            nbparlaws <- c(nbparlaws,length(na.omit(parlaws[[s]])))
            parlawtmp <- c(parlawtmp,c(parlaws[[s]],rep(0.0,4-nbparlaws[s])))
        } else {
            if (!all(is.na(parlaws[[s]]))) {
                npartmp <- length(unlist(formals(eval(parse(text=tmp[s])))[-1]))
                if (sum(!is.na(parlaws[[s]])) !=  npartmp) stop(paste("The number of law parameters set for ",i,"th component of 'parlaws' should be ",npartmp,sep=""))
                nbparlaws <- c(nbparlaws,npartmp)
                parlawtmp <- c(parlawtmp,c(parlaws[[s]],rep(0.0,4-nbparlaws[s])))
            }
        }
    }
} else {
    for (s in 1:laws.len) {
        if (law.indices[s] != 0) {
            tmp <- law.cstr(law.indices[s])
            nbparlaws <- c(nbparlaws,tmp$nbparams)
            parlawtmp <- c(parlawtmp,c(tmp$law.pars,rep(0.0,4-nbparlaws[s])))
        } else {
            npartmp <- length(unlist(formals(eval(parse(text=tmp[s])))[-1]))
            nbparlaws <- c(nbparlaws,npartmp)
            parlawtmp <- c(parlawtmp,c(unlist(formals(eval(parse(text=tmp[s])))[-1]),rep(0.0,4-nbparlaws[s])))
        }
    }
}
  parlaws <- parlawtmp

  
# Management of parstats     pas tres bien gere quand stat.indices[s] = 0 ? .....
  nbparstats <- rep(0, length(stat.indices))
  nbparstats[stat.indices != 0] <- getnbparstats(stat.indices[stat.indices != 0])
  parstatstmp <- c()
  if (!is.null(parstats)) {
    if (!is.list(parstats)) stop("'parstats' should be a list")
    if (is.null(names(parstats))) stop("'parstats' should be a named list")
    if (any(is.na(names(parstats)))) stop("'parstats' names should all be defined")
    if (length(parstats) != length(stat.indices)) stop("'parstats' and 'stat.indices' should have the same length")
    for (s in 1:stats.len) {
      if (stat.indices[s] != 0) {
        if (names(parstats)[s] != paste("stat", stat.indices[s], sep="")) stop(paste("Name of 'parstats'[[", s, "]] should be equal to 'stat", stat.indices[s], sep = ""))
        if (!is.na(parstats[[s]]) && (nbparstats[s] == 0)) stop(paste("'parstats[['", s, "]] should be equal to NA", sep = ""))
        if ((nbparstats[s] != 0) && (length(parstats[[s]]) != nbparstats[s])) stop(paste("The length of parstats[[", s, "]] should be ", nbparstats[s], sep = ""))
        parstatstmp <- c(parstatstmp, parstats[[s]])
      }
    }
  } else {
    for (s in 1:stats.len) {
      if (stat.indices[s] != 0) {parstatstmp <- c(parstatstmp, stat.cstr(stat.indices[s])$stat.pars)}
    }
  }
  parstats <- parstatstmp
#  parstats[is.na(parstats)] <- 0
  parstats <- parstats[!is.na(parstats)]


  
# Management of model
  ## FAUDRA REGARDER CA DE PLUS PRES QUAND JE RENDRAI MODEL FONCTIONNEL!!  
  if (is.double(model) || is.integer(model)) {
    modelnum <- model
    funclist <- list(function(){})
    thetavec <- 0
    xvec <- 0
    p <- length(thetavec)
    np <- length(xvec)
  } else {
    if (is.null(model)) {
      modelnum <- 1
      funclist <- list(function(){})
      thetavec <- 0
      xvec <- 0
      p <- length(thetavec)
      np <- length(xvec)
    } else { # model should be a list (function(x,thetavec,xvec),theta,xvec)
      modelnum <- 0
      funclist <- list(model[[1]])
      thetavec <- model[[2]]
      xvec <- model[[3]]
      p <- length(thetavec)
      np <- length(xvec)     
    }
  }



# We perform the computations

  decision.len <- stats.len*vectn.len*laws.len*nblevel
  decision <- rep(0,decision.len)

    if (Rcpp | any(stat.indices == 0)) {
        if (nbclus > 1) { # We start the cluster
    # makeCluster = Create a set of copies of R running in parallel and communicating over sockets or using MPI.
            cl <- parallel::makeCluster(nbclus, type = cluster.type)		
#    clusterSetupSPRNG(cl)
                                        
      
            myfunc <- function(M) {
                require(PoweR)
                .Call("powcompfastRcpp",M=as.integer(M),law.indices=as.integer(law.indices),laws.len=as.integer(laws.len),vectn=as.integer(vectn),vectn.len=as.integer(vectn.len),stat.indices=as.integer(stat.indices),stats.len=as.integer(stats.len),decision=as.integer(decision),decision.len=as.integer(decision.len),levels=as.double(levels),nblevel=as.integer(nblevel),
                   cL=as.double(critvalL),cR=as.double(critvalR),usecrit=as.integer(usecrit),alter=as.integer(alter),nbparlaws=as.integer(nbparlaws),parlaws=as.double(parlaws),nbparstats=as.integer(nbparstats),parstats=as.double(parstats),as.integer(modelnum), funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),as.list(Rlaws),Rstats,as.integer(center), as.integer(scale), as.integer(compquant), PACKAGE="PoweR",NAOK=TRUE)
            }
    
            out <- parallel::clusterCall(cl, myfunc, round(M/nbclus)) # M/nbclus iterations are performed on each core
          
      # We stop the cluster
            parallel::stopCluster(cl)
    
        } else {

  #b) or without a cluster

            out <- list(.Call("powcompfastRcpp",M=as.integer(M),law.indices=as.integer(law.indices),laws.len=as.integer(laws.len),vectn=as.integer(vectn),vectn.len=as.integer(vectn.len),stat.indices=as.integer(stat.indices),stats.len=as.integer(stats.len),decision=as.integer(decision),decision.len=as.integer(decision.len),levels=as.double(levels),nblevel=as.integer(nblevel),
                           cL=as.double(critvalL),cR=as.double(critvalR),usecrit=as.integer(usecrit),alter=as.integer(alter),nbparlaws=as.integer(nbparlaws),parlaws=as.double(parlaws),nbparstats=as.integer(nbparstats),parstats=as.double(parstats),as.integer(modelnum), funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),as.list(Rlaws),Rstats,as.integer(center), as.integer(scale), as.integer(compquant), PACKAGE="PoweR",NAOK=TRUE))
            
        }

    } else {
  # a) Using a cluster
        if (nbclus > 1) { # We start the cluster
    # makeCluster = Create a set of copies of R running in parallel and communicating over sockets or using MPI.
            cl <- parallel::makeCluster(nbclus, type = cluster.type)		
#    clusterSetupSPRNG(cl)
                                        
      
            myfunc <- function(M) {
                require(PoweR)
                .C("powcompfast",M=as.integer(M),law.indices=as.integer(law.indices),laws.len=as.integer(laws.len),vectn=as.integer(vectn),vectn.len=as.integer(vectn.len),stat.indices=as.integer(stat.indices),stats.len=as.integer(stats.len),decision=as.integer(decision),decision.len=as.integer(decision.len),levels=as.double(levels),nblevel=as.integer(nblevel),
                   cL=as.double(critvalL),cR=as.double(critvalR),usecrit=as.integer(usecrit),alter=as.integer(alter),nbparlaws=as.integer(nbparlaws),parlaws=as.double(parlaws),nbparstats=as.integer(nbparstats),parstats=as.double(parstats),as.integer(modelnum), funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),as.integer(center), as.integer(scale), as.integer(compquant) ,PACKAGE="PoweR",NAOK=TRUE)
            }
    
            out <- parallel::clusterCall(cl, myfunc, round(M/nbclus)) # M/nbclus iterations are performed on each core
          
      # We stop the cluster
            parallel::stopCluster(cl)
    
        } else {

  #b) or without a cluster

            out <- list(.C("powcompfast",M=as.integer(M),law.indices=as.integer(law.indices),laws.len=as.integer(laws.len),vectn=as.integer(vectn),vectn.len=as.integer(vectn.len),stat.indices=as.integer(stat.indices),stats.len=as.integer(stats.len),decision=as.integer(decision),decision.len=as.integer(decision.len),levels=as.double(levels),nblevel=as.integer(nblevel),
                           cL=as.double(critvalL),cR=as.double(critvalR),usecrit=as.integer(usecrit),alter=as.integer(alter),nbparlaws=as.integer(nbparlaws),parlaws=as.double(parlaws),nbparstats=as.integer(nbparstats),parstats=as.double(parstats),as.integer(modelnum), funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),as.integer(center),as.integer(scale), as.integer(compquant), PACKAGE="PoweR",NAOK=TRUE))
            
        }
    }


    out[[1]] <- out[[1]][c("M","law.indices","vectn","stat.indices","decision","levels","cL","cR","usecrit","alter","nbparlaws","parlaws","nbparstats","parstats")]

  if (nbclus > 1) {
      out[[1]]$M <- nbclus * round(M / nbclus)
    for (clus in 2:nbclus) {
      out[[1]]$decision <- out[[1]]$decision + out[[clus]]$decision
    }
  }

    out[[1]]$nbclus <- nbclus
    out[[1]]$Rlaws <- Rlaws
    if (Rcpp) names(out[[1]]$Rlaws) <- tmpnames
  
  k <- 1
  for (i in 1:length(out[[1]]$nbparstats)) {
    if (out[[1]]$nbparstats[i] == 0) {out[[1]]$parstats[k]<- NA ; k <- k + 1} else k <- k + out[[1]]$nbparstats[i]
  }

  return(structure(out[[1]], class = c("power","list")))
  
}

