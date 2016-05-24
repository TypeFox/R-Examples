###{{{ multigroup

##' @export
multigroup <- function(models, datasets, fix, exo.fix=TRUE, keep=NULL, missing=FALSE, ...) {
  nm <- length(models)
  if (nm!=length(datasets)) stop("Supply dataset for each model")
  if (nm<2) stop("Two or more groups neeeded")
  mynames <- names(models)

  ## Check for random slopes
  xfix <- list()
  for (i in seq_len(nm)) {
    x0 <- models[[i]]
    data0 <- datasets[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x0,exo=TRUE))]
    xfix <- c(xfix, list(xfix0))
  }
  if (missing(fix)) {
    fix <- !any(unlist(lapply(xfix, function(x) length(x)>0)))
  }

  for (i in seq_len(nm)) {
    x0 <- models[[i]]
    data0 <- datasets[[i]]
    if (length(exogenous(x0)>0)) {
      catx <- categorical2dummy(x0,data0)
      models[[i]] <- catx$x; datasets[[i]] <- catx$data
    }
    if (!lava.options()$exogenous) exogenous(models[[i]]) <- NULL
  }

  models.orig <- NULL
######################
### MLE with MAR mechanism
######################
  if (missing) {

    parcount <- 0
    reservedpars <- c()
    mynpar <- c()
    for (i in seq_len(nm)) {
      ## Fix some parameters (predictors,latent variables,...)

      d0 <- datasets[[i]][1,,drop=FALSE]; d0[,] <- 1
      if (fix)
        models[[i]] <- fixsome(models[[i]], exo.fix=exo.fix, measurement.fix=fix, data=d0)
      ## Find named/labelled parameters
      rpar <- unique(parlabels(models[[i]]))
      reservedpars <- c(reservedpars, rpar)
      mynpar <- c(mynpar, with(index(models[[1]]), npar+npar.mean+npar.ex))
    }; reservedpars <- unique(reservedpars)
    nonamepar <- sum(mynpar)
    ## Find unique parameter-names for all parameters
    newpars <- c()
    i <- 0
    pos <- 1
    while(pos<=nonamepar) {
      i <- i+1
      newname <- paste0("par",i)
      if (!(newname%in%reservedpars)) {
        newpars <- c(newpars,newname)
        pos <- pos+1
      }
    }

    pos <- 0
    models0 <- list()
    datasets0 <- list()
    complidx <- c()
    nmodels <- 0
    modelclass <- c()
    nmis <- c()
    for (i in seq_len(nm)) {
      myvars <- unlist(intersect(colnames(datasets[[i]]),c(vars(models[[i]]),xfix[[i]],keep)))
      mydata <- datasets[[i]][,myvars]
      if (any(is.na(mydata))) {
        if (i>1) pos <- pos+mynpar[i-1]
        models[[i]] <- baptize(models[[i]],newpars[pos+seq_len(mynpar[i])] ,overwrite=FALSE)
        val <- missingModel(models[[i]],mydata,fix=FALSE,keep=keep,...)
        nmodels <- c(nmodels,length(val$models))
        complidx <- c(complidx,val$pattern.allcomp+nmodels[i]+1)
        nmis0 <- rowSums(val$patterns);
        allmis <- which(nmis0==ncol(val$patterns))
        if (length(allmis)>0) nmis0 <- nmis0[-allmis]
        nmis <- c(nmis,nmis0)
        datasets0 <- c(datasets0, val$datasets)
        models0 <- c(models0, val$models)
        modelclass <- c(modelclass,rep(i,length(val$models)))
      } else {
        datasets0 <- c(datasets0, list(mydata))
        models0 <- c(models0, list(models[[i]]))
        modelclass <- c(modelclass,i)
        nmis <- c(nmis,0)
      }
    }

    models.orig <- models

    suppressWarnings(val <- multigroup(models0,datasets0,fix=FALSE,missing=FALSE,exo.fix=TRUE,...))
    val$models.orig <- models.orig; val$missing <- TRUE
    val$complete <- complidx-1
    val$mnames <- mynames
    attributes(val)$modelclass <- modelclass
    attributes(val)$nmis <- nmis
    return(val)
  }


######################
### Usual analysis:
######################
  warned <- FALSE
  for (i in seq_len(nm)) {
    if (inherits(datasets[[i]],c("data.frame","matrix"))) {
      myvars <- intersect(colnames(datasets[[i]]),c(vars(models[[i]]),xfix[[i]],keep))
      if (any(is.na(datasets[[i]][,myvars]))) {
        if (!warned) warning(paste0("Missing data encountered. Going for complete-case analysis"))
        warned  <- TRUE
        datasets[[i]] <- na.omit(datasets[[i]][,myvars,drop=FALSE])
      }
    }
  }

  exo <- exogenous(models)
  means <- lvms <- As <- Ps <- ps <- exs <- datas <- samplestat <- list()
  for (i in seq_len(nm)) {

    if (!is.null(exogenous(models[[i]]))) {
      if (any(is.na(exogenous(models[[i]])))) {
        exogenous(models[[i]]) <- exo
      }
    }

    mydata <- datasets[[i]]
    mymodel <- fixsome(models[[i]], data=mydata, measurement.fix=fix, exo.fix=exo.fix)
    mymodel <- updatelvm(mymodel,zeroones=TRUE,deriv=TRUE)

    P <- index(mymodel)$P1; P[P==0] <- NA
    P[!is.na(P) & !is.na(mymodel$covpar)] <- mymodel$covpar[!is.na(P) & !is.na(mymodel$covpar)]

    A <- index(mymodel)$M1; A[A==0] <- NA
    A[!is.na(A) & !is.na(mymodel$par)] <- mymodel$par[!is.na(A) & !is.na(mymodel$par)]

    mu <- unlist(mymodel$mean)[which(index(mymodel)$v1==1)]
    ex <- names(mymodel$expar)[which(index(mymodel)$e1==1)]

    p <- pars(mymodel, A, P, e=ex)
    p[p=="1"] <- NA

    means <- c(means, list(mu))
    lvms <- c(lvms, list(mymodel))
    datas <- c(datas, list(mydata))
    samplestat <- c(samplestat, list(procdata.lvm(models[[i]],data=mydata)))
    As <- c(As, list(A))
    Ps <- c(Ps, list(P))
    ps <- c(ps, list(p))
    exs <- c(exs, list(ex))
  };

######
  pp <- unlist(ps)
  parname <- unique(pp[!is.na(pp)])
  opt <- options(warn=-1); pidx <- is.na(as.numeric(parname)); options(opt)
  parname <- unique(unlist(pp[!is.na(pp)]));
  nfree <- sum(is.na(pp)) + length(parname)

  if (nfree>0) {
    pp0 <- lapply(ps, is.na)
    usedname <- cbind(parname, rep(NA,length(parname)))
    counter <- 1
    pres <- pres0 <- pp0
    for (i in seq_len(length(pp0))) {
      if (length(pp0[[i]]>0))
      for (j in seq_len(length(pp0[[i]]))) {
        pidx <- match(ps[[i]][j],parname)
        if (pp0[[i]][j]) {
          pres[[i]][j] <- paste0("p",counter)
          pres0[[i]][j] <- counter
          counter <- counter+1
        } else if (!is.na(pidx)) {
          if (!is.na(usedname[pidx,2])) {
            pres[[i]][j] <- usedname[pidx,2]
            pres0[[i]][j] <- as.numeric(substr(pres[[i]][j],2,nchar(pres[[i]][j])))
          } else {
            val <- paste0("p",counter)
            pres[[i]][j] <- val
            pres0[[i]][j] <- counter
            usedname[pidx,2] <- val
            counter <- counter+1
          }
        } else {
          pres[[i]][j] <- NA
        }
      }
    }
    mypar <- paste0("p",seq_len(nfree))
    myparPos <- pres0
    myparpos <- pres
    myparlist <- lapply(pres, function(x) x[!is.na(x)])
  } else {
    myparPos <- NULL
    mypar <- NULL
    myparpos <- NULL
    myparlist <- NULL
  }

  ### Mean parameter

  mm <- unlist(means)
  meanparname <- unique(mm[!is.na(mm)])
  opt <- options(warn=-1); midx <- is.na(as.numeric(meanparname)); options(opt)
  meanparname <- meanparname[midx]
  any.mean <- sum(is.na(mm)) + length(meanparname)
  nfree.mean <- sum(is.na(mm)) + length(setdiff(meanparname,parname))
  ## mean.fixed <- na.omit(match(parname,mm))
  mean.omit <- lapply(means,function(x) na.omit(match(parname,x)))
  nmean <- lapply(means,length)

  if (any.mean>0) {
    mm0 <- lapply(means, is.na)
    usedname <- cbind(meanparname, rep(NA,length(meanparname)))
    counter <- 1
    res0 <- res <- mm0
    for (i in seq_len(length(mm0))) {
      if (length(mm0[[i]])>0)
      for (j in seq_len(length(mm0[[i]]))) {
        midx <- match(means[[i]][j],meanparname)
        if (mm0[[i]][j]) {
          res[[i]][j] <- paste0("m",counter)
          res0[[i]][j] <- counter
          counter <- counter+1
        } else if (!is.na(midx)) {
          pidx <- match(meanparname[midx],pp)
          if (!is.na(pidx)) {
            res[[i]][j] <- unlist(myparlist)[pidx]
            res0[[i]][j] <- as.numeric(substr(res[[i]][j],2,nchar(res[[i]][j]))) +
              nfree.mean
              ##nmean[[i]]
          } else {
            if (!is.na(usedname[midx,2])) {
              res[[i]][j] <- usedname[midx,2]
              res0[[i]][j] <- as.numeric(substr(res[[i]][j],2,nchar(res[[i]][j])))
            } else {
              val <- paste0("m",counter)
              res[[i]][j] <- val
              res0[[i]][j] <- counter
              usedname[midx,2] <- val
              counter <- counter+1
            }
          }
        } else {
          res[[i]][j] <- NA
        }
      }
    }
    mymeanPos <- res0
    mymeanpos <- res
    mymeanlist <- lapply(res, function(x) x[!is.na(x)])
    mymean <- unique(unlist(mymeanlist))
  } else {
    mymeanPos <- NULL
    mymean <- NULL
    mymeanpos <- NULL
    mymeanlist <- NULL
  }

### Extra parameters

  N <- nfree+nfree.mean
  m0 <- p0 <- c()
  coefs <- coefsm <- mm0 <- mm <- pp0 <- pp <- c()
  for (i in seq_len(length(myparPos))) {
    mi <- mymeanPos[[i]]
    nmi <- length(mi)
    pi <- myparPos[[i]]
    p1 <- setdiff(pi,p0)
    p0 <- c(p0,p1)
    ##    pp0 <- c(pp0,list(match(p1,pi)+nfree.mean))
    pp0 <- c(pp0,list(match(p1,pi)))
    if (length(mean.omit[[i]])>0) mi <- mi[-mean.omit[[i]]]
    m1 <- setdiff(mi,m0)
    m0 <- c(m0,m1)
    mm0 <- c(mm0,list(match(m1,mi)))
    pp <- c(pp,list(c(m1,p1+nfree.mean)))
    if (length(p1)>0)
      coefs <- c(coefs,paste(i,coef(lvms[[i]],fix=FALSE,mean=FALSE)[pp0[[i]]],sep="@"))
    if (length(m1)>0) {
      coefsm0 <- paste(i,coef(lvms[[i]],fix=FALSE,mean=TRUE)[mm0[[i]]],sep="@")
      coefsm <- c(coefsm,coefsm0)
    }
  }
  coefs <- c(coefsm,coefs)

  res <- list(npar=nfree, npar.mean=nfree.mean,
              ngroup=length(lvms), names=mynames,
              lvm=lvms, data=datas, samplestat=samplestat,
              A=As, P=Ps, expar=exs,
              meanpar=names(mu), name=coefs, coef=pp, coef.idx=pp0,
              par=mypar, parlist=myparlist,  parpos=myparpos,
              mean=mymean, meanlist=mymeanlist, meanpos=mymeanpos,
              parposN=myparPos,
              meanposN=mymeanPos,
              models.orig=models.orig, missing=missing
              )
  class(res) <- "multigroup"
  checkmultigroup(res)
  return(res)
}

###}}}

###{{{ checkmultigroup

checkmultigroup <- function(x) {
    ## Check validity:
  for (i in seq_len(x$ngroup)) {
    if (nrow(x$data[[i]])<2) {
      warning("With only one observation in the group, all parameters should be inherited from another a group!")
    }
  }
}

###}}} checkmultigroup
