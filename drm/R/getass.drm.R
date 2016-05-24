"getass.drm" <-
  function (dep, assoc, parms, nrep, npar, nclass, data, subset,cluster) 
{
  asslength <- sapply(assoc, function(i) length(as.character(i)))
  assmodel <- assoc[asslength > 2]
  tfunction <- asspos <- alab <- parm <- equal <- constant <- NULL
  if (any(asslength > 2)) {
    parm <- sapply(assmodel, paste)[2, ]
    tpos <- grep("tau", parm)
    if (length(grep("M", dep)) > 0 & length(tpos)>0) {
      pnam <- names(parms)
      tfunc <- as.list(sapply(assmodel[tpos], paste)[3, ])
      names(tfunc) <- parm[tpos]
      tfn <- lapply(tfunc, function(i) {
        if (exists(i)) 
          get(i)
        else{
          if(!is.na(charmatch("function", i))){
            eval(parse(text = i))}
          else(NULL) }})
      tpos <- tpos[!sapply(tfn, is.null)]
      tfn <- tfn[!sapply(tfn, is.null)]
      if (length(tpos) > 0) {
        if (!all(sapply(tfn, is.function))) 
          stop("Argument passed to tau** has to be a function")
        tlen <- lapply(tfn, function(tau.fcn) tau.fcn())
        if (length(grep("M2", dep)) > 0) {
          tfunction <- rep(list(function(x) x), 3)
        }
        else {
          if (any(sapply(tlen, length) != (nrep - 1))) 
            stop(paste("Number of adjacent overlapping bivariate distributions and\n", 
                       "length of the resulting object passed to tau in 'dep' do not match"))
          tfunction <- rep(list(function(x) x), (nclass -  1)^2)
        }
        names(tfunction) <- names(parms)[grep("tau", names(parms))]
        tfunction[match(names(tfn), names(tfunction))] <- tfn
        tdef <- lapply(tfn, formals)
        tparm <- match(names(tdef), pnam)
        parms <- as.list(parms)
        names(parms) <- pnam
        parms[tparm] <- tdef
        parm <- parm[-tpos]
        assmodel <- assmodel[-tpos]
      }
    }
    if (length(parm) > 0) {
      lp <- length(parm)
      lps <- length(unlist(parms))
      data <- cbind(data,
                    matrix(rep(rep(1, nrow(data)), lp), ncol = lp))
      names(data)[(ncol(data) - lp + 1):ncol(data)] <- parm
      assmod <- lapply(assmodel, function(i, data, subset, nrep, 
                                          cluster) {
        m <- model.matrix(terms(i, data = data),
                          model.frame(i,data = data),
                          contrasts.arg = c("contr.treatment","contr.treatment"))
        if(length(attributes(m)$contrasts)==0)
          stop("Only factors allowed as explanatory variables in the association model")
        dimnames(m)[[2]][1] <- paste(i)[2]
        dimnames(m)[[1]] <- cluster
        m <- m[order(cluster), -2, drop = FALSE]
        m
      }, data = data, nrep = nrep, cluster = cluster, subset=subset)
      assmodel <- do.call("cbind", assmod)
      asslab <- dimnames(assmodel)[[2]]
      parms <- c(parms, rep(0, length(asslab)))
      names(parms)[(length(parms) - length(asslab) + 1):length(parms)] <- asslab
      parms <- parms[!duplicated(names(parms))]
      Alab <- apply(matrix(apply(assmodel, 1, paste, collapse = ""), 
                           ncol = nrep, byrow = TRUE), 1, paste, collapse = "")
      alab <- match(dimnames(assmodel)[[1]], grep("FALSE", 
                                                  paste(duplicated(Alab))))
      assmodel <- assmodel[!is.na(alab), drop = FALSE, ]
      tba <- tapply(dimnames(assmodel)[[1]], apply(assmodel, 
                                                   1, paste, collapse = ""), table)
      if (any(unlist(tba) != nrep)) 
        stop("Time-varying covariate can't be used in the association model")
      assmodel <- assmodel[!duplicated(dimnames(assmodel)[[1]]), 
                           drop = FALSE, ]
      alab <- c(1:length(unique(Alab)))[match(Alab, unique(Alab))]
      assmodel <- lapply(parm,
                         function(i, assmodel, asslab)
                         assmodel[, grep(i, asslab), drop = FALSE], assmodel = assmodel, 
                         asslab = asslab)
      ##asspos <- match(parm, names(parms))
    }
  }
  assmatch <- lapply(parm, function(i, parms, npar) {
    grep(i, names(unlist(parms))) - npar
  }, parms = parms, npar = npar)
  asspos <- match(parm, names(unlist(parms))) - npar
  if (any(asslength == 2)) {
    asseq <- lapply(assoc[asslength < 3], function(i) lapply(i, 
                                                             paste))
    assoc <- lapply(asseq, unlist)
    if (length(grep("==", do.call("rbind", assoc)[, 2])) < 
        length(asseq)) 
      stop("Can only handle \"==\" : Use e.g. dep=list(\"M\",~tau23==tau32)")
    assoc <- lapply(assoc, function(i) i[-c(1, 2)])
    pos <- lapply(assoc, function(i, parms) match(i, names(parms)), 
                  parms = parms)
    num <- lapply(seq(pos), function(i, assoc, pos) assoc[[i]][is.na(pos[[i]])], 
                  assoc = assoc, pos = pos)
    assnum <- assoc[sapply(num, length) > 0]
    assoc <- assoc[sapply(num, length) == 0]
    posnum <- sapply(pos[sapply(num, length) > 0], function(i) i[1])
    pos <- pos[sapply(num, length) == 0]
    num <- sapply(assnum, function(i) eval(parse(text = i[2])))
    if (any(is.na(num))) 
      stop(paste("Couldn't find parameter", unlist(assoc)[is.na(unlist(pos))]))
    constant <- list(posnum, num)
    if (length(unlist(pos)) > 0) {
      for (i in 1:length(pos)) {
        names(parms)[pos[[i]][2]] <- names(parms)[pos[[i]][1]]
        parms[pos[[i]][2]] <- parms[pos[[i]][1]]
        if (!is.null(tfunction)) 
          tfunction[assoc[[i]][2]] <- tfunction[assoc[[i]][1]]
      }
    }
    pos <- lapply(assnum, function(i, parms) match(i, names(parms)), 
                  parms = unlist(parms))
    posnum <- sapply(pos[sapply(num, length) > 0], function(i) i[1])
    constant <- list(posnum, num)
    if (length(constant[[1]]) > 0) {
      parms <- unlist(parms)[-constant[[1]]]
      constant[[1]] <- (constant[[1]] - npar)[(constant[[1]] - 
                                               npar) > 0]
    }
    else (parms <- unlist(parms))
    equal <- names(parms)
    parms <- parms[!duplicated(equal)]
    equal <- match(equal, names(parms))
    equal <- (equal - npar)[(equal - npar) > 0]
  }
  list(parms = unlist(parms), equal = equal,
       ass = list(assmodel = assmodel, asspos = asspos,
         assmatch = assmatch, alab = alab), tfunction = tfunction, 
       constant = constant)
}

