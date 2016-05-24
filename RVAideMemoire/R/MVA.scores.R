MVA.scores <- function(x,xax=1,yax=2,scaling=2,set=c(12,1,2),space=1,...) {
  if (length(set)!=1) {set <- 12}
  if (!set %in% c(12,1,2)) {stop("wrong 'set'")}
  x <- MVA.ident(x)
  coord.temp <- if (inherits(x,c("PCA.ade4","PCA.prcomp","PCA.princomp","PCA.mixOmics",
    "PCA.labdsv","sPCA.mixOmics","IPCA.mixOmics","sIPCA.mixOmics","PCoA.ade4","PCoA.ape",
    "PCoA.labdsv","PCoA.stats","PCoA.vegan","nMDS.MASS","nMDS.mono.vegan","nMDS.iso.vegan",
    "LDA.MASS","LDA.ade4","PLSDA.mixOmics","sPLSDA.mixOmics","CPPLS.pls","PLSR.pls",
    "PLSR.mixOmics","PLSR.plsRglm","sPLSR.mixOmics","PLSGLR.plsRglm","PCR.pls","CDA.ade4",
    "NSCOA.ade4","MCA.ade4","Mix.ade4","RDA.ade4","RDAortho.ade4",
    "Multilevel.sPLSDA.mixOmics","Multilevel.sPLSR.mixOmics"))) {MVA.get.scores(x)} else
    if (inherits(x,c("PCA.vegan"))) {MVA.get.scores(x,xax,yax,scaling)} else
    if (inherits(x,c("COA.ade4","DCOA.ade4","PCIA.ade4","CCA.ade4","DPCoA.ade4"))) {MVA.get.scores(x,set)} else
    if (inherits(x,c("COA.vegan"))) {MVA.get.scores(x,xax,yax,scaling,set)} else
    if (inherits(x,c("RDA.vegan","dbRDA.vegan"))) {MVA.get.scores(x,xax,yax,scaling,space)} else
    if (inherits(x,c("CCA.vegan"))) {MVA.get.scores(x,xax,yax,scaling,set,space)} else
    if (inherits(x,c("CCorA.vegan","rCCorA.mixOmics","2BPLS.mixOmics","2BsPLS.mixOmics",
    "Multilevel.2BsPLS.mixOmics","rGCCA.RGCCA","rGCCA.mixOmics","sGCCA.RGCCA",
    "sGCCA.mixOmics"))) {MVA.get.scores(x,space)} else
    if (inherits(x,c("CIA.ade4"))) {MVA.get.scores(x,set,space)} else
    {MVA.get.scores(x)}
  coord <- if (is.data.frame(coord.temp)) {coord.temp} else {coord.temp[[1]]}
  if (!inherits(x,c("PCA.vegan","COA.vegan","RDA.vegan","dbRDA.vegan","CCA.vegan"))) {
    if (!xax %in% c(1:ncol(coord))) {stop("wrong 'xax'")}
    if (ncol(coord)==1) {
	xax <- 1
	yax <- NULL
    }
    if (!is.null(yax) && !yax %in% c(1:ncol(coord))) {
	warning("wrong 'yax', only 'xax' used")
	yax <- NULL
    }
    coordx <- coord[,xax]
    coordy <- NULL
    if (!is.null(yax) && ncol(coord)>1) {coordy <- coord[,yax]}
  } else {
    coordx <- coord[,1]
    coordy <- NULL
    if (!is.null(yax) && ncol(coord)>1) {coordy <- coord[,2]}
  }
  res.temp <- as.data.frame(cbind(coordx,coordy))
  rownames(res.temp) <- rownames(coord)
  if (inherits(x,c("RDA.vegan","dbRDA.vegan","CCA.vegan"))) {
    colnames(res.temp)[1] <- paste(ifelse(space==1,"Constr. comp.","Unconstr. comp."),xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste(ifelse(space==1,"Constr. comp.","Unconstr. comp."),yax)}
  } else if (inherits(x,c("RDA.ade4","CCA.ade4"))) {
    colnames(res.temp)[1] <- paste("Constr. comp.",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Constr. comp.",yax)}
  } else if (inherits(x,"RDAortho.ade4")) {
    colnames(res.temp)[1] <- paste("Unconstr. comp.",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Unconstr. comp.",yax)}
  } else if (inherits(x,c("CCorA.vegan","rCCorA.mixOmics"))) {
    colnames(res.temp)[1] <- paste("Canonical axis",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Canonical axis",yax)}
  } else if (inherits(x,c("CIA.ade4"))) {
    colnames(res.temp)[1] <- paste("Coinertia axis",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Coinertia axis",yax)}
  } else {
    colnames(res.temp)[1] <- paste("Comp.",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Comp.",yax)}
  }
  res <- list(coord=res.temp)
  if (!is.data.frame(coord.temp)) {
    res[[2]] <- coord.temp[[2]]
    names(res)[2] <- names(coord.temp)[2]
  }
  return(res)
}

MVA.get.scores <- function(x,...) {
  UseMethod("MVA.get.scores")
}

MVA.get.scores.default <- MVA.get.scores.unknown <- function(x,...) {
  stop("unknown multivariate analysis")
}

MVA.get.scores.PCA.ade4 <- MVA.get.scores.PCoA.ade4 <-
MVA.get.scores.LDA.ade4 <- MVA.get.scores.CDA.ade4 <-
MVA.get.scores.NSCOA.ade4 <- MVA.get.scores.MCA.ade4 <-
MVA.get.scores.Mix.ade4 <- MVA.get.scores.RDA.ade4 <-
MVA.get.scores.RDAortho.ade4 <- function(x,...) {as.data.frame(x$li)}

MVA.get.scores.PCA.prcomp <- function(x,...) {
  if ("x" %in% names(x)) {res <- as.data.frame(x$x)} else
    {stop("no scores available, compute the analysis with 'retx=TRUE'")}
}

MVA.get.scores.PCA.princomp <- function(x,...) {
  if (!is.null(x$scores)) {res <- as.data.frame(x$scores)} else
    {stop("no scores available, compute the analysis with 'scores=TRUE'")}
}

MVA.get.scores.PCA.mixOmics <- MVA.get.scores.sPCA.mixOmics <-
MVA.get.scores.IPCA.mixOmics <- MVA.get.scores.sIPCA.mixOmics <- function(x,...) {as.data.frame(x$x)}

MVA.get.scores.PCA.labdsv <- function(x,...) {as.data.frame(x$scores)}

MVA.get.scores.PCA.vegan <- function(x,xax,yax,scaling,...) {
  sumev <- x$tot.chi
  slam <- sqrt(x$CA$eig[c(xax,yax)]/sumev)
  nr <- nrow(x$CA$u)
  const <- sqrt(sqrt((nr-1)*sumev))
  if (length(const)==1) {const <- c(const,const)}
  wa <- x$CA$u[,c(xax,yax),drop=FALSE]
  scal <- list(slam,1,sqrt(slam))[[abs(scaling)]]
  wa <- sweep(wa,2,scal,"*")
  wa <- const[2]*wa
  return(as.data.frame(wa))
}

MVA.get.scores.PCoA.ape <- function(x,...) {
  if (!is.null(x$vectors.cor)) {res <- as.data.frame(x$vectors.cor)} else
    {res <- as.data.frame(x$vectors)}
  return(res)
}

MVA.get.scores.PCoA.labdsv <- MVA.get.scores.PCoA.stats <-
MVA.get.scores.nMDS.MASS <- MVA.get.scores.nMDS.mono.vegan <-
MVA.get.scores.nMDS.iso.vegan <- MVA.get.scores.nMDS.labdsv <- function(x,...) {as.data.frame(x$points)}

MVA.get.scores.PCoA.vegan <- function(x,...) {
  if (inherits(x,"wcmdscale")) {res <- as.data.frame(x$points)} else
    {res <- as.data.frame(x$CA$u)}
  return(res)
}

MVA.get.scores.LDA.MASS <- function(x,...) {as.data.frame(LDA.format(x)$li)}

MVA.get.scores.PLSDA.mixOmics <- MVA.get.scores.sPLSDA.mixOmics <-
MVA.get.scores.PLSR.mixOmics <- MVA.get.scores.sPLSR.mixOmics <-
MVA.get.scores.Multilevel.sPLSDA.mixOmics <- MVA.get.scores.Multilevel.sPLSR.mixOmics <- 
function(x,...) {as.data.frame(x$variates$X)}

MVA.get.scores.CPPLS.pls <- MVA.get.scores.PLSR.pls <- 
MVA.get.scores.PCR.pls <- function(x,...) {as.data.frame(cbind(x$scores))}

MVA.get.scores.PLSR.plsRglm <- MVA.get.scores.PLSGLR.plsRglm <- function(x,...) {as.data.frame(x$tt)}

MVA.get.scores.COA.ade4 <- MVA.get.scores.DCOA.ade4 <-
MVA.get.scores.CCA.ade4 <- function(x,set,...) {
  if (set==1) {
    tab <- as.data.frame(x$li)
  } else if (set==2) {
    tab <- as.data.frame(x$co)
  } else {
    X <- x$li
    Y <- x$co
    colnames(Y) <- colnames(X)
    tab <- as.data.frame(rbind(X,Y))
  }
  if (set==12) {
    res <- list(coord=tab)
    res$set <- factor(rep(c("rows","columns"),c(nrow(x$li),nrow(x$co))))
  } else {res <- tab}
  return(res)
}

MVA.get.scores.DPCoA.ade4 <- function(x,set,...) {
  if (set==1) {
    tab <- as.data.frame(x$dls)
  } else if (set==2) {
    tab <- as.data.frame(x$li)
  } else {
    X <- x$dls
    Y <- x$li
    colnames(Y) <- colnames(X)
    tab <- as.data.frame(rbind(X,Y))
  }
  if (set==12) {
    res <- list(coord=tab)
    res$set <- factor(rep(c("categories","collections"),c(nrow(x$dls),nrow(x$li))))
  } else {res <- tab}
  return(res)
}

MVA.get.scores.COA.vegan <- function(x,xax,yax,scaling,set,...) {
  slam <- sqrt(x$CA$eig[c(xax,yax)])
  wa <- x$CA$u[,c(xax,yax),drop=FALSE]
  scal <- list(slam,1,sqrt(slam))[[abs(scaling)]]
  wa <- sweep(wa,2,scal,"*")
  if (scaling<0) {
   scal <- sqrt(1/(1-slam^2))
    wa <- sweep(wa,2,scal,"*")
  }
  ro <- as.data.frame(wa)
  v <- x$CA$v[,c(xax,yax),drop=FALSE]
  scal <- list(1,slam,sqrt(slam))[[abs(scaling)]]
  v <- sweep(v,2,scal,"*")
  if (scaling<0) {
    scal <- sqrt(1/(1-slam^2))
    v <- sweep(v,2,scal,"*")
  }
  co <- v
  if (set==1) {
    tab <- as.data.frame(ro)
  } else if (set==2) {
    tab <- as.data.frame(co)
  } else {
    tab <- as.data.frame(rbind(ro,co))
  }
  if (set==12) {
    res <- list(coord=tab)
    res$set <- factor(rep(c("rows","columns"),c(nrow(ro),nrow(co))))
  } else {res <- tab}
  return(res)
}

MVA.get.scores.PCIA.ade4 <- function(x,set,...) {
  tab <- if (set==1) {
    as.data.frame(x$scorX)
  } else if (set==2) {
    as.data.frame(x$scorY)
  } else {
    as.data.frame(rbind(x$scorX,x$scorY))
  }
  if (set==12) {
    res <- list(coord=tab)
    res$set <- gl(2,nrow(x$scorX),labels=c("X","Y"))
  } else {res <- tab}
  return(res)
}

MVA.get.scores.RDA.vegan <- MVA.get.scores.dbRDA.vegan <- function(x,xax,yax,scaling,space,...) {
  if (!space%in%c(1,2)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  tab <- if (space==1) {x$CCA$eig} else {x$CA$eig}
  if (!xax %in% 1:length(tab)) {stop("wrong 'xax'")}
  if (length(tab)==1) {
    xax <- 1
    yax <- NULL
  }
  if (!is.null(yax) && !yax%in%c(1:length(tab))) {
    warning("wrong 'yax', only 'xax' used")
    yax <- NULL
  }
  sumev <- x$tot.chi
  slam <- sqrt(tab[c(xax,yax)]/sumev)
  nr <- nrow(x$CCA$u)
  const <- sqrt(sqrt((nr-1)*sumev))
  if (length(const)==1) {const <- c(const,const)}
  wa <- if (space==1) {
    x$CCA$wa[,c(xax,yax),drop=FALSE]
  } else {
    x$CA$u[,c(xax,yax),drop=FALSE]
  }
  scal <- list(slam,1,sqrt(slam))[[abs(scaling)]]
  wa <- sweep(wa,2,scal,"*")
  wa <- const[2]*wa
  return(as.data.frame(wa))
}

MVA.get.scores.CCA.vegan <- function(x,xax,yax,scaling,set,space,...) {
  if (!space%in%c(1,2)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  tab <- if (space==1) {x$CCA$eig} else {x$CA$eig}
  if (!xax %in% 1:length(tab)) {stop("wrong 'xax'")}
  if (length(tab)==1) {
    xax <- 1
    yax <- NULL
  }
  slam <- sqrt(tab[c(xax,yax)])
  wa <- if (space==1) {
    x$CCA$wa[,c(xax,yax),drop=FALSE]
  } else {
    x$CA$u[,c(xax,yax),drop=FALSE]
  }
  scal <- list(slam,1,sqrt(slam))[[abs(scaling)]]
  wa <- sweep(wa,2,scal,"*")
  if (scaling<0) {
   scal <- sqrt(1/(1-slam^2))
    wa <- sweep(wa,2,scal,"*")
  }
  ro <- wa
  v <- if (space==1) {
    x$CCA$v[,c(xax,yax),drop=FALSE]
  } else {
    x$CA$v[,c(xax,yax),drop=FALSE]
  }
  scal <- list(1,slam,sqrt(slam))[[abs(scaling)]]
  v <- sweep(v,2,scal,"*")
  if (scaling<0) {
    scal <- sqrt(1/(1-slam^2))
    v <- sweep(v,2,scal,"*")
  }
  co <- v
  if (set==1) {
    tab <- as.data.frame(ro)
  } else if (set==2) {
    tab <- as.data.frame(co)
  } else {
    tab <- as.data.frame(rbind(ro,co))
  }
  if (set==12) {
    res <- list(coord=tab)
    res$set <- factor(rep(c("rows","columns"),c(nrow(ro),nrow(co))))
  } else {res <- tab}
  return(res)
}

MVA.get.scores.CCorA.vegan <- function(x,space,...) {
  if (!space%in%c(1,2)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  res <- if (space==1) {
    as.data.frame(x$Cx)
  } else {
    as.data.frame(x$Cy)
  }
  return(res)
}

MVA.get.scores.rCCorA.mixOmics <- MVA.get.scores.2BPLS.mixOmics <- 
MVA.get.scores.2BsPLS.mixOmics <- MVA.get.scores.Multilevel.2BsPLS.mixOmics <- function(x,space,...) {
  if (!space %in% c(1,2,3)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 3}
  res <- if (space==1) {
    as.data.frame(x$variates$X)
  } else if (space==2) {
    as.data.frame(x$variates$Y)
  } else {
    as.data.frame((x$variates$X+x$variates$Y)/2)
  }
  return(res)
}

MVA.get.scores.CIA.ade4 <- function(x,set,space,...) {
  if (!space %in% c(1,2,3)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  if (space==1) {
    tab <- as.data.frame(x$lX)
  } else if (space==2) {
    tab <- as.data.frame(x$lY)
  } else {
    if (set==1) {
	tab <- as.data.frame(x$mX)
    } else if (set==2) {
	tab <- as.data.frame(x$mY)
    } else {
	X <- x$mX
	Y <- x$mY
	colnames(Y) <- colnames(X)
	tab <- as.data.frame(rbind(X,Y))
	rownames(tab) <- c(paste0("X.",rownames(x$mX)),paste0("Y.",rownames(x$mY)))
    }
  }
  if (space==3 && set==12) {
    res <- list(coord=tab)
    res$set <- factor(rep(c("X","Y"),c(nrow(x$mX),nrow(x$mY))))
  } else {res <- tab}
  return(res)
}

MVA.get.scores.rGCCA.RGCCA <- MVA.get.scores.sGCCA.RGCCA <- function(x,space,...) {
  if (!space %in% 1:length(x$Y)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  res <- as.data.frame(x$Y[[space]])
  return(res)
}

MVA.get.scores.rGCCA.mixOmics <- MVA.get.scores.sGCCA.mixOmics <- function(x,space,...) {
  if (!space %in% 1:length(x$variates)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  res <- as.data.frame(x$variates[[space]])
  return(res)
}

