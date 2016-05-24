MVA.load <- function(x,xax=1,yax=2,set=c(12,1,2),space=1,...) {
  if (length(set)!=1) {set <- 12}
  if (!set %in% c(12,1,2)) {stop("wrong 'set'")}
  x <- MVA.ident(x)
  loads.temp <- if (inherits(x,c("PCA.ade4","PCA.prcomp","PCA.princomp","PCA.mixOmics","PCA.labdsv","PCA.vegan",
    "sPCA.mixOmics","IPCA.mixOmics","sIPCA.mixOmics","LDA.MASS","LDA.ade4","PLSDA.mixOmics","sPLSDA.mixOmics",
    "Multilevel.sPLSDA.mixOmics","CPPLS.pls","PLSR.pls","PLSR.mixOmics","PLSR.plsRglm","sPLSR.mixOmics",
    "Multilevel.sPLSR.mixOmics","PLSGLR.plsRglm","PCR.pls","CDA.ade4","NSCOA.ade4","MCA.ade4","Mix.ade4",
    "RDA.ade4","RDAortho.ade4"))) {MVA.get.loads(x)} else
    if (inherits(x,c("PCIA.ade4"))) {MVA.get.loads(x,set)} else
    if (inherits(x,c("RDA.vegan","CIA.ade4","rCCorA.mixOmics","2BPLS.mixOmics","2BsPLS.mixOmics",
    "Multilevel.2BsPLS.mixOmics","rGCCA.RGCCA","rGCCA.mixOmics","sGCCA.RGCCA","sGCCA.mixOmics")))
    {MVA.get.loads(x,space)} else
    {MVA.get.loads(x)}
  loads <- if (is.data.frame(loads.temp)) {loads.temp} else {loads.temp[[1]]}
  keep <- apply(abs(as.data.frame(loads[,c(xax,yax)])),1,sum)>0
  loads <- loads[keep,]
  if (!is.data.frame(loads.temp)) {loads.temp[[2]] <- loads.temp[[2]][keep]}
  if (!xax %in% c(1:ncol(loads))) {stop("wrong 'xax'")}
  if (ncol(loads)==1) {
    xax <- 1
    yax <- NULL
  }
  if (!is.null(yax) && !yax %in% c(1:ncol(loads))) {
    warning("wrong 'yax', only 'xax' used")
    yax <- NULL
  }
  loadsx <- loads[,xax]
  loadsy <- NULL
  if (!is.null(yax) && ncol(loads)>1) {loadsy <- loads[,yax]}
  res.temp <- as.data.frame(cbind(loadsx,loadsy))
  rownames(res.temp) <- rownames(loads)
  if (inherits(x,c("RDA.vegan"))) {
    colnames(res.temp)[1] <- paste(ifelse(space==1,"Constr. comp.","Unconstr. comp."),xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste(ifelse(space==1,"Constr. comp.","Unconstr. comp."),yax)}
  } else if (inherits(x,c("RDA.ade4"))) {
    colnames(res.temp)[1] <- paste("Constr. comp.",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Constr. comp.",yax)}
  } else if (inherits(x,"RDAortho.ade4")) {
    colnames(res.temp)[1] <- paste("Unconstr. comp.",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Unconstr. comp.",yax)}
  } else if (inherits(x,c("rCCorA.mixOmics"))) {
    colnames(res.temp)[1] <- paste("Canonical axis",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Canonical axis",yax)}
  } else if (inherits(x,c("CIA.ade4"))) {
    colnames(res.temp)[1] <- paste("Coinertia axis",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Coinertia axis",yax)}
  } else {
    colnames(res.temp)[1] <- paste("Comp.",xax)
    if (ncol(res.temp)==2) {colnames(res.temp)[2] <- paste("Comp.",yax)}
  }
  res <- list(loads=res.temp)
  if (!is.data.frame(loads.temp)) {
    res[[2]] <- loads.temp[[2]]
    names(res)[2] <- names(loads.temp)[2]
  }
  return(res)
}

MVA.get.loads <- function(x,...) {
  UseMethod("MVA.get.loads")
}

MVA.get.loads.default <- MVA.get.loads.unknown <- function(x,...) {
  stop("unknown multivariate analysis or no loading available")
}

MVA.get.loads.PCA.ade4 <- MVA.get.loads.NSCOA.ade4 <-
MVA.get.loads.MCA.ade4 <- MVA.get.loads.RDA.ade4 <-
MVA.get.loads.RDAortho.ade4 <- function(x,...) {as.data.frame(x$c1)}

MVA.get.loads.PCA.prcomp <- MVA.get.loads.PCA.mixOmics <-
MVA.get.loads.sPCA.mixOmics <- function(x,...) {as.data.frame(x$rotation)}

MVA.get.loads.PCA.princomp <- MVA.get.loads.PCA.labdsv <- 
MVA.get.loads.IPCA.mixOmics <- MVA.get.loads.sIPCA.mixOmics <- 
MVA.get.loads.CPPLS.pls <- MVA.get.loads.PLSR.pls <-
MVA.get.loads.PCR.pls <- function(x,...) {as.data.frame(unclass(x$loadings))}

MVA.get.loads.PCA.vegan <- function(x,xax,yax,...) {as.data.frame(x$CA$v)}

MVA.get.loads.LDA.MASS <- function(x,...) {as.data.frame(x$scaling)}

MVA.get.loads.LDA.ade4 <- MVA.get.loads.CDA.ade4 <- function(x,...) {as.data.frame(x$fa)}

MVA.get.loads.PLSDA.mixOmics <- MVA.get.loads.sPLSDA.mixOmics <-
MVA.get.loads.Multilevel.sPLSDA.mixOmics <- MVA.get.loads.PLSR.mixOmics <-
MVA.get.loads.sPLSR.mixOmics <- MVA.get.loads.Multilevel.sPLSR.mixOmics <- function(x,...) {as.data.frame(x$loadings$X)}

MVA.get.loads.PLSR.plsRglm <- MVA.get.loads.PLSGLR.plsRglm <- function(x,...) {as.data.frame(x$pp)}

MVA.get.loads.Mix.ade4 <- function(x,...) {
  tab <- as.data.frame(x$c1)
  index <- x$index
  if (any("o" %in% index)) {
    index[which(index=="o")] <- "f"
    index <- droplevels(index)
  }
  names(index) <- unique(unlist(lapply(strsplit(rownames(tab),split="[.]"),function(y) y[1])))
  res <- list(tab,index=index)
  return(res)
}

MVA.get.loads.PCIA.ade4 <- function(x,set,...) {
  if (set==1) {
    tab <- as.data.frame(x$loadX)
  } else if (set==2) {
    tab <- as.data.frame(x$loadY)
  } else {
    X <- x$loadX
    Y <- x$loadY
    colnames(Y) <- colnames(X)
    if (any(rownames(X) %in% rownames(Y))) {
	rownames(X) <- paste0("X.",rownames(X))
	rownames(Y) <- paste0("Y.",rownames(Y))
    }
    tab <- as.data.frame(rbind(X,Y))
  }
  if (set==12) {
    res <- list(coord=tab)
    res$set <- gl(2,nrow(x$loadX),labels=c("X","Y"))
  } else {res <- tab}
  return(res)
}

MVA.get.loads.RDA.vegan <- function(x,space,...) {
  if (!space %in% c(1,2)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  tab <- if (space==1) {
    as.data.frame(x$CCA$v)
  } else {
    as.data.frame(x$CA$v)
  }
  return(tab)
}

MVA.get.loads.CIA.ade4 <- function(x,space,...) {
  if (!space %in% c(1,2)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  tab <- if (space==1) {
    as.data.frame(x$c1)
  } else {
    as.data.frame(x$l1)
  }
  return(tab)
}

MVA.get.loads.rCCorA.mixOmics <- MVA.get.loads.2BPLS.mixOmics <-
MVA.get.loads.2BsPLS.mixOmics <- MVA.get.loads.Multilevel.2BsPLS.mixOmics <- function(x,space,...) {
  if (!space %in% c(1,2)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  tab <- if (space==1) {
    as.data.frame(x$loadings$X)
  } else {
    as.data.frame(x$loadings$Y)
  }
  return(tab)
}

MVA.get.loads.rGCCA.RGCCA <- MVA.get.loads.sGCCA.RGCCA <- function(x,space,...) {
  if (!space %in% 1:length(x$a)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  res <- as.data.frame(x$a[[space]])
  return(res)
}

MVA.get.loads.rGCCA.mixOmics <- MVA.get.loads.sGCCA.mixOmics <- function(x,space,...) {
  if (!space %in% 1:length(x$variates)) {stop("wrong 'space'")}
  if (length(space)!=1) {space <- 1}
  res <- as.data.frame(x$loadings[[space]])
  return(res)
}

