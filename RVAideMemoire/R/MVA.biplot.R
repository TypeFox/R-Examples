MVA.biplot <- function(x,xax=1,yax=2,scaling=2,sco.set=c(12,1,2),
  cor.set=c(12,1,2),space=1,ratio=0.9,weights=1,constraints=c("nf","n","f",NULL),
  sco.args=list(),cor.args=list(),f.col=1,f.cex=1) {
  coord <- MVA.scores(x,xax,yax,scaling,sco.set,space)$coord
  if (ncol(coord)==1) {stop("biplots cannot be drawn with only one dimension")}
  x <- MVA.ident(x)
  arg.sco <- c(list(as.name("MVA.scoreplot"),x=x,xax=xax,yax=yax,
    scaling=scaling,space=space,set=sco.set,keepmar=TRUE),sco.args)
  call.sco <- as.call(arg.sco)
  lims <- eval(call.sco)
  constr <- inherits(x,c("CCA.vegan","CCA.ade4","RDA.vegan","RDA.ade4"))
  if (constr & !is.null(constraints)) {
    constr.type <- MVA.constr(x)
    if (constr.type=="n") {
	constraints <- "n"
    } else if (constr.type=="f") {
	constraints <- "f"
	cor.set <- 2
    }
  }
  if (!constr | (constr & is.null(constraints))) {
    corr <- suppressWarnings(MVA.cor(x,xax,yax,cor.set,space)$corr)
    left <- min(corr[,1])/lims$xlim[1]
    right <- max(corr[,1])/lims$xlim[2]
    bottom <- min(corr[,2])/lims$ylim[1]
    top <- max(corr[,2])/lims$ylim[2]
    ratios <- c(left,right,bottom,top)
    rmax <- max(ratios)
    const <- 1/rmax*ratio
    arg.cor <- c(list(as.name("MVA.corplot"),x=x,xax=xax,yax=yax,
	space=space,set=cor.set,add=TRUE,add.const=const,xlab="",
	ylab="",circle=FALSE,intcircle=NULL,drawintaxes=FALSE),cor.args)
    call.cor <- as.call(arg.cor)
    eval(call.cor)
  } else {
    if (constraints=="n") {
	corr <- suppressWarnings(MVA.cor(x,xax,yax,cor.set,space)$corr)
	left <- min(corr[,1])/lims$xlim[1]
	right <- max(corr[,1])/lims$xlim[2]
	bottom <- min(corr[,2])/lims$ylim[1]
	top <- max(corr[,2])/lims$ylim[2]
	ratios <- c(left,right,bottom,top)
	rmax <- max(ratios)
	const <- 1/rmax*ratio
	arg.cor <- c(list(as.name("MVA.corplot"),x=x,xax=xax,yax=yax,
	  space=space,set=cor.set,add=TRUE,add.const=const,xlab="",
	  ylab="",circle=FALSE,intcircle=NULL,drawintaxes=FALSE),cor.args)
	call.cor <- as.call(arg.cor)
	eval(call.cor)
    }
    if (constraints %in% c("nf","f")) {
	corr <- suppressWarnings(MVA.cor(x,xax,yax,cor.set,space)$corr)
	left <- min(corr[,1])/lims$xlim[1]
	right <- max(corr[,1])/lims$xlim[2]
	bottom <- min(corr[,2])/lims$ylim[1]
	top <- max(corr[,2])/lims$ylim[2]
	ratios <- c(left,right,bottom,top)
	rmax <- max(ratios)
	const <- 1/rmax*ratio
	arg.cor <- c(list(as.name("MVA.corplot"),x=x,xax=xax,yax=yax,
	  space=space,set=cor.set,add=TRUE,add.const=const,xlab="",
	  ylab="",circle=FALSE,intcircle=NULL,drawintaxes=FALSE),cor.args)
	call.cor <- as.call(arg.cor)
	eval(call.cor)
	if (space==1) {
	  if (length(weights)!=nrow(coord)) {
	    if (length(weights)==1) {
		weights=rep(weights,nrow(coord))
	    } else {stop("non-convenient 'weights' argument")}
	  }
	  centr <- MVA.centr(x)
	  if (length(f.col)!=length(centr$lev)) {
	    if (length(f.col)==1) {f.col <- rep(f.col,length(centr$lev))} else
	    if (length(f.col)==ncol(centr$fac)) {
		f.col <- rep(f.col,apply(centr$fac,2,function(z) nlevels(factor(z))))
	    } else {stop("non convenient 'f.col' argument")}
	  }
	  if (length(f.cex)!=length(centr$lev)) {
	    if (length(f.cex)==1) {f.cex <- rep(f.cex,length(centr$lev))} else
	    if (length(f.cex)==ncol(centr$fac)) {
		f.cex <- rep(f.cex,apply(centr$fac,2,function(z) nlevels(factor(z))))
	    } else {stop("non convenient 'f.cex' argument")}
	  }
	  bar.x <- bar.y <- list()
	  length(bar.x) <- length(bar.y) <- ncol(centr$fac)
	  for (i in 1:length(bar.x)) {
	    f <- centr$fac[,i]
	    nlev <- nlevels(f)
	    bar.x[[i]] <- bar.y[[i]] <- integer(nlev)
	    for (j in 1:nlev) {
		bar.x[[i]][j] <- wmean(coord[as.numeric(f)==j,1],weights[as.numeric(f)==j])
		bar.y[[i]][j] <- wmean(coord[as.numeric(f)==j,2],weights[as.numeric(f)==j])
	    }
	  }
	  bar.x <- unlist(bar.x)
	  bar.y <- unlist(bar.y)
	  text(bar.x,bar.y,centr$lev,col=f.col,cex=f.cex)
	}
    }
  }
}

MVA.constr <- function(x,...) {
  UseMethod("MVA.constr")
}

MVA.constr.CCA.vegan <- MVA.constr.RDA.vegan <- function(x,...) {
  constr <- if ("formula" %in% names(x$call)) {
    as.data.frame(model.frame(x))
  } else {
    as.data.frame(eval(x$call$Y))
  }
  type <- logical(ncol(constr))
  for (i in 1:ncol(constr)) {type[i] <- is.factor(constr[,i])}
  res <- if (all(!type)) {"n"} else
    if (all(type)) {"f"} else {"nf"}
  return(res)
}

MVA.constr.CCA.ade4 <- function(x,...) {
  constr <- as.data.frame(eval(x$call$sitenv))
  type <- logical(ncol(constr))
  for (i in 1:ncol(constr)) {type[i] <- is.factor(constr[,i])}
  res <- if (all(!type)) {"n"} else
    if (all(type)) {"f"} else {"nf"}
  return(res)}

MVA.constr.RDA.ade4 <- function(x,...) {
  constr <- as.data.frame(eval(x$call$df))
  type <- logical(ncol(constr))
  for (i in 1:ncol(constr)) {type[i] <- is.factor(constr[,i])}
  res <- if (all(!type)) {"n"} else
    if (all(type)) {"f"} else {"nf"}
  return(res)
}


MVA.centr <- function(x,...) {
  UseMethod("MVA.centr")
}

MVA.centr.CCA.vegan <- MVA.centr.RDA.vegan <- function(x,...) {
  constr <- if ("formula" %in% names(x$call)) {
    as.data.frame(model.frame(x))
  } else {
    as.data.frame(eval(x$call$Y))
  }
  type <- logical(ncol(constr))
  for (i in 1:ncol(constr)) {type[i] <- is.factor(constr[,i])}
  constr <- as.data.frame(constr[,type])
  res <- list(fac=constr,lev=rownames(x$CCA$centroids))
  return(res)
}

MVA.centr.CCA.ade4 <- function(x,...) {
  constr <- as.data.frame(eval(x$call$sitenv))
  type <- logical(ncol(constr))
  for (i in 1:ncol(constr)) {type[i] <- is.factor(constr[,i])}
  constr <- as.data.frame(constr[,type])
  lev <- character()
  for (i in 1:ncol(constr)) {
    lev <- c(lev,paste0(colnames(constr)[i],levels(constr[,i])))
  }
  res <- list(fac=constr,lev=lev)
  return(res)
}

MVA.centr.RDA.ade4 <- function(x,...) {
  constr <- as.data.frame(eval(x$call$df))
  type <- logical(ncol(constr))
  for (i in 1:ncol(constr)) {type[i] <- is.factor(constr[,i])}
  constr <- as.data.frame(constr[,type])
  lev <- character()
  for (i in 1:ncol(constr)) {
    lev <- c(lev,paste0(colnames(constr)[i],levels(constr[,i])))
  }
  res <- list(fac=constr,lev=lev)
  return(res)
}

