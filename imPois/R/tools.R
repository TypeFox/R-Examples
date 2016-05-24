#' @rdname characterize_imprecise_prior
#' @title Characterize Imprecise Prior
#' @description A set of linear inequalities is used for characterizing a polygonal convex set.  The \code{chull} function in the \code{grDevices} package and the \code{convhulln} function in the \code{geometry} package are used to search for extreme points constructing a convex hull. 
#' @param ui constraint matrix (k x p), see below.
#' @param ci constrain vector of length k, see below.
#' @param pmat matrix (k x p) containig coordinate information in d-dimensions.
#' @examples 
#' # lc0 <- list(lhs=rbind(diag(2), -diag(2)), rhs=c(0,0,-1,-1))
#' # op <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
#' # op <- iprior(ui=rbind(c(1,0),c(0,1),c(-1,-1)), ci=c(0,0,-5)) 
#' op <- iprior(ui=rbind(c(1,0),c(0,1),c(0,-1),c(1,1),c(-2,-1)), 
#'    ci=c(1,2,-8,5,-14)) # (3,8),(1,8), (1,4),(3,2)(6,2)
#' 
#' @author Chel Hee Lee <gnustats@@gmail.com>
#' @export 
iprior <- function(ui, ci, pmat){ # in 2 dimensions
	
	if (missing(ui)) ui <- NULL
	if (missing(ci)) ci <- NULL
	if (missing(pmat)) pmat <- NULL
	
	robj <- list()
	
	if (!is.null(pmat) & (is.null(ui) | is.null(ci))) {
		if (!is.null(ui) | !is.null(ci)) message("ui or ci is ignored; instead, uses NA.")
		if (ncol(pmat) == 1) stop("Matrix with at least two columns is needed.")
		if (ncol(pmat) == 2) {
			imat <- grDevices::chull(pmat)
			vtx <- pmat[imat,]
		}
		if (ncol(pmat) >= 3) {
			robj$vtx0 <- pmat
			imat <- t(geometry::convhulln(pmat, options="Tv"))
			robj$imat <- imat
			vtx <- pmat[imat,]
			vtx <- vtx[!duplicated(vtx),]
		}
	}
	
	if(!is.null(ui) & !is.null(ci)){
		if( !is.null(pmat) ) message("pmat is ignored.")
		pmat <- NULL
	
		allsys <- utils::combn(nrow(ui), ncol(ui))
		sol <- list()
		for (i in 1:ncol(allsys)) {
			idx <- c(allsys[,i])
			A <- ui[idx,]
			b <- ci[idx]
			sol[[i]] <- tryCatch(solve(A,b), error=function(e) return(NaN))
		}
		roots <- do.call(cbind, sol)
		idx <- which(colSums(ui %*% roots >= ci) == nrow(ui))
		vtx <- t(roots[,idx])
		
		if (ncol(vtx) == 2) {
			vtx <- vtx[grDevices::chull(vtx),]
		}
		
		if (ncol(vtx) >= 3) {
			imat <- t(geometry::convhulln(vtx, options="Tv"))
			robj$imat <- imat
			robj$vtx0 <- vtx
			vtx <- vtx[imat,]
			vtx <- vtx[!duplicated(vtx),]
		}
		robj$ui <- ui
		robj$ci <- ci
		robj$sol <- sol
		robj$roots <- roots
		robj$idx <- idx 
	} # end of if ui and ci
	
	rownames(vtx) <- paste("p", 1:nrow(vtx), sep="")
	colnames(vtx) <- paste("d", 1:ncol(vtx), sep="")
	
	robj$vtx <- vtx
	class(robj) <- c("impinf", "iprior")
	return(robj)
}


#' @rdname imprecise_learning_from_observation
#' @title Applying Bayes Rule
#' @param object an object for which the Bayes rule is needed. 
#' @param y a vector of observations
#' @param ztrunc logical value indicating truncation at zero.
#' @param ... further arguments passed to methods
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export 
update.impinf <- function(object, y, ztrunc=FALSE, ...){ 
	
	stopifnot(inherits(object, "impinf"))
	
	vtx <- object$vtx # matrix
	if (ncol(vtx) == 2) vtx1 <- cbind(vtx[,1]+sum(y), vtx[,2]+length(y))
	else vtx1 <- cbind(vtx[,1], vtx[,2]+sum(y), vtx[,3]+length(y))
	colnames(vtx1) <- colnames(vtx)
	
	if(!ztrunc){ 
		op <- lapply( as.list(as.data.frame(t(vtx))), function(x){
			if (length(x) == 2) tryCatch(evfn(y=y, pars=c(0,x[1],x[2]))$value, error=function(e) return(NaN))
			else tryCatch(evfn(y=y, pars=c(x[1],x[2],x[3]))$value, error=function(e) return(NaN))
		})
	} 
	
	if (ztrunc) {
		op <- lapply(as.list(as.data.frame(t(vtx))), function(x){
			if (length(x) == 2) tryCatch(evfn.ztrunc(y=y, pars=c(0,x[1],x[2]))$value, error=function(e) return(NaN))
			else tryCatch(evfn.ztrunc(y=y,pars=c(x[1],x[2],x[3]))$value, error=function(e) return(NaN))
		})
	}
	robj <- list(y=y, ztrunc=ztrunc, vtx=vtx, vtx1=vtx1, ev=unlist(op))
	class(robj) <- c("impinf", "update")
	return(robj)
}

#' @rdname summary_imprecise_object
#' @title Summary of \code{impinf} object
#' @param object an object for which a summary is needed.
#' @param ... additional arguments affecting the summary produced.
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export
summary.impinf <- function(object, ...){
	stopifnot(inherits(object, "impinf"))
	
	ztrunc <- object$ztrunc
	vtx <- object$vtx
	vtx1 <- object$vtx1
	ev <- object$ev
	y <- object$y
	
	inf <- ev[which.min(ev)]
	sup <- ev[which.max(ev)]
	delta <- sup - inf
	names(delta) <- NULL
	
	inf.p <- vtx[which.min(ev), ]
	sup.p <- vtx[which.max(ev), ]
	inf.p1 <- vtx1[which.min(ev), ]
	sup.p1 <- vtx1[which.max(ev), ]
	
	robj <- list(y=y, ztrunc=ztrunc, inf=inf, sup=sup, delta=delta, inf.p=inf.p, sup.p=sup.p, inf.p1=inf.p1, sup.p1=sup.p1, ev=ev)
	class(robj) <- c("summary.impinf")
	invisible(robj)
}

#' @rdname print_imprecise_objects
#' @title Print Imprecise Objects
#' @param x an object used to select a method 
#' @param ... further arguments passed to or from other methods 
#' @export
print.summary.impinf <- function(x, ...){
	# stopifnot(inherits(x, "impinf"), inherits(x, "summary"))
	stopifnot(inherits(x, "summary.impinf"))
	
	cat("  zero-truncated Poisson = ", x$ztrunc, "\n\n")
	m <- c(x$inf, x$sup, x$delta)	
	names(m) <- c("min.E(X|y)", "max.E(X|y)", "imprecision")
	print(m)
	
	cat("\n  min.E(X|y) and max.E(X|y) are found \nat", names(x$inf), "and", names(x$sup), "respectively.\n\n  Coordinates:\n" )
	m <- rbind(x$inf.p, x$sup.p)
	rownames(m) <- c(names(x$inf), names(x$sup))
	print(m)
	
	if(any(is.na(x$ev))) message("\nNotes:\nNaN is produced from one of evaluation points.\nmin.E(X|y) or max.E(X|y) may be incorrect.\n")
	
}

# TODO:
# model <- function(formula, data, dist=c("pois", "binom", "geo", "exp", "multinom"), ztrunc=FALSE){
# }
