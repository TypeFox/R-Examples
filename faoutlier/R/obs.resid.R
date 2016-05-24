#' Model predicted residual outliers
#'
#' Compute model predicted residuals for each variable using regression
#' estimated factor scores.
#'
#' @aliases obs.resid
#' @param data matrix or data.frame
#' @param model if a single numeric number declares number of factors to extract in
#'   exploratory factor analysis. If \code{class(model)} is a sem (semmod), or lavaan (character),
#'   then a confirmatory approach is performed instead
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#'   \code{\link{gCD}}, \code{\link{LD}}, \code{\link{robustMD}}
#' @references
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data screening and assumption testing for
#' exploratory and confirmatory factor analysis. \emph{Frontiers in Psychology, 3}, 1-21.
#' @keywords covariance
#' @export obs.resid
#' @examples
#'
#' \dontrun{
#' data(holzinger)
#' data(holzinger.outlier)
#'
#' #Exploratory
#' nfact <- 3
#' (ORresult <- obs.resid(holzinger, nfact))
#' (ORresult.outlier <- obs.resid(holzinger.outlier, nfact))
#' plot(ORresult)
#' plot(ORresult.outlier)
#'
#' #-------------------------------------------------------------------
#' #Confirmatory with sem
#' model <- sem::specifyModel()
#'    F1 -> Remndrs,    lam11
#' 	  F1 -> SntComp,    lam21
#' 	  F1 -> WrdMean,    lam31
#' 	  F2 -> MissNum,    lam41
#' 	  F2 -> MxdArit,    lam52
#' 	  F2 -> OddWrds,    lam62
#' 	  F3 -> Boots,      lam73
#'	  F3 -> Gloves,     lam83
#' 	  F3 -> Hatchts,    lam93
#' 	  F1 <-> F1,   NA,     1
#' 	  F2 <-> F2,   NA,     1
#' 	  F3 <-> F3,   NA,     1
#'
#' (ORresult <- obs.resid(holzinger, model))
#' (ORresult.outlier <- obs.resid(holzinger.outlier, model))
#' plot(ORresult)
#' plot(ORresult.outlier)
#'
#' #-------------------------------------------------------------------
#' #Confirmatory with lavaan
#' model <- 'F1 =~  Remndrs + SntComp + WrdMean
#' F2 =~ MissNum + MxdArit + OddWrds
#' F3 =~ Boots + Gloves + Hatchts'
#'
#' (obs.resid2 <- obs.resid(holzinger, model, orthogonal=TRUE))
#' (obs.resid2.outlier <- obs.resid(holzinger.outlier, model, orthogonal=TRUE))
#' plot(obs.resid2)
#' plot(obs.resid2.outlier)
#'
#' }
obs.resid <- function(data, model, ...)
{
	ret <- list()
	rownames(data) <- 1:nrow(data)
	N <- nrow(data)
	if(any(is.na(data)))
	    stop('All routines require complete datasets (no NA\'s)')
    if(!all(sapply(data, function(x) is.numeric(x) || is.integer(x))))
        stop('data must conist of numeric elements')
	if(is.numeric(model)){
		R <- cor(data)
		mod <- factanal(data, model, rotation='none', scores = 'regression', ...)
		scores <- mod$scores
		ret$fascores <- scores
		Lambda <- unclass(mod$loadings)
		Theta <- diag(mod$uniquenesses)
		e <- data - scores %*% t(Lambda)
		VAR <- Theta %*% solve(R) %*% Theta
		eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(scale(e, scale = FALSE)))
		colnames(eji) <- colnames(e) <- colnames(data)
		ret$res <- e
		ret$std_res <- eji
	} else if(class(model) == "semmod"){
        C <- cov(data)
        vnames <- colnames(C)
        mod <- sem::sem(model, C, N, ...)
        scores <- sem::fscores(mod, data)
        ret$fascores <- scores
        lnames <- setdiff(colnames(mod$P), vnames)
        Lambda <- mod$A[1:length(vnames), (length(vnames)+1):ncol(mod$A)]
        Theta <- mod$P[1:length(vnames),1:length(vnames)]
        e <- data - scores %*% t(Lambda)
        VAR <- Theta %*% solve(C) %*% Theta
        eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(scale(e, scale = FALSE)))
        colnames(eji) <- colnames(e) <- colnames(data)
        ret$res <- e
        ret$std_res <- eji
	} else if(class(model) == "character"){
        dots <- list(...)
        if(!is.null(dots$ordered))
            stop('observed residuals only defined for continuous data')
        mod <- lavaan::sem(model, data=data, ...)
        C <- mod@SampleStats@cov[[1L]]
        scores <- lavaan::predict(mod)
        ret$fascores <- scores
        Lambda <- mod@Model@GLIST$lambda
        Theta <- mod@Model@GLIST$theta
        Psi <- mod@Model@GLIST$psi
        dat <- data[,mod@Data@ov.names[[1L]]]
        e <- dat - scores %*% t(Lambda)
        VAR <- Theta %*% solve(C) %*% Theta
        eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(scale(e, scale = FALSE)))
        colnames(eji) <- colnames(e) <- colnames(dat)
        ret$res <- e
        ret$std_res <- eji
	} else {
	    stop('model class not supported')
	}
	ret$id <- rownames(data)
	class(ret) <- 'obs.resid'
	ret
}

#' @rdname obs.resid
#' @param x an object of class \code{obs.resid}
#' @param restype type of residual used, either \code{'obs'} for observation value
#'   (inner product), \code{'res'} or \code{'std_res'} for unstandardized and standardized
#'   for each variable, respectively
#' @param ... additional parameters to be passed
#' @export
print.obs.resid <- function(x, restype = 'obs', ...)
{
    if(restype == 'res') return(print(x$res))
    if(restype == 'std_res') return(print(x$std_res))
	stat <- c()
	for(i in 1:length(x$id))
		stat[i] <- x$std_res[i, ] %*% x$std_res[i, ]
    if(restype == 'obs') return(stat)
}

#' @rdname obs.resid
#' @param y a \code{NULL} value ignored by the plotting function
#' @param main the main title of the plot
#' @param type type of plot to use, default displays points and lines
#' @export
plot.obs.resid <- function(x, y = NULL, main = 'Observed Residuals',
	type = c('p','h'), restype = 'obs', ...)
{
	ylab <- switch(restype,
		obs = 'Observed residuals',
		res = 'Observed variable residuals',
		std_res = 'Standardized variable residuals')
	ID <- as.numeric(x$id)
	stat <- c()
	for(i in 1:length(x$id))
		stat[i] <- x$std_res[i, ] %*% x$std_res[i, ]
	if(restype == 'obs'){
		dat <- data.frame(stat=stat,ID=ID)
		ret <- lattice::xyplot(stat~ID, dat, type = type, main = main, ylab = ylab, ...)
	}
	if(restype == 'res' || restype == 'std_res'){
		if(restype == 'res') dat <- data.frame(ID=ID,x$res)
			else dat <- data.frame(ID=ID,x$std_res)
		rownames(dat) <- ID
		dat2 <- reshape(dat, v.names='res', direction = 'long', varying=2:ncol(dat),
                        times = colnames(dat)[-1], timevar='variable', sep='')
		ret <- lattice::xyplot(res~ID|as.factor(variable), dat2, type = type, main = main,
                               ylab = ylab, ...)
	}
	return(ret)
}
