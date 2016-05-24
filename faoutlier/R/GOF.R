#' Goodness of Fit Distance
#'
#' Compute Goodness of Fit distances between models when removing the \eqn{i_{th}} case.
#' If mirt is used, then the values will be associated with the unique response patterns instead.
#'
#' Note that \code{GOF} is not limited to confirmatory factor analysis and
#' can apply to nearly any model being studied
#' where detection of influential observations is important.
#'
#' @aliases GOF
#' @param data matrix or data.frame
#' @param model if a single numeric number declares number of factors to extract in
#'   exploratory factor analysis (requires complete dataset, i.e., no missing).
#'   If \code{class(model)} is a sem (semmod), or lavaan (character),
#'   then a confirmatory approach is performed instead. Finally, if the model is defined with
#'   \code{mirt::mirt.model()} then distances will be computed for categorical data with the
#'   mirt package
#' @param M2 logical; use the M2 statistic for when using mirt objects instead of G2?
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#'   \code{\link{gCD}}, \code{\link{LD}}, \code{\link{obs.resid}},
#'   \code{\link{robustMD}}, \code{\link{setCluster}}
#' @references
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data screening and assumption testing for
#' exploratory and confirmatory factor analysis. \emph{Frontiers in Psychology, 3}, 1-21.
#' @keywords cooks
#' @export GOF
#' @examples
#'
#' \dontrun{
#'
#' #run all GOF functions using multiple cores
#' setCluster()
#'
#' #Exploratory
#' nfact <- 3
#' (GOFresult <- GOF(holzinger, nfact))
#' (GOFresult.outlier <- GOF(holzinger.outlier, nfact))
#' plot(GOFresult)
#' plot(GOFresult.outlier)
#'
#' #-------------------------------------------------------------------
#' #Confirmatory with sem
#' model <- sem::specifyModel()
#'	  F1 -> Remndrs,    lam11
#' 	  F1 -> SntComp,    lam21
#' 	  F1 -> WrdMean,    lam31
#' 	  F2 -> MissNum,    lam42
#' 	  F2 -> MxdArit,    lam52
#' 	  F2 -> OddWrds,    lam62
#' 	  F3 -> Boots,      lam73
#'	  F3 -> Gloves,     lam83
#' 	  F3 -> Hatchts,    lam93
#' 	  F1 <-> F1,   NA,     1
#' 	  F2 <-> F2,   NA,     1
#' 	  F3 <-> F3,   NA,     1
#'
#' (GOFresult <- GOF(holzinger, model))
#' (GOFresult.outlier <- GOF(holzinger.outlier, model))
#' plot(GOFresult)
#' plot(GOFresult.outlier)
#'
#' #-------------------------------------------------------------------
#' #Confirmatory with lavaan
#' model <- 'F1 =~  Remndrs + SntComp + WrdMean
#' F2 =~ MissNum + MxdArit + OddWrds
#' F3 =~ Boots + Gloves + Hatchts'
#'
#' (GOFresult <- GOF(holzinger, model, orthogonal=TRUE))
#' (GOFresult.outlier <- GOF(holzinger.outlier, model, orthogonal=TRUE))
#' plot(GOFresult)
#' plot(GOFresult.outlier)
#'
#'
#' # categorical data with mirt
#' library(mirt)
#' data(LSAT7)
#' dat <- expand.table(LSAT7)
#' model <- mirt.model('F = 1-5')
#' result <- GOF(dat, model)
#' plot(result)
#'
#' }
GOF <- function(data, model, M2 = TRUE, ...)
{
    G2 <- function(mirt_mod){
        F <- mirt_mod@Data$Freq[[1L]]
        F <- F[F != 0L]
        N <- sum(F)
        P <- mirt_mod@Internals$Pl
        2 * sum(F * log(F / (N*P)))
    }
    fM2 <- function(mirt_mod) mirt::M2(mirt_mod)$M2

    f_numeric <- function(ind, data, model, ...){
        tmp <- factanal(data[-ind, ], model, ...)
        tmp$STATISTIC
    }
    f_sem <- function(ind, data, model, objective, ...){
        tmp <- sem::sem(model, data=data[-ind, ], objective=objective, ...)
        tmp$criterion * (tmp$N - 1)
    }
    f_lavaan <- function(ind, data, model, ...){
        tmp <- lavaan::sem(model, data=data[-ind, ], ...)
        tmp@Fit@test[[1L]]$stat
    }
    f_mirt <- function(ind, data, model, large, sv, M2, ...){
        large$Freq[[1L]][ind] <- large$Freq[[1L]][ind] - 1L
        tmp <- mirt::mirt(data=data, model=model, verbose=FALSE, large=large, pars=sv, ...)
        ret <- if(M2) fM2(tmp)
            else G2(tmp)
        ret
    }

	rownames(data) <- 1:nrow(data)
    index <- matrix(1L:nrow(data))
	if(is.numeric(model)){
	    if(any(is.na(data)))
	        stop('Numeric model requires complete dataset (no NA\'s)')
		MLmod <- factanal(data, model, ...)$STATISTIC
		LR <- myApply(index, MARGIN=1L, FUN=f_numeric, data=data, model=model, ...)
	} else if(class(model) == "semmod"){
        objective <- if(any(is.na(data))) sem::objectiveFIML else sem::objectiveML
	    MLmod <- sem::sem(model, data=data, objective=objective, ...)
        MLmod <- MLmod$criterion * (MLmod$N - 1)
	    LR <- myApply(index, MARGIN=1L, FUN=f_sem, data=data, model=model, objective=objective, ...)
	} else if(class(model) == "character"){
        MLmod <- lavaan::sem(model, data=data, ...)
        MLmod <- MLmod@Fit@test[[1]]$stat
        LR <- myApply(index, MARGIN=1L, FUN=f_lavaan, data=data, model=model, ...)
	} else if(class(model) == "mirt.model"){
        large <- MLmod_full <- mirt::mirt(data=data, model=model, large = TRUE)
        index <- matrix(1L:length(large$Freq[[1L]]))
        MLmod_full <- mirt::mirt(data=data, model=model, verbose = FALSE, large=large, ...)
        sv <- mirt::mod2values(MLmod_full)
        MLmod <- if(M2) fM2(MLmod_full)
            else G2(MLmod_full)
        LR <- myApply(index, MARGIN=1L, FUN=f_mirt, data=data, model=model,
                      sv=sv, large=large, M2=M2, ...)
	} else {
        stop('model class not supported')
    }
	deltaX2 <- LR - MLmod
	if(class(model) != "mirt.model")
	    names(deltaX2) <- rownames(data)
	class(deltaX2) <- 'GOF'
	deltaX2
}

#' @rdname GOF
#' @param x an object of class \code{GOF}
#' @param ncases number of extreme cases to display
#' @param digits number of digits to round in the printed result
#' @param ... additional parameters to be passed
#' @export
print.GOF <- function(x, ncases = 10, digits = 5, ...)
{
	sorted <- x[order(abs(x), decreasing = TRUE)]
	if(ncases %% 2 != 0) ncases <- ncases + 1
	sorted <- c(sorted[1:(ncases/2)],
		sorted[(length(sorted)-(ncases/2 + 1)):length(sorted)])
	ret <- matrix(sorted)
	rownames(ret) <- names(sorted)
	colnames(ret) <- 'GOF'
	return(print(round(ret, digits)))
}

#' @rdname GOF
#' @param y a \code{NULL} value ignored by the plotting function
#' @param type type of plot to use, default displays points and lines
#' @param main the main title of the plot
#' @param ylab the y label of the plot
#' @param absolute logical; use absolute values instead of deviations?
#' @export
plot.GOF <- function(x, y = NULL, main = 'Goodness of Fit Distance',
	type = c('p','h'), ylab = 'GOF', absolute = FALSE, ...)
{
	GOF <- if(absolute) abs(as.numeric(x)) else as.numeric(x)
	ID <- 1:length(x)
	dat <- data.frame(GOF,ID)
	ret <- lattice::xyplot(GOF~ID, dat, type = type, main = main, ylab = ylab,
	                       panel = function(...){
                               panel.xyplot(ID, GOF, type='h')
                               panel.xyplot(ID, GOF, type='p')
                               panel.abline(0)
                               }, ...)
	return(ret)
}

