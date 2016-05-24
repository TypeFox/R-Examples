#' Summarizing InSilicoVA Model Fits
#' 
#' This function is the summary method for class \code{insilico}.
#' 
#' \code{summary.insilico} formats some basic information about the InSilicoVA
#' fitted object on screen and show the several top CSMFs of user's choice. See
#' below for more detail.
#' 
#' @aliases summary.insilico print.insilico_summary
#' @param object Fitted \code{"insilico"} object.
#' @param CI.csmf Confidence interval for CSMF estimates.
#' @param CI.cond Confidence interval for conditional probability estimates
#' @param file Optional .csv file to write to. If it is specified, individual
#' cause of death distribution will be saved to the file.
#' @param top Number of top causes to display on screen.
#' @param id ID of specific death to display on screen.
#' @param \dots Not used.
#' @return \item{id.all}{ all IDs of the deaths.} 
#' \item{indiv}{ individual Cause of Death distribution matrix.} 
#' \item{csmf}{ CSMF distribution and confidence interval for each cause.} 
#' \item{csmf.ordered}{ CSMF distribution and confidence interval for each cause, ordered by mean.} 
#' \item{condprob}{ Conditional probability matrix and confidence intervals.}
#' \item{updateCondProb}{Component of \code{"insilico"} object.}
#' \item{keepProbbase.level}{Component of \code{"insilico"} object.}
#' \item{datacheck}{Component of \code{"insilico"} object.}
#' \item{Nsim}{Component of \code{"insilico"} object.}
#' \item{thin}{Component of \code{"insilico"} object.} \item{burnin}{Component
#' of \code{"insilico"} object.} 
#' \item{jump.scale}{Component of \code{"insilico"} object.}
#' \item{levels.prior}{Component of \code{"insilico"} object.}
#' \item{levels.strength}{Component of \code{"insilico"} object.}
#' \item{trunc.min}{Component of \code{"insilico"} object.}
#' \item{trunc.max}{Component of \code{"insilico"} object.}
#' \item{subpop_counts}{Component of \code{"insilico"} object.}
#' \item{showTop}{Component of \code{"insilico"} object.}
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{insilico}}, \code{\link{plot.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death
#' assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042}
#' \url{http://arxiv.org/abs/1411.3042} (2014)
#' @examples
#' \dontrun{
#' # load sample data together with sub-population list
#' data(RandomVA1)
#' # extract InterVA style input data
#' data <- RandomVA1$data
#' # extract sub-population information. 
#' # The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
#' subpop <- RandomVA1$subpop
#' 
#' # run without subpopulation
#' fit1<- insilico( data, subpop = NULL, 
#'               Nsim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE)
#' summary(fit1)
#' summary(fit1, top = 10)
#'
#' # save individual COD distributions to files
#' summary(fit1, file = "results.csv")
#' }
#' @export summary.insilico
summary.insilico <- function(object, CI.csmf = 0.95, CI.cond = 0.95, file = NULL, top = 10, id = NULL, ...){
	id.all <- object$id
	prob <- object$indiv.prob
	csmf <- object$csmf
	if(is.null(object$conditional.probs)){
		cond.prob = FALSE
	}else{
		cond.prob = TRUE
	}
	if(cond.prob){
		cond <- object$conditional.probs
		if(length(dim(cond)) == 2){
			isLevel = TRUE
		}else{
			isLevel = FALSE
		}
	}

	indiv <- cbind(id.all, prob)
	## write individual COD distribution to file
	if(!is.null(file)){
		colnames(indiv)[1] <- "ID"
		write.csv(indiv, file = file, row.names = FALSE)
	}

	# organize CSMF
	ci.low <- (1 - CI.csmf) / 2
	ci.up <- 1 - ci.low
	if(class(object$csmf) == "list"){
		csmf.out <- vector("list", length(csmf))
		csmf.out.ordered <- vector("list", length(csmf))
		for(i in 1:length(csmf)){
			mean <- apply(csmf[[i]], 2, mean)
			median <- apply(csmf[[i]], 2, median)
			sd <- apply(csmf[[i]], 2, sd)
			low <- apply(csmf[[i]], 2, function(object){quantile(object, prob = ci.low)})
			up <- apply(csmf[[i]], 2, function(object){quantile(object, prob = ci.up)})
			csmf.out[[i]]  <- cbind(mean, sd,  low, median, up)
			colnames(csmf.out[[i]]) <- cbind("Mean","Std.Error", "Lower", "Median", "Upper")
			csmf.out.ordered[[i]] <- csmf.out[[i]][order(csmf.out[[i]][,1],
												decreasing = TRUE), ]
			names(csmf.out)[i] <- names(csmf)[i]
			names(csmf.out.ordered)[i] <- names(csmf)[i]
		}
	}else{
		mean <- apply(csmf, 2, mean)
		median <- apply(csmf, 2, median)
		sd <- apply(csmf, 2, sd)
		low <- apply(csmf, 2, function(object){quantile(object, prob = ci.low)})
		up <- apply(csmf, 2, function(object){quantile(object, prob = ci.up)})
		csmf.out <- cbind(mean, sd, low, median, up)
		colnames(csmf.out) <- cbind("Mean","Std.Error","Lower",  "Median", "Upper")
		csmf.out.ordered <- csmf.out[order(csmf.out[,1], decreasing = TRUE), ]
	}
	
	

	# organize conditional probability matriobject
	ci.low <- (1 - CI.cond) / 2
	ci.up <- 1 - ci.low
	
	if(cond.prob){
		if(isLevel){
			mean <- apply(cond, 2, mean)
			median <- apply(cond, 2, median)
			sd <- apply(cond, 2, sd)
			low <- apply(cond, 2, function(object){quantile(object, prob = ci.low)})
			up <- apply(cond, 2, function(object){quantile(object, prob = ci.up)})
			cond.out <- cbind(mean, sd, low, median, up)
			colnames(cond.out) <- cbind("Mean","Std.Error",  "Lower", "Median","Upper")
			rownames(cond.out) <- colnames(cond)
		}else{
			mean <- t(apply(cond, c(2,3), mean))
			median <- t(apply(cond, c(2,3), median))
			sd <- t(apply(cond, c(2,3), sd))
			low <- t(apply(cond, c(2,3), function(object){quantile(object, probs =ci.low)}))
			up <- t(apply(cond, c(2,3), function(object){quantile(object, probs =ci.up)}))
			cond.out <- list(mean, sd, low,median, up)
			names(cond.out) <- c("Mean","Std.Error", "Lower", "Median", "Upper")
		}
	}else{
		cond.out <- NULL
	}

	# get subpopulation counts table
	if(!is.null(object$subpop)){
		subpop_counts <- table(object$subpop)
		subpop_counts <- subpop_counts[match(names(object$csmf), 
											 names(subpop_counts))]
	}else{
		subpop_counts <- NULL
	}

	if(!is.null(id)){
		if(!id %in% rownames(object$indiv.prob)){
			stop("Invalid ID, not exist in data.")
		}
		whichtoprint <- order(object$indiv.prob[id, ], decreasing = TRUE)[1:top]

		if(is.null(object$indiv.prob.lower)){
			warning("C.I. for individual probabilities have not been calculated. Please use get.indiv() function to update C.I. first.\n")
			indiv.prob <- cbind(object$indiv.prob[id, whichtoprint], NA, NA, NA)

		}else{		
		indiv.prob <- cbind(object$indiv.prob[id, whichtoprint], 
							object$indiv.prob.lower[id, whichtoprint], 
							object$indiv.prob.median[id, whichtoprint], 
							object$indiv.prob.upper[id, whichtoprint])

		}
		colnames(indiv.prob) = c("Mean", "Lower", "Median", "Upper")	
	}else{
		indiv.prob <- NULL
	}
	out <- list( id.all = id.all, 
				 indiv = indiv, 
				 csmf = csmf.out, 
				 csmf.ordered = csmf.out.ordered,
				 condprob = cond.out, 				
				 updateCondProb = object$updateCondProb, 
				 keepProbbase.level = object$keepProbbase.level, 
				 datacheck = object$datacheck,
				 Nsim = object$Nsim, 
				 thin = object$thin, 
				 burnin = object$burnin, 
				 HIV = object$HIV, 
				 Malaria = object$Malaria, 
				 jump.scale = object$jump.scale, 
				 levels.prior = object$levels.prior, 
				 levels.strength = object$levels.strength, 
				 trunc.min = object$trunc.min, 
				 trunc.maobject = object$trunc.max, 
			     subpop_counts = subpop_counts,
				 showTop = top, 
				 id = id, 
				 indiv.prob = indiv.prob, 
				 indiv.CI = object$indiv.CI)
	class(out) <- "insilico_summary"
	return(out)
}
