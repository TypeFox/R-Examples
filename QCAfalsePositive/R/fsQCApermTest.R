

#' A Simple Permutation Test for Type I Error in fsQCA
#'
#' A permutation test for fuzzy-set qualitative comparative analysis (fsQCA),
#' designed to calculate the probability of a false positive given the number of
#' hypotheses implicitly tested and the distribution of the data.
#' @param y The outcome variable of interest.
#' @param configs A list of configurations to be tested against \code{y}.
#' @param total.configs The total number of configurations used in the original
#' fsQCA analysis. This will generally equal the number of lines in the
#' truth table used for Boolean minimization.
#' @param num.iter The number of iterations to use for the permutation test.
#' Larger numbers of iterations result in more precise p-values.
#' @param my.seed The seed used to generate random numbers.
#' @param adj.method The method used to calculate adjusted p-values (see 
#' \code{\link{p.adjust}} for details).
#' @return An object containing the aggregate results of the permutation test
#' as well as the individual permutations.
#' @keywords permutation fsQCA multiple inference p-value adjust
#' @export
#' @examples
#' data(social.revolutions)
#' attach(social.revolutions)
#'   
#' intersect <- pmin(breakdown, pop.ins)
#' intersect2 <- pmin(breakdown, (1-pop.ins))
#' intersect3 <- pmin((1-breakdown), pop.ins)
#' intersect4 <- pmin((1-breakdown), (1-pop.ins))
#'   
#' test <- fsQCApermTest(y=soc.rev, configs=list(BI=intersect, Bi=intersect2, 
#'    bI=intersect3, bi=intersect4), total.configs=4)
#' summary(test)
#' plot(test)





fsQCApermTest <- function(y, configs, total.configs, num.iter=10000, my.seed=123, adj.method="holm"){
	if(typeof(configs)!="list"){stop("configs must be a list of the form list(name1=configuration1, name2=configuration2, ...)")}
	if(is.null(names(configs))){stop("configs must be a list of the form list(name1=configuration1, name2=configuration2, ...)")}
	if(total.configs < length(configs)){stop("Value of total.configs cannot be less than the length of configs.")}
	out.mat <- matrix(, nrow = length(configs), ncol = 14)
	permutes.cex <- list()
	permutes.con <- list()
	for(j in 1:length(configs)){
		counterex <- NULL
		consistency <- NULL
		necessary.condition <- sum(y<configs[[j]])>(length(y)/2)
		set.seed(my.seed)
		message(c("Generating permutations for ", names(configs)[j], "..."))
		if(necessary.condition){
			message("  Lower-triangular configuration detected...")
			} else{
				message("   Upper-triangular configuration detected...")
				}
		for(i in 1:num.iter){
			new.y <- sample(y, size=length(y), replace=F)
			if(necessary.condition){
			consistency <- c(consistency, sum(pmin(configs[[j]], new.y))/sum(new.y))
			counterex <- c(counterex, sum(configs[[j]] < new.y))
			} else{
			consistency <- c(consistency, sum(pmin(configs[[j]], new.y))/sum(configs[[j]]))
			counterex <- c(counterex, sum(configs[[j]] > new.y))
			}
			}
		
		if(necessary.condition){
		obs.counterexamples <- sum(configs[[j]] < y)
		} else{
		obs.counterexamples <- sum(configs[[j]] > y)
		}
		counterex.p <- sum(counterex<=obs.counterexamples)/length(counterex)
	
		if(necessary.condition){
			obs.consistency <- sum(pmin(configs[[j]], y))/sum(y)
		} else{
			obs.consistency <- sum(pmin(configs[[j]], y))/sum(configs[[j]])
		}
		consistency.p <- sum(consistency>=obs.consistency)/length(consistency)
		
		out.mat[j, 1] <- obs.counterexamples
		out.mat[j, 2] <- sort(counterex)[length(counterex)*0.05]
		out.mat[j, 3] <- max(counterex)
		out.mat[j, 4] <- counterex.p
		out.mat[j, 8] <- obs.consistency
		out.mat[j, 9] <- min(consistency)
		out.mat[j, 10] <- sort(consistency)[length(consistency)*0.95]
		out.mat[j, 11] <- consistency.p
		permutes.cex[[j]] <- counterex
		permutes.con[[j]] <- consistency
	}
	out.mat[,5] <- p.adjust(c(out.mat[,4], rep(1,(total.configs-length(configs)))), method=adj.method)[1:length(configs)]
	out.mat[,6] <- p.threshold.adjust(total.configs, adj.method)[rank(out.mat[,4])][1:length(configs)]
	out.mat[,7] <- sqrt((out.mat[,5]*(1-out.mat[,5]))/length(configs[[1]]))
	out.mat[,12] <- p.adjust(c(out.mat[,11], rep(1,(total.configs-length(configs)))), method=adj.method)[1:length(configs)]
	out.mat[,13] <- p.threshold.adjust(total.configs, adj.method)[rank(out.mat[,10])][1:length(configs)]
	out.mat[,14] <- sqrt((out.mat[,12]*(1-out.mat[,12]))/length(configs[[1]]))
	
	for(j in 1:length(configs)){
	out.mat[j, 2] <- sort(permutes.cex[[j]])[length(permutes.cex[[j]])*out.mat[j,6]]
	out.mat[j, 10] <- sort(permutes.con[[j]])[length(permutes.con[[j]])*(1-out.mat[j,13])]
		}
	
	rownames(out.mat) <- names(configs)
	colnames(out.mat) <- c("Observed", "Lower c.i.", "Upper Bound", "p-raw", "p-adj", "p-val-adj", "se(p-adj)", "Observed", "Lower Bound", "Upper c.i.", "p-raw", "p-adj", "p-val-adj", "se(p-adj)")
	names(permutes.cex) <- names(configs)
	names(permutes.con) <- names(configs)
	return.obj <- list(call = match.call(), total.configurations = total.configs, config.names = names(configs), p.adj.method = adj.method, result.cex = out.mat[,1:7], result.con = out.mat[,8:14], permutations.cex = permutes.cex, permutations.con = permutes.con)
	class(return.obj) <- "fsQCApt"
	return.obj
}
