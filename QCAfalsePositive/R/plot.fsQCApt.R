
#' Plot Results of fsQCA Permutation Test
#'
#' Plots distributions of consistencies and counterexamples from permutation tests 
#' of fsQCA data, including confidence intervals adjusted to account for multiple 
#' inference. Also prints observed consistency values and number of counterexamples 
#' as black dots along the x-axis, for comparison.
#' @param x Object returned by \code{\link{fsQCApermTest}}.
#' @param y A vector of configurations to examine. Default behavior is 
#' to examine all configurations.
#' @param statistic The statistic to examine (\code{consistency}, 
#' \code{counterexamples}, or \code{both}).
#' @param ... Additional parameters to pass on.
#' @return Plots of distributions of consistencies, counterexamples, or both.
#' @keywords fsQCA permutation test distribution
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
#' plot(test)
#' plot(test, "bi", statistic="consistency")
#' plot(test, c("BI", "Bi"), statistic="both")
#' plot(test, statistic="consistency")
#' plot(test, "BI")





plot.fsQCApt <- function(x, y=x$config.names, statistic="both", ...){
	if(!all(y %in% x$config.names)){
		stop("Unknown configuration of variables.")
		}
	if(!(statistic %in% c("consistency", "counterexamples", "both"))){
		stop("statistic must equal 'consistency', 'counterexamples', or 'both'.")
		}
	if(dev.cur()==1){
		dev.new(width=(5+(2*(statistic=="both"))), height=min(7, length(y)*3), units="in", res=360)
		}
	par(mar = c(4,2,2,1))
	config.num <- match(y, x$config.names)
	config.num <- config.num[!is.na(config.num)]
	if(statistic=="both"){
		par(mfrow=c(length(y),2))
		}
		else {
			par(mfrow=c(length(y),1))
			}
	for(m in config.num){
	if(m==max(config.num)){
		conlab <- "Consistency"
		cexlab <- "Counterexamples"
		} else{
		conlab <- "   "
		cexlab <- "   "
		}
	if((statistic=="consistency") | (statistic=="both")){
		pctile <- sort(x$permutations.con[[m]])[length(x$permutations.con[[m]])*(1-x$result.con[m,6])]
		plot(density(x$permutations.con[[m]], cut=0), axes=FALSE, col="white", xlab=conlab, ylab=" ", xlim=c(min(unlist(lapply(x$permutations.con, FUN=min))), 1), main=x$config.names[m])
		axis(1)
		polygon(c(density(x$permutations.con[[m]], cut=0)$x, rev(density(x$permutations.con[[m]], cut=0)$x)), c(density(x$permutations.con[[m]], cut=0)$y, rep(0, length(density(x$permutations.con[[m]], cut=0)$y))), col="#6baed6", border="white", lwd=0.1)
		polygon(c(density(x$permutations.con[[m]], cut=0)$x[density(x$permutations.con[[m]], cut=0)$x>pctile], rev(density(x$permutations.con[[m]], cut=0)$x[density(x$permutations.con[[m]], cut=0)$x>pctile])), c(density(x$permutations.con[[m]], cut=0)$y[density(x$permutations.con[[m]], cut=0)$x>pctile], rep(0, length(density(x$permutations.con[[m]], cut=0)$y[density(x$permutations.con[[m]], cut=0)$x>pctile]))), col="#08519c", border="white", lwd=0.1)	
		points(x$result.con[m,1], par("usr")[3], pch=19, xpd=TRUE)
	}
	if((statistic=="counterexamples") | (statistic=="both")){
		histpct <- sort(x$permutations.cex[[m]])[length(x$permutations.cex[[m]])*(x$result.cex[m,6])]
		hist(x$permutations.cex[[m]], breaks=seq(0,max(unlist(lapply(x$permutations.cex, FUN=max)))+1,1)-0.5, axes=FALSE, col=c(rep("#08519c", histpct), rep("#6baed6", max(x$permutations.cex[[m]])-histpct+1)), border="white", xlab=cexlab, ylab=" ", xlim=c(-1, max(unlist(lapply(x$permutations.cex, FUN=max)))+1), main=x$config.names[m])
		axis(1)
		points(x$result.cex[m,1], par("usr")[3], pch=19, xpd=TRUE)

	}
	}
}
