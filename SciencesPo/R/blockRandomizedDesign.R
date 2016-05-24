#' @encoding UTF-8
#' @title Create Block-randomized designs
#'
#' @description Generate block-randomized designs based on the number of units \code{n} and block size, where the block size is the number of experimental conditions. The number of Independent Variables and the number of levels in each IV are specified as input. The output is a the block randomized design. This function is intended for planning randomized trails.
#'
#' @param blocksize is the number of control blocks or n per block/group.
#' @param n is the total number of subjects or units.
#' @param seed the random number generation seed.
#'
#' @references
#' Alan S Gerber, Donald P Green (2012). \emph{Field experiments: Design, analysis, and interpretation}. WW Norton.
#'
#' RB Morton, KC Williams (2010). \emph{Experimental political science and the study of causality: From nature to the lab}. Cambridge University Press.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' blk <- blockRandomizedDesign(blocksize = 20, n = 80, seed = 51)
#' blk;
#' table(blk$block, blk$condition)
#'
#' # let's do some analysis
#' set.seed(51);
#' blk$y <- rnorm(n = 80, mean = 20, sd = 5)
#'
#' # Let's look at some descriptives:
#' tapply(blk$y, list(blk$condition, blk$block), mean)
#' tapply(blk$y, list(blk$condition, blk$block), sd)
#'
#' # Do the ANOVA and make some graphs
#' # This formula describes the response `y` by both the treatment factor `condition` and the block control `block`. Note that aov() treats `block` as a random error component of the variance, while lm() treats `block` as a fixed effect.
#'
#' fit.aov <- aov(y ~ factor(condition) + factor(block), data=blk)
#' summary(fit.aov) # display Type I ANOVA table
#' drop1(fit.aov,~.,test="F") # type III SS and F Tests
#'
#' # Since the p-value of 0.254 is much greater than the .05 significance level, we cannot reject the null hypothesis that the mean of `y` for each treatment conditions are all equal.
#'
#' model.tables(fit.aov, "means", se=TRUE) # SE for differences, NOT for means
#' # Calculate the pooled standard error of the means.
#' pooled.se = sqrt(1688.1/4)
#'
#' block <- c(1,2,3,4) # the values of the x axis
#' outcome <- c(19.76, 20.03, 18.44, 18.16) # the results from the means output
#' plot(block, outcome, type = "b", ylab = "outcome", xlab = "blocks of experimental conditions", ylim = c(0, 30) )
#'
#' fit.lm <- lm(y ~ factor(condition) + factor(block), data = blk)
#' anova(fit.aov)
#'
#'@export
#'
#' @importFrom stats runif
#'
`blockRandomizedDesign` = function(blocksize, n, seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  # blocking factor
  block = rep(1:ceiling(n/blocksize), each = blocksize)
  a1 = data.frame(id= 1: length(block), block, rand=runif(length(block)))
  a2 = a1[order(a1$block,a1$rand),]
  # matching treatment
  a2$condition = rep(c("Treat", "Control"),times = length(block)/2)
  assign = a2[order(a2$id),]
  class(assign) <- c("SciencesPo", "randomize", "data.frame")
  return(assign)
}
NULL
