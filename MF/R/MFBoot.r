#' Estimates bootstrap confidence intervals for the mitigated fraction.
#' 
#' Resamples the data and produces bootstrap confidence intervals. Equal tailed intervals are estimated by the 
#' percentile method. Highest density intervals are estimated by selecting the shortest of all possible intervals. 
#' For BCa intervals, see Efron and Tibshirani section 14.3.
#' 
#' @title Bootstrap MF CI
#' @usage MFBoot(formula, data, compare = c("con", "vac"), b = 100, B = 100, 
#'    alpha = 0.05, hpd = TRUE, bca = FALSE, return.boot = FALSE, trace.it = FALSE)
#' @param formula Formula of the form \code{y ~ x}, where y is a continuous response and x is a factor with two levels
#' @param data Data frame
#' @param compare Text vector stating the factor levels - \code{compare[1]} is the control or reference group to which \code{compare[2]} is compared
#' @param b Number of bootstrap samples to take with each cycle
#' @param B Number of cycles, giving the total number of samples = B * b
#' @param alpha Complement of the confidence level
#' @param hpd Estimate highest density intervals? Default TRUE.
#' @param bca Estimate BCa intervals? Default FALSE.
#' @param return.boot Save the bootstrap sample of the MF statistic? Default FALSE.
#' @param trace.it Verbose tracking of the cycles? Default FALSE.
#' @return a \code{\link{mfboot-class}} data object
#' @seealso \code{\link{mfboot-class}} 
#' @export
#' @references Siev D. (2005). An estimator of intervention effect on disease severity. \emph{Journal of Modern Applied Statistical Methods.} \bold{4:500--508} \cr \cr
#' Efron B, Tibshirani RJ. \emph{An Introduction to the Bootstrap.} Chapman and Hall, New York, 1993.
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @examples 
#' MFBoot(lesion~group, calflung)
#'
#' # 10000 bootstrap samples
#' # 95% confidence interval
#' #
#' # Comparing vac to con 
#' #                 observed median  lower  upper
#' # Equal Tailed        0.44 0.4464 0.1360 0.7056
#' # Highest Density     0.44 0.4464 0.1456 0.7088
##
##--------------------------------------------------------------------
## Bootstrap simple MF
##--------------------------------------------------------------------
MFBoot <- function(formula, data, compare = c("con", "vac"), b = 100, B = 100, 
	alpha = 0.05, hpd = TRUE, bca = FALSE, return.boot = FALSE, trace.it = FALSE){
    # bootstrap confidence intervals for MF
    # 11/17/99 initial coding
    # 2/24/04 added BC.a interval
    # 5/25/10 added empirical HPD interval
    # takes b bootstrap samples B times, so nboot = B * b

   

    A <- data.frame(model.frame(formula = formula, data = data))
    resp <- A[, 1]
    tx <- A[, 2]
    x <- resp[tx == compare[1]]
    y <- resp[tx == compare[2]]

    rng <- 'Mersenne-Twister'
    RNGkind(rng)
    seed <- .Random.seed
    nboot <- b * B
    n.x <- length(x)
    n.y <- length(y)
    w <- function(xy, n.x){sum(rank(xy)[1:n.x])}
    W <- rep(NA, b * B)
    for(i in 1:B) {
        x.b <- matrix(sample(x, size = b * n.x, replace = T), b, n.x)
        y.b <- matrix(sample(y, size = b * n.y, replace = T), b, n.y)
        W[(i - 1) * b + (1:b)] <- apply(cbind(x.b, y.b), 1, w, n.x)
        if(trace.it){ cat("bootstrap samples", (i - 1) * b + 1, "to", (i - 1) * b + b, "\n")}
    }
    MF <- ((2 * W - n.x * (1 + n.x + n.y))/(n.x * n.y))
    mf <- ((2 * w(c(x, y), n.x) - n.x * (1 + n.x + n.y))/(n.x * n.y))
    qprob <- c(0.5, alpha/2, 1 - alpha/2)
    qmf <- quantile(MF, prob = qprob)
    stat <- matrix(c(mf, qmf), 1, 4, dimnames = list(c('Equal Tailed'), c("observed", "median", "lower", "upper")))
    if(hpd){
        hpdmf <- emp.hpd(MF, alpha=alpha)
        stat <- rbind(stat, 'Highest Density' = c(mf, median(MF), hpdmf))
    }
    if(bca) {
        z0 <- qnorm(sum(MF < mf)/(b * B))
        xy <- c(x, y)
        nx.i <- rep(n.x - (1:0), c(n.x, n.y))
        ny.i <- rep(n.y - (0:1), c(n.x, n.y))
        theta <- rep(NA, n.x + n.y)
        for(i in 1:(n.x + n.y)){
			theta[i] <- ((2 * w(xy[ - i], nx.i[i]) - nx.i[i] * (1 + nx.i[i] + ny.i[i]))/(nx.i[i] * ny.i[i]))
		}
        theta.hat <- mean(theta)
        acc <- sum((theta.hat - theta)^3)/(6 * sum((theta.hat - theta)^2)^(3/2))
        z1 <- qnorm(alpha/2)
        z2 <- qnorm(1 - alpha/2)
        a1 <- pnorm(z0 + (z0 + z1)/(1 - acc * (z0 + z1)))
        a2 <- pnorm(z0 + (z0 + z2)/(1 - acc * (z0 + z2)))
        a5 <- pnorm(z0 + (z0 + 0)/(1 - acc * (z0 + 0)))
        qprob <- c(a5, a1, a2)
        stuff <- c(acc, z0, a1, a2)
        names(stuff) <- c("acc", "z0", "a1", "a2")
        if(trace.it){
			print(round(stuff, 4))
		}
        qmf <- quantile(MF, prob = qprob)
        stat <- rbind(stat, 'BC.a'= c(mf, qmf))
    }
    if(return.boot){
        out <- mfboot$new(stat = stat, nboot = nboot, alpha = alpha, seed = seed, 
			rng = rng, compare = compare, sample = MF)
	} else {
		out <- mfboot$new(stat = stat, nboot = nboot, alpha = alpha, seed = seed, 
			rng = rng, compare = compare, sample = NULL)
	}

    return(out)
}
