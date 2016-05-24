#' False discovery rate estimation for massively parallel restricted likelihood
#' ratio tests
#' 
#' Given a set of RLRT results and a threshold, this function outputs an
#' estimate of the FDR (in the empirical Bayes sense of Efron, 2010) when the
#' given threshold is used to determine which null hypotheses to reject.
#' 
#' 
#' @param rlrt.obj an RLRT object obtained from \code{\link{rlrt.mp}} or
#' \code{\link{rlrt4d}}.
#' @param threshold threshold at which the null hypothesis is rejected.
#' @return A list with elements \item{MoM}{FDR based on method of moments
#' estimator of RLRT parameters (Greven et al., 2008).} \item{ML}{FDR based on
#' maximum likelihood estimation of RLRT parameters, as described in Greven et
#' al. (2008).}
#' @author Philip Reiss \email{phil.reiss@@nyumc.org}
#' @seealso \code{\link{rlrt.mp}}, \code{\link{rlrt4d}}
#' @references Efron, B. (2010).  \emph{Large-Scale Inference: Empirical Bayes
#' Methods for Estimation, Testing, and Prediction}.  New York: Cambridge
#' University Press.
#' 
#' Greven, S., Crainiceanu, C. M., Kuechenhoff, H., and Peters, A. (2008).
#' Restricted likelihood ratio testing for zero variance components in linear
#' mixed models.  \emph{Journal of Computational and Graphical Statistics},
#' 17(4), 870--891.
#' @examples
#' 
#' # See example for rlrt.mp
#' @export
Fdr.rlrt <- function(rlrt.obj, threshold) {
	rstat = pmax(0, rlrt.obj$stat)
	rsim = rlrt.obj$sim
	tbar = mean(rsim)
	p.mom = 1 - 3 * tbar^2 / mean(rsim^2)
	a.mom = tbar / (1 - p.mom)
	pi0 = mean(rstat == 0) / p.mom
	cat("Estimated null proportion", pi0, "\n")
	Fdr1 = pi0 * (1 - p.mom) * (1 - pchisq(threshold / a.mom, 1)) / mean(rstat >= threshold)
	Fdr2 = pi0 * mean(rsim >= threshold) / mean(rstat >= threshold)
	otpt = c(Fdr1, Fdr2)
	names(otpt) = c("MoM", "ML")
	otpt
}

