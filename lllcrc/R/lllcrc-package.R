

#' Local Log-linear Models for Capture-Recapture
#' 
#' Applies local log-linear capture-recapture models (LLLMs) a for closed
#' populations, as described in the ongoing doctoral thesis work of Zachary
#' Kurtz in the department of statistics at Carnegie Mellon University.  The
#' method is relevant when there are 3-4 capture occasions, with auxiliary
#' covariates available for all capture occasions.  As part of estimating the
#' number of missing population units, the method estimates the "rate of
#' missingness" as it varies over the covariate space.  In addition,
#' user-friendly functions are provided to recreate (approximately) the method
#' of Zwane and van der Heijden (2004), which applied the \code{VGAM} package
#' in a way that is closely related to LLLMs.
#' 
#' \tabular{ll}{ Package: \tab lllcrc\cr Type: \tab Package\cr Version: \tab
#' 1.1\cr Date: \tab 2013-12-16\cr License: GPL-2
#' <http://cran.r-project.org/web/licenses/GPL-2> \cr }
#' 
#' @name lllcrc-package
#' @docType package
#' @author Zach Kurtz 
#' Maintainer: Zach \email{zkurtz@@gmail.com}
#' @references 
#' Kurtz ZT (2013). "Smooth Post-Stratification for Multiple
#' Capture-Recapture." \emph{arXiv preprint arXiv:1302.0890}.
#' @examples
#' 
#' ########################################################################
#' ## Simulate a dataset much like the French multiple sclerosis data
#' ##	(El Adssi et. al., 2012)
#' set.seed(100)
#' pop = poptop(k = 3, n=5000, heters = "french.1", covs = "age.sex.zip")
#' dt = formatdata(x = pop, y = "y", cont.vars = c("age"), categ.vars = c("sex", "zip"))
#' # Notice the new variable notation:
#' attributes(dt)$var.keys
#' 
#' 
#' ########################################################################
#' ## Implement basic log-linear models (LLM)
#' llmod = flat.IC(pop, rasch = FALSE, ic = "AICc", adjust = TRUE)
#' boot.list = list(n.reps = 25, est = llmod$pred, pop = pop, ic = "AICc", rasch = FALSE, 
#' 	adjust = TRUE, averaging = FALSE)
#' llmod.boots = llcrc.flat.boots(boot.list)
#' hist(llmod.boots); abline(v = llmod$pred, lwd = 3, col = "red")
#' 
#' ########################################################################
#' ## Use local log-linear models (LLLM)
#' # LLL Step 1: Use your excellent researcher intuition to set a bandwidth for each dimension:
#' x.covs = names(dt)[substr(names(dt), 1,1) == "x"]
#' bw.key = data.frame(bw = matrix(NA, nrow = length(x.covs), ncol = 1))
#' row.names(bw.key) = x.covs
#' bw.key[1,1] = 5
#' bw.key[2:7,1] = sqrt(2)
#' bw.key = as.matrix(bw.key)
#' # LLL Step 2: Obtain point estimates and bootstrap the variance for LLLMs:
#' #	 (For more-thorough results, increase the number n.reps of bootstrap replications)
#' lmod = lllcrc(dat = dt, kfrac = 0.2, bw = bw.key, ic = "AICc", round.vars = c("x.con.1"), 
#' 	rounding.scale = 20, boot.control = list(n.reps = 10, seed = 13))
#' summary(lmod)
#' # View rates of missingness by category
#' rates.by.category(x=lmod, reference.levels = c("x.dis.1.F","x.dis.2.a"))
#' # View a pointwise C.I. for the rate of missingness across age in a chosen category:
#' ci = extract.CI(mod = lmod, probs = c(0.025, 0.975), cont.var = "x.con.1", 
#' 	selection = c("x.dis.1.F", "x.dis.2.a"))
#' plot(ci$x, ci$pi0, ylim = c(0,max(ci$upper)), bty = "n", xlab = "age", type = "l", 
#' 	lwd = 2, ylab = "Rate of missingness", main = "LLLM")
#' lines(ci$x, ci$lower, lty = 2)
#' lines(ci$x, ci$upper, lty = 2)
#' mtext("females in zip code 57")
#' # Alternatively, view the relative capture pattern frequencies in a stacked form,
#' #	along with the confidence interval:
#' plot(lmod, selection = c("x.dis.1.F", "x.dis.2.a"), main = "LLLM")
#' 
#' \dontrun{
#' ########################################################################
#' ## Use a VGAM for CRC
#' # For better results, increase the number n.reps of bootstrap replications
#' vg = vgam.crc(dat = dt, sdf = 4, round.vars = c("x.con.1"), rounding.scale = 10, 
#' 	boot.control = list(n.reps=10, seed = 4))
#' summary(vg)
#' # View rates of missingness by category
#' rates.by.category(x=vg, reference.levels = c("x.dis.1.F","x.dis.2.a"))
#' # View a pointwise C.I. for the rate of missingness across age in a chosen category:
#' ci = extract.CI(mod = vg, probs = c(0.025, 0.975), cont.var = "x.con.1", 
#' 	selection = c("x.dis.1.F", "x.dis.2.a"))
#' plot(ci$x, ci$pi0, ylim = c(0,max(ci$upper)), bty = "n", xlab = "age", type = "l", 
#' 	lwd = 2, ylab = "Rate of missingness", main = "VGAM for CRC")
#' lines(ci$x, ci$lower, lty = 2)
#' lines(ci$x, ci$upper, lty = 2)
#' mtext("females in zip code 57")
#' # Alternatively, view the relative capture pattern frequencies in a stacked form,
#' #	along with the confidence interval:
#' plot(vg, selection = c("x.dis.1.F", "x.dis.2.a"), main = "VGAM for CRC")
#' }
NULL



