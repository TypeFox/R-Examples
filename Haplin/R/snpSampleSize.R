snpSampleSize <- function(fam.cases, fam.controls, fraction = 0.5, RR, MAF, alpha = 0.05, power = 0.80){
##
##
## Sample size calculations for a single SNP. Allows for sample size calculations of several combinations at once.
##
## fam.cases is a character vector of the case family design. The possible family designs are "mfc" (full triad), "mc" (mother-child dyad), "fc" (father-child dyad) and "c" (a single case child).
## fam.controls is a character vector of the control family design. The possible family designs are "mfc" (full triad), "mc" (mother-child-dyad), "fc" (father-child dyad), "mf" (mother-father dyad), "c" (a single control child), "m" (a single control mother), "f" (a single control father) or "no_controls" (no control families). 
## fraction is a numeric vector of the proportion of case families. Equals 0.5 by default, i.e. there are as many case families as control families. If fam.controls equals "no_controls", fraction is automatically set to 1.
## RR is a numeric vector of the relative risks.
## MAF is a numeric vector of the minor allele frequencies.
## alpha is a numeric vector of the Type I Errors. Equals 0.05 by default.
## power is a numeric vector of the desired probability of identifying a difference in the relative risks. Default is 0.80.
##
## All arguments must be of equal length or have length equal to one.
##
##
#
## Function for calculating the power, based on normal approximation on log(OR)
or.power <- function(p2, OR, n1, fraction, alpha = 0.05, power = 0.80){
	n2 <- n1*(1-fraction)/fraction
	p1 <- p2*OR/(1-p2+p2*OR)
	.sd <- sqrt((1/(n1*p1*(1-p1)))+1/(n2*p2*(1-p2)))
	return(1+pnorm(-qnorm(1-alpha/2)-log(OR)/.sd)-pnorm(qnorm(1-alpha/2)-log(OR)/.sd) - power)
}
#
## Function for computing the upper limit for uniroot()
upper.limit <- function(p2, OR, fraction, alpha = 0.05, power = 0.80){
	.power <- power
	z.alpha <- qnorm(1-alpha/2)
   	z.beta <- qnorm(.power)
	p1 <- p2*OR/(1-p2+p2*OR)
	k <- fraction/(1-fraction)
	n1 <- (z.beta+z.alpha)^2/(log(OR))^2*(1/(p1*(1-p1))+k/(p2*(1-p2)))
}
#
.power <- power
#
## Possible family designs
.case.des <- c("mfc", "mc", "fc", "c")
.control.des <- c("mfc", "mc", "fc", "mf", "c", "m", "f", "no_controls")
#
## Checking if all arguments have the same length or have length equal to one
.check.length <- unique(c(length(fam.cases), length(fam.controls), length(fraction), length(RR), length(MAF), length(alpha), length(.power)))
if(length(.check.length[.check.length != 1]) >= 2) stop("All arguments must be of equal length or have length equal to 1", call. = F)
#
## Misc errors
if(!all(fam.cases %in% .case.des)) stop("Argument \"fam.cases\" is not specified correctly", call. = F)
if(!all(fam.controls %in% .control.des)) stop("Argument \"fam.controls\" is not specified correctly", call. = F)
if(!is.numeric(fraction) || any(fraction <= 0 | fraction > 1)) stop("Invalid fraction value(s)", call. = F)
if(!is.numeric(RR) || any(RR <= 0)) stop("The relative risk(s) must be numeric and larger than 0", call. = F)
if(!is.numeric(MAF) || any(MAF <= 0 | MAF >= 1)) stop("\"MAF\" contains invalid value(s)", call. = F)
if(!is.numeric(alpha) || any(alpha <= 0 | alpha >= 1)) stop("\"alpha\" contains invalid Type I Error(s)", call. = F)
if(!is.numeric(.power) || any(.power <= 0 | .power >= 1) | any(.power <= alpha)) stop("\"power\" is not specified correctly", call. = F) 
#
## Weights for the different family designs
.case.w <- c(2, 2, 2, 2)
.pseudo.w <- c(2, 1, 1, 0)
.control.w <- c(4, 3, 3, 4, 2, 2, 2, 0)
names(.case.w) <- .case.des
names(.pseudo.w) <- .case.des
names(.control.w) <- .control.des
#
.allele.weights <- as.data.frame(cbind(cases = .case.w[fam.cases], controls = .control.w[fam.controls], pseudo = .pseudo.w[fam.cases], fraction = fraction))
.allele.weights$fraction[which(fam.controls == "no_controls")] <- 1
rownames(.allele.weights) <- NULL
#
## Finding the fraction of case alleles
.case.alleles <- .allele.weights$cases*.allele.weights$fraction
.control.alleles <- .allele.weights$controls*(1-fraction)+.allele.weights$pseudo*.allele.weights$fraction
.allele.fraction <- .case.alleles/(.control.alleles+.case.alleles)
#
if(any(.allele.fraction == 1)) stop("There are no control alleles present", call. = F)
#
.tab <- cbind(RR = RR, MAF = MAF, alpha = alpha, power = .power, allele.fraction = .allele.fraction)
#
## Calculating the upper limit of the interval in uniroot()
.upper.limit <- apply(.tab, 1, function(x){
 	ceiling(upper.limit(p2 = x["MAF"], OR = x["RR"], fraction = x["allele.fraction"], alpha = x["alpha"], power = x["power"]))
})
#
.tab <- cbind(.tab, upper.limit = .upper.limit)
#
## Sample size calculations
#
## Normal approximation on log(OR)
.or.sample.size <- apply(.tab, 1, function(x){
 	uniroot(f = or.power, interval = c(0,x["upper.limit"]), p2 = x["MAF"], OR = x["RR"], fraction = x["allele.fraction"], alpha = x["alpha"], power = x["power"], tol = 0.5)$root
})
# 
## Sample sizes
.n1 <- .or.sample.size/.allele.weights$cases
.n2 <- .n1*(1-.allele.weights$fraction)/.allele.weights$fraction
#
## Return output
.output <- as.data.frame(cbind(fam.cases, fam.controls, .tab, case.families = ceiling(.n1), control.families = ceiling(.n2)))
.output[, c("allele.fraction", "p2", "upper.limit")] <- list(NULL, NULL, NULL)
rownames(.output) <- NULL
return(.output)
}
