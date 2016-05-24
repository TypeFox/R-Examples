`summary.hglm` <-
	function(object, ...) {

object$nRand <- cumsum(object$RandC)
Call <- object$call
Method <- object$method
if (object$Converge == "did not converge") stop("There is no valid estimate to produce summary statistics")
if (is.null(object$fixef)) stop("There in no valid estimate to produce summary statistics")
if (is.null(object$SeFe)) stop("There in no valid standard error estimate to produce summary statistics")

### Fixed effects summary ###
Nfix <- length(object$fixef)
FixCoefMat <- matrix(numeric(4*Nfix), nrow = Nfix, byrow = TRUE)
dimnames(FixCoefMat) <- list(names(object$fixef), c("Estimate", "Std. Error", "t-value", "Pr(>|t|)"))
FixCoefMat[,1] <- as.numeric(object$fixef)
FixCoefMat[,2] <- as.numeric(object$SeFe)
FixCoefMat[,3] <- as.numeric(FixCoefMat[,1]/FixCoefMat[,2])
FixCoefMat[,4] <- 2*pt(abs(as.numeric(FixCoefMat[,3])), object$dfReFe, lower.tail = FALSE)

### Random effects summary ###

if (!is.null(object$ranef)) {
	if (length(object$RandC) == 1) {
		Nran <- object$RandC
		RandCoefMat <- matrix(numeric(2*Nran), nrow = Nran, byrow = TRUE)
		dimnames(RandCoefMat) <- list(names(object$ranef), c("Estimate", "Std. Error"))
		RandCoefMat[,1] <- as.numeric(object$ranef)
		RandCoefMat[,2] <- as.numeric(object$SeRe)
	} else {
		RandCoefMat <- c()
		for (J in length(object$RandC):2) {
			Nran <- object$RandC[J]
			RandCoefMat[[J]] <- matrix(numeric(2*Nran), nrow = Nran, byrow = TRUE)
			dimnames(RandCoefMat[[J]]) <- list(names(object$ranef[(object$nRand[J - 1] + 1):object$nRand[J]]), c("Estimate", "Std. Error"))
			RandCoefMat[[J]][,1] <- as.numeric(object$ranef[(object$nRand[J - 1] + 1):object$nRand[J]])
			RandCoefMat[[J]][,2] <- as.numeric(object$SeRe[(object$nRand[J - 1] + 1):object$nRand[J]])
		}
		Nran <- object$RandC[1]
		RandCoefMat[[1]] <- matrix(numeric(2*Nran), nrow = Nran, byrow = TRUE)
		dimnames(RandCoefMat[[1]]) <- list(names(object$ranef[1:object$nRand[1]]), c("Estimate", "Std. Error"))
		RandCoefMat[[1]][,1] <- as.numeric(object$ranef[1:object$nRand[1]])
		RandCoefMat[[1]][,2] <- as.numeric(object$SeRe[1:object$nRand[1]])
	}
} else RandCoefMat <- NULL

smst <- list(Method = Method, FixCoefMat = FixCoefMat, RandCoefMat = RandCoefMat, RandC = object$RandC,
             nRand = object$nRand, SummVC1 = object$SummVC1, SummVC2 = object$SummVC2, iter = object$iter,
             converge = object$Converge, call = Call, ProfLogLik = object$ProfLogLik, devdf = object$dfReFe,
             LogLik = object$LogLik, varFix = object$varFix, varRanef = object$varRanef, link.disp = object$link.disp,
			 likelihood = object$likelihood, call.rand.family = object$call.rand.family, bad = object$bad,
			 CAR.tau = object$CAR.tau, CAR.rho = object$CAR.rho, SAR.tau = object$SAR.tau, SAR.rho = object$SAR.rho)
  
class(smst) <- "summary.hglm"
  
return(smst)

}

