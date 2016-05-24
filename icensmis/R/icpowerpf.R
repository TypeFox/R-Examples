#' Study design in the presence of interval censored outcomes (assuming perfect
#' diagnostic tests)
#' 
#' This function implements power and sample size calculations for interval
#' censored time-to-event outcomes, when the diagnostic tests are assumed to be
#' perfect (i.e. sensitivity=1 and specificity=1). This is a special case of the
#' more general study design function \code{\link{icpower}}. However, for the
#' special case of perfect diagnostic tests, this function can be used with
#' significantly improved computational efficiency.
#' 
#' @param HR hazard ratio under the alternative hypothesis.
#' @param survivals a vector of survival function at each test time for 
#'   baseline(reference) group. Its length determines the number of tests.
#' @param N a vector of sample sizes to calculate corresponding powers. If one 
#'   needs to calculate sample size, then set to NULL.
#' @param power a vector of powers to calculate corresponding sample sizes. If 
#'   one needs to calculate power, then set to NULL.
#' @param rho proportion of subjects in baseline(reference) group.
#' @param alpha type I error.
#' @param pmiss a value or a vector (must have same length as survivals) of the 
#'   probabilities of each test being randomly missing at each test time. If 
#'   pmiss is a single value, then each test is assumed to have an identical 
#'   probability of missingness.
#'   
#' @note See \code{\link{icpower}} for more details in a general situation.
#' 
#' @return same form as returned value of \code{\link{icpower}}
#' 
#' @examples
#' powpf1 <- icpowerpf(HR =2 , survivals = seq(0.9, 0.1, by=-0.1), N = NULL,
#'    power = 0.9, pmiss = 0)
#'    
#' powpf2 <- icpowerpf(HR =2 , survivals = seq(0.9, 0.1, by=-0.1), N = NULL,
#'    power = 0.9, pmiss = 0.7)
#'    
#' ## Different missing probabilities at each test time
#' powpf3 <- icpowerpf(HR =2 , survivals = seq(0.9, 0.1, -0.1), N = NULL, 
#'    power = 0.9, pmiss = seq(0.1, .9, 0.1))    
#'    
#' @export

icpowerpf <- function(HR, survivals, N = NULL, power = NULL, rho = 0.5, alpha = 0.05, pmiss = 0) {

   	## Basic input check
   	if (!(HR>0)) stop("Check input for HR")
   	if (sum(diff(survivals)>0)>0 | sum(survivals<0)>0 | sum(survivals>1)>0) stop("Check input for survivals")
   	if (rho<0 | rho>1) stop("Check input for rho")
   	if (alpha<0 | alpha>1) stop("Check input for alpha")
   	if (sum(pmiss<0)>0 | sum(pmiss>1)>0) stop("Check input for pmiss")
   	if (length(pmiss)!=1 & length(pmiss)!=length(survivals))
       stop("length of pmiss must be 1 or same length as survivals")
   	if ((is.null(N)+is.null(power))!=1) stop("Need input for one and only one of power and N")
   	if (!is.null(power)) {if (any(power > 1) | any(power < 0)) stop("Check input for power")}
   	
	miss <- ifelse(sum(pmiss!=0)>0, TRUE, FALSE)     	         	   	
	beta <- log(HR)
	surv <- c(1, survivals, 0)
	J <- length(survivals)
	if (length(pmiss)==1) {pmiss <- rep(pmiss, J)}
	
	## Enumerate all possible results
	L <- rep(1:(J+2), each=J+2)
	R <- rep(1:(J+2), times=J+2)
	keep <- L<R & !(L==1 & R==J+2)
	L <- L[keep]
	R <- R[keep]
	## calculate probability of each missing pattern
	pmissE <- c(0, pmiss, 0)
	compphi <- function(x, y) {
		if (y==x+1) {
			return((1-pmissE[x])*(1-pmissE[y]))
		} else {
			return ((1-pmissE[x])*(1-pmissE[y])*prod(pmissE[(x+1):(y-1)]))
		}
	}
	phi <- sapply(1:length(L), function(i) compphi(L[i], R[i]))		

	## Group 1 Fisher Information (updated 10-18)
	survL <- surv[L]
	survR <- surv[R]
	survDif <- 1/(survL-survR)*phi
	I1 <- I2 <- matrix(NA, nrow=J+1, ncol=J+1)
	I1[1, ] <- I1[, 1] <- 0
	I1[2:(J+1), 2:(J+1)] <- sapply(2:(J+1), function(x) sum(survDif[L==x])+sum(survDif[R==x]))	
	i <- (L!=1 & L!=J+2 &R!=1 &R!=J+2 & L!=R)
	I1[cbind(L[i], R[i])] <- I1[cbind(R[i], L[i])] <- -survDif[i]

	## Group 2 Fisher Information (updated 10-18)
	a <- HR
	sLexp <- survL^(a)
	sLexp1 <- survL^(a-1)
	sLexp2 <- survL^(a-2)
	sLlog <- log(survL)
	sLexplog <- sLexp*sLlog
	sRexp <- survR^(a)
	sRexp1 <- survR^(a-1)
	sRexp2 <- survR^(a-2)
	sRlog <- ifelse(R==J+2, 0, log(survR))
	sRexplog <- sRexp*sRlog
	survD <- 1/(sLexp-sRexp)
	## for beta
	I2[1, 1] <- sum(phi*(survD*(sLexplog-sRexplog)^2*a^2 - sLexplog*(1+sLlog*a)*a + sRexplog*(1+sRlog*a)*a))
	## for Sj and beta
	ssL <- phi*(survD*(sLexp1*a)^2 - a*(a-1)*sLexp2)
	ssR <- phi*(survD*(sRexp1*a)^2 + a*(a-1)*sRexp2)
	I2[2:(J+1), 2:(J+1)] <- sapply(2:(J+1), function(x) sum(ssL[L==x])+sum(ssR[R==x]))
	sbL <- phi*(survD*(sLexplog-sRexplog)*sLexp1*a^2 - sLexp1*a^2*(1+sLlog))
	sbR <- phi*(-survD*(sLexplog-sRexplog)*sRexp1*a^2 + sRexp1*a^2*(1+sRlog))
	I2[1, 2:(J+1)] <- I2[2:(J+1), 1] <- sapply(2:(J+1), function(x) sum(sbL[L==x])+sum(sbR[R==x]))
	## for Sj and Sj'
	I2[cbind(L[i], R[i])] <- I2[cbind(R[i], L[i])] <- -phi[i]*survD[i]*a^2*sLexp1[i]*sRexp1[i]


   	If       <- rho*I1+(1-rho)*I2
   	inv.If   <- solve(If)
   	beta.var <- inv.If[1]

   	## Calculate Sample size or Power
   	if (!is.null(power)) {
    	N  <- ceiling((qnorm(1-alpha/2)+qnorm(power))^2*beta.var/beta^2)
        N1 <- round(rho*N)
        N2 <- N-N1
        return(list(result=data.frame(N,N1,N2,power), I1=I1, I2=I2))
        } else {
        beta.varN <- beta.var/N
        za    <- qnorm(1-alpha/2)
        zb    <- beta/sqrt(beta.varN)
        power <- pnorm(-za,zb,1)+1-pnorm(za,zb,1)        
        N1 <- round(rho*N)
        N2 <- N-N1        
        return(list(result=data.frame(N,N1,N2,power), I1=I1, I2=I2))
        }
}
