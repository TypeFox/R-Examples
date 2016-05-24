#' Study design in the presence of error-prone diagnostic tests and 
#' self-reported outcomes
#' 
#' This function calculates the power and sample in the presence of error-prone 
#' diagnostic tests and self-reported outcomes.  Two missing mechanisms can be 
#' assumed, namely MCAR and NTFP. The MCAR setting assumes that each test is 
#' subject to a constant, independent probability of missingness. The NTFP 
#' mechanism includes two types of missingness - (1) incorporates a constant, 
#' independent, probability of missing for each test prior to the first positive
#' test result; and (2) all test results after first positive are missing.
#' 
#' @param HR hazard ratio under the alternative hypothesis.
#' @param sensitivity the sensitivity of test.
#' @param specificity the specificity of test
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
#' @param pcensor a value or a vector (must have same length as testtimes) of
#' the probability of censoring at each visit, assuming censoring process
#' is independent on other missing mechanisms.
#' @param design missing mechanism: "MCAR" or "NTFP".
#' @param negpred baseline negative predictive value, i.e. the probability of 
#'   being truely disease free for those who were tested (reported) as disease 
#'   free at baseline. If baseline screening test is perfect, then negpred = 1.
#'   
#' @details To calculate sample sizes for a vector of powers, set N = NULL. To 
#'   calculate powers for a vector of sample sizes, set power = NULL. One and 
#'   only one of power and N should be specified, and the other set to NULL. 
#'   This function uses an enumeration algorithm to calculate the expected 
#'   Fisher information matrix. The expected Fisher information matrix is used 
#'   to obtain the variance of the coefficient corresponding to the treatment 
#'   group indicator.
#'   
#' @return \itemize{ \item result: a data frame with calculated sample size and 
#'   power \item I1 and I2: calculated unit Fisher information matrices for each
#'   group, which can be used to calculate more values of sample size and power 
#'   for the same design without the need to enumerate again }
#'   
#' @note   When diagnostic test is perfect, i.e. sensitivity=1 and 
#'   specificity=1, use \code{\link{icpowerpf}} instead to obtain significantly 
#'   improved computational efficiency.
#'   
#' @seealso \code{\link{icpowerpf}}
#'   
#' @examples
#' ## First specificy survivals. Assume test times are 1:8, with survival function 
#' ## at the end time 0.9  				  
#' surv <- exp(log(0.9)*(1:8)/8)					  
#'
#' ## Obtain power vs. N					  				
#' pow1 <- icpower(HR = 2, sensitivity = 0.55, specificity = 0.99, survivals = surv, 
#'    N = seq(500, 1500, 50), power = NULL, rho = 0.5, alpha = 0.05, 
#'    pmiss = 0, design = "MCAR", negpred = 1)
#' 
#' plot(pow1$result$N, pow1$result$power, type="l", xlab="N", ylab="power")
#' 
#' ## Calculate sample size, assuming desired power is 0.9
#' pow2 <- icpower(HR = 2, sensitivity = 0.55, specificity = 0.99, survivals = surv,
#'    N = NULL, power = 0.9, rho = 0.5, alpha = 0.05, pmiss = 0, design = "MCAR",
#'    negpred = 1)
#' 
#' ## When missing test is present with MCAR
#' pow3 <- icpower(HR = 2, sensitivity = 0.55, specificity = 0.99, survivals = surv,
#'    N = NULL, power = 0.9, rho = 0.5, alpha = 0.05, pmiss = 0.4, design = "MCAR",
#'    negpred = 1)
#' 
#' ## When missing test is present with NTFP
#' pow4 <- icpower(HR = 2, sensitivity = 0.55, specificity = 0.99, survivals = surv,
#'    N = NULL, power = 0.9, rho = 0.5, alpha = 0.05, pmiss = 0.4, design = "NTFP",
#'    negpred = 1)
#' 
#' ## When baseline misclassification is present
#' pow5 <- icpower(HR = 2, sensitivity = 0.55, specificity = 0.99, survivals = surv,
#'    N = NULL, power = 0.9, rho = 0.5, alpha = 0.05, pmiss = 0, design = "MCAR",
#'    negpred = 0.98)		 
#' 
#' ## When test is  perfect and no missing test		 
#' pow6 <- icpower(HR = 2, sensitivity = 1, specificity = 1, survivals = surv,
#'    N = NULL, power = 0.9, rho = 0.5, alpha = 0.05, pmiss = 0, design = "MCAR",
#'    negpred = 1)	
#' 
#' ## Different missing probabilities at each test time
#' pow7 <- icpower(HR = 2, sensitivity = 0.55, specificity = 0.99, survivals = surv,
#'    N = NULL, power = 0.9, rho = 0.5, alpha = 0.05, pmiss = seq(0.1, 0.8, 0.1),
#'    design = "MCAR")	  
#'    
#' @export

icpower <- function(HR, sensitivity, specificity, survivals, N = NULL, power = NULL,
                    rho = 0.5, alpha = 0.05, pmiss = 0, pcensor = 0, design = "MCAR", negpred = 1) {
  
  ## Basic input check
  if (!(HR>0)) stop("Check input for HR")
  if (!(specificity<=1 & specificity>0)) stop("Check input for specificity")
  if (sum(diff(survivals)>0)>0 | sum(survivals<0)>0 | sum(survivals>1)>0) stop("Check input for survivals")
  if (rho<0 | rho>1) stop("Check input for rho")
  if (alpha<0 | alpha>1) stop("Check input for alpha")
  if (any(pmiss < 0) | any(pmiss > 1)) stop("Check input for pmiss")
  if ((is.null(N)+is.null(power))!=1) stop("Need input for one and only one of power and N")
  if (!is.null(power)) {if (any(power > 1) | any(power < 0)) stop("Check input for power")}
  if (length(pmiss)!=1 & length(pmiss)!=length(survivals))
    stop("length of pmiss must be 1 or same length as survivals")
  if (!(design %in% c("MCAR", "NTFP"))) stop("invalid design")
  if (negpred<0 | negpred>1) stop("Check input for negpred")  	
  
  ## Use dedicated function for perfect test, so when perfect test use alternative instead
  if (sensitivity==1 & specificity==1) {
    cat("For perfect test, use icpowerpf for improved computational efficiency... \n \n")
    design <- "NTFP"
  }
  
  miss 	<- any(pmiss != 0) || any(pcensor != 0)
  beta  	<- log(HR)
  J 		<- length(survivals)
  Surv  	<- c(1, survivals)  	                   	
  if (length(pmiss)==1) pmiss <- rep(pmiss, J)
  if (length(pcensor) == 1) pcensor <- rep(pcensor, J)
  
  if (design=="MCAR" & !miss) {	
    Dm <- powerdmat1(sensitivity, specificity, J, negpred)
    prob <- 1	
  } else if (design=="MCAR" & miss) {
    Dm <- powerdmat2(sensitivity, specificity, J, negpred, pmiss, pcensor)
    prob <- Dm[[2]]
    Dm <- Dm[[1]]
  } else if (design=="NTFP" & !miss) {
    Dm <- powerdmat3(sensitivity, specificity, J, negpred)
    prob <- 1
  } else if (design=="NTFP" & miss) {
    Dm <- powerdmat4(sensitivity, specificity, J, negpred, pmiss, pcensor)
    prob <- Dm[[2]]
    Dm <- Dm[[1]]
  }
  
  I1 <- I2 <- matrix(NA, nrow=J+1, ncol=J+1)
  
  ## Group 1 information matrix
  DS         	<- c(prob/(Dm%*%Surv))
  I1[1, ] 	<- I1[, 1] 	<- 0
  I1[-1, -1] 	<- (t(DS*Dm)%*%Dm)[-1, -1]
  
  ## Group 2 information matrix	
  a 		<- HR
  Surv2 	<- Surv^a
  Surv21	<- Surv^(a-1)
  lSurv 	<- log(Surv)
  DS2		<- c(1/(Dm%*%Surv2))
  
  I2[1, 1] <- sum(prob*((Dm%*%(Surv2*lSurv))^2*a^2*DS2 - Dm%*%(Surv2*lSurv*(1+a*lSurv))*a))
  I2[-1, -1] <- (t(Dm*prob*DS2)%*%Dm*outer(Surv21, Surv21, "*")*a^2)[-1, -1] 	
  I2[cbind(2:(J+1), 2:(J+1))] <- colSums(prob*t(Surv21^2*t(Dm^2*DS2)*a^2 - t(Dm)*Surv^(a-2)*a*(a-1)))[-1]	
  I2[1, 2:(J+1)] <- I2[2:(J+1), 1] <- colSums(prob*t(t(c(Dm%*%(Surv2*lSurv))*DS2*Dm)*Surv21*a^2-t(Dm)*Surv21*(1+a*lSurv)*a))[-1]
  
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
