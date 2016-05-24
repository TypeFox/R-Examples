#' Study design in the presence of error-prone diagnostic tests and 
#' self-reported outcomes when sensitivity and specificity are unkonwn and a 
#' validation set is used
#' 
#' This function calculates the power and sample size in the presence of 
#' error-prone diagnostic tests and self-reported outcomes when both sensitivity
#' and specificity are unknown. In this case, a subject of the subjects receive 
#' both gold standard test and error-prone test at each non-missing visit. The 
#' remaining subjects receive only error-prone test. Here, for the validation 
#' set, NTFP refers to no test after first positive result from the gold 
#' standard test. Both sensitivity and specificity are treated as unknown 
#' parameters in this setting.
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
#' @param rhoval proportion of subjects in validation set.
#' @param rho proportion of subjects in baseline(reference) group.
#' @param alpha type I error.
#' @param pmiss a value or a vector (must have same length as survivals) of the 
#'   probabilities of each test being randomly missing at each test time. If 
#'   pmiss is a single value, then each test is assumed to have an identical 
#'   probability of missingness.
#' @param design missing mechanism: "MCAR" or "NTFP".
#' @param designval missing mechanism of validation set: "MCAR" or "NTFP".
#' @param negpred baseline negative predictive value, i.e. the probability of 
#'   being truely disease free for those who were tested (reported) as disease 
#'   free at baseline. If baseline screening test is perfect, then negpred = 1.
#'   
#' @return \itemize{ \item result: a data frame with calculated sample size and
#' power \item IR1 and IR2: calculated unit Fisher information matrices for each
#' group in non-validation set \item IV1 and IV2: calculated unit Fisher
#' information matrices for each group in validation set }
#' 
#' @examples
#' surv <- exp(log(0.9)*(1:8)/8)
#' pow <- icpower.val(HR = 2, sensitivity = 0.55, specificity = 0.99, 
#'    survivals = surv, power = 0.9, rhoval=0.05, design= "NTFP", designval = "NTFP")
#' pow$result
#' 
#' @export

icpower.val <- function(HR, sensitivity, specificity, survivals, N = NULL, power = NULL, rhoval, 
                   rho = 0.5, alpha = 0.05, pmiss = 0, design = "MCAR", designval= "MCAR", negpred = 1) {

   	## Basic input check
   	if (!(HR>0)) stop("Check input for HR")
   	if (!(sensitivity<=1 & sensitivity>0)) stop("Check input for sensitivity")
   	if (!(specificity<=1 & specificity>0)) stop("Check input for specificity")
   	if (sum(diff(survivals)>0)>0 | sum(survivals<0)>0 | sum(survivals>1)>0) stop("Check input for survivals")
   	if (rho<0 | rho>1) stop("Check input for rho")
   	if (alpha<0 | alpha>1) stop("Check input for alpha")
   	if (sum(pmiss<0)>0 | sum(pmiss>1)>0) stop("Check input for pmiss")
   	if ((is.null(N)+is.null(power))!=1) stop("Need input for one and only one of power and N")
   	if (!is.null(power)) {if (any(power > 1) | any(power < 0)) stop("Check input for power")}
   	if (length(pmiss)!=1 & length(pmiss)!=length(survivals))
       stop("length of pmiss must be 1 or same length as survivals")
  	if (!(design %in% c("MCAR", "NTFP"))) stop("invalid design")
  	if (negpred<0 | negpred>1) stop("Check input for negpred")

	beta <- log(HR)

	RegFisher <- function(HR, phi1, phi0, survivals, pmiss = 0, design = "MCAR", negpred = 1) {
		miss 	<- any(pmiss!=0)  
   		beta  	<- log(HR)
   		J 		<- length(survivals)
   		Surv  	<- c(1, survivals)  	                   	
		if (length(pmiss)==1) pmiss <- rep(pmiss, J)
		if (design=="MCAR" & !miss) {	
			## MCAR design with no missing tests
			data 	<- as.matrix(do.call("expand.grid", rep(list(0:1), J))) 
			prob 	<- 1	
		} else if (design=="MCAR" & miss) {
			## MCAR design with missing tests
			data 	<- as.matrix(do.call("expand.grid", rep(list(0:2), J)))[-3^J, ]
			prob 	<- apply(pmiss^(t(data==2))*(1-pmiss)^(t(data!=2)), 2, prod)
		} else if (design=="NTFP" & !miss) {
			## NTFP without missing tests
			data <- sapply(0:(J-1), function(x) c(rep(0, x), 1, rep(2, J-x-1)))
			data <- t(cbind(data, rep(0, J)))
			prob	<- 1
		} else if (design=="NTFP" & miss) {
			## NTFP with missing tests
			crdata <- function(n) {
				if (n==J+1) {
					data <- as.matrix(do.call("expand.grid", rep(list(c(0, 2)), J)))[-2^J, ]
					probs <-  apply(pmiss^(t(data==2))*(1-pmiss)^(t(data!=2)), 2, prod) 
				} else {
					data <- as.matrix(do.call("expand.grid", c(rep(list(c(0, 2)), n-1), 1, rep(list(2), J-n))))
					p <- pmiss[1:n]
					probs <-  apply(p^(t(data[, 1:n]==2))*(1-p)^(t(data[, 1:n]!=2)), 2, prod) 	
				}
				return(cbind(data, probs))
			}
			data <- do.call("rbind", sapply(1:(J+1), crdata))
			prob <-  data[, J+1]
			data <- data[, 1:J] 
		}

		o <- outer(1:J, 1:(J+1), ">=") 
		mat1 <- apply(o, 2, function(x) rowSums(data==1 & matrix(x, nrow=nrow(data), ncol=J, byrow=T)))
		mat2 <- apply(o, 2, function(x) rowSums(data==0 & matrix(x, nrow=nrow(data), ncol=J, byrow=T)))
		mat3 <- apply(o, 2, function(x) rowSums(data==1 & !matrix(x, nrow=nrow(data), ncol=J, byrow=T)))
		mat4 <- apply(o, 2, function(x) rowSums(data==0 & !matrix(x, nrow=nrow(data), ncol=J, byrow=T)))
		Cm <- phi1^mat1*(1-phi1)^mat2*(1-phi0)^mat3*phi0^mat4
		## Create transformation matrix and obtain D matrix
   		mat     <- c(1, rep(0, J), rep(c(-1, 1, rep(0, J)), J-1), -1, 1)
   		S2theta <- matrix(mat, ncol=J+1)
   		C2C     <- negpred*diag(J+1) + (1-negpred)*matrix(c(1, rep(0,J)),
   		       nrow=J+1, ncol=J+1)   ## Adjust for baseline misclassification  
		Tm <- C2C%*%S2theta
		Dm <- Cm%*%Tm
		## Get Derivative matrices
		Dm1 <- ((ifelse(mat1>0, mat1*phi1^(mat1-1), 0)*(1-phi1)^mat2 - 
         	ifelse(mat2>0, mat2*(1-phi1)^(mat2-1), 0)*phi1^mat1)*(1-phi0)^mat3*phi0^mat4)%*%Tm
		Dm0 <- ((ifelse(mat4>0, mat4*phi0^(mat4-1), 0)*(1-phi0)^mat3 - 
         	ifelse(mat3>0, mat3*(1-phi0)^(mat3-1), 0)*phi0^mat4)*(1-phi1)^mat2*phi1^mat1)%*%Tm
		Dm11 <-  ((ifelse(mat1>1, mat1*(mat1-1)*phi1^(mat1-2), 0)*(1-phi1)^mat2 +
         	ifelse(mat2>1, mat2*(mat2-1)*(1-phi1)^(mat2-2), 0)*phi1^mat1 -
         	ifelse(mat1 > 0 & mat2 > 0, 2*mat1*mat2*phi1^(mat1-1)*(1-phi1)^(mat2-1), 0))*
         	(1-phi0)^mat3*phi0^mat4)%*%Tm
		Dm00 <- ((ifelse(mat4>1, mat4*(mat4-1)*phi0^(mat4-2), 0)*(1-phi0)^mat3 +
         	ifelse(mat3>1, mat3*(mat3-1)*(1-phi0)^(mat3-2), 0)*phi0^mat4 -
         	ifelse(mat4 > 0 & mat3 > 0, 2*mat4*mat3*phi0^(mat4-1)*(1-phi0)^(mat3-1), 0))*
         	(1-phi1)^mat2*phi1^mat1)%*%Tm
		Dm01 <- ((ifelse(mat1>0, mat1*phi1^(mat1-1), 0)*(1-phi1)^mat2 - 
         	ifelse(mat2>0, mat2*(1-phi1)^(mat2-1), 0)*phi1^mat1)*
        	(ifelse(mat4>0, mat4*phi0^(mat4-1), 0)*(1-phi0)^mat3 - 
         	ifelse(mat3>0, mat3*(1-phi0)^(mat3-1), 0)*phi0^mat4))%*%Tm
         
		I1 <- I2 <- matrix(NA, nrow=J+3, ncol=J+3)
		## Group 1 information
		DS1         	<- c(prob/(Dm%*%Surv))
		I1[1, ] 	<- I1[, 1] 	<- 0
		I1[2:(J+1), 2:(J+1)] 	<- (t(DS1*Dm)%*%Dm)[-1, -1]
		I1[J+2, J+2] <- colSums((Dm1%*%Surv)^2*DS1 - Dm11%*%Surv*prob)
		I1[J+3, J+3] <- colSums((Dm0%*%Surv)^2*DS1 - Dm00%*%Surv*prob)
		I1[J+2, J+3] <- I1[J+3, J+2] <- colSums((Dm1%*%Surv)*(Dm0%*%Surv)*DS1 - Dm01%*%Surv*prob)
		I1[J+2, 2:(J+1)] <- I1[2:(J+1), J+2] <- colSums(c((Dm1%*%Surv)*DS1)*Dm - Dm1*prob)[-1]
		I1[J+3, 2:(J+1)] <- I1[2:(J+1), J+3] <- colSums(c((Dm0%*%Surv)*DS1)*Dm - Dm0*prob)[-1]                
		## Group 2 information
		a 		<- HR
		Surv2 	<- Surv^a
		Surv21	<- Surv^(a-1)
		lSurv 	<- log(Surv)
		DS2		<- c(1/(Dm%*%Surv2))
		I2[1, 1] <- sum(prob*((Dm%*%(Surv2*lSurv))^2*a^2*DS2 - Dm%*%(Surv2*lSurv*(1+a*lSurv))*a))
		I2[2:(J+1), 2:(J+1)] <- (t(Dm*prob*DS2)%*%Dm*outer(Surv21, Surv21, "*")*a^2)[-1, -1]
		I2[cbind(2:(J+1), 2:(J+1))] <- colSums(prob*t(Surv21^2*t(Dm^2*DS2)*a^2 - t(Dm)*Surv^(a-2)*a*(a-1)))[-1]	
		I2[1, 2:(J+1)] <- I2[2:(J+1), 1] <- colSums(prob*t(t(c(Dm%*%(Surv2*lSurv))*DS2*Dm)*Surv21*a^2-t(Dm)*Surv21*(1+a*lSurv)*a))[-1]
		I2[J+2, J+2] <- colSums(prob*((Dm1%*%Surv2)^2*DS2 - Dm11%*%Surv2))
		I2[J+3, J+3] <- colSums(prob*((Dm0%*%Surv2)^2*DS2 - Dm00%*%Surv2))
		I2[J+2, J+3] <- I2[J+3, J+2] <- colSums(prob*((Dm1%*%Surv2)*(Dm0%*%Surv2)*DS2 - Dm01%*%Surv2))
		I2[J+2, 2:(J+1)] <- I2[2:(J+1), J+2] <- colSums(prob*(t(t(c((Dm1%*%Surv2)*DS2)*Dm - Dm1)*a*Surv21)))[-1]
		I2[J+3, 2:(J+1)] <- I2[2:(J+1), J+3] <- colSums(prob*(t(t(c((Dm0%*%Surv2)*DS2)*Dm - Dm0)*a*Surv21)))[-1]
		I2[1, J+2] <- I2[J+2, 1] <- colSums(prob*((c((Dm1%*%Surv2)*DS2)*Dm - Dm1)%*%(Surv2*lSurv*a)))
		I2[1, J+3] <- I2[J+3, 1] <- colSums(prob*((c((Dm0%*%Surv2)*DS2)*Dm - Dm0)%*%(Surv2*lSurv*a)))                  
		return(list(I1, I2))	
	}	
	
	VSetFisher <- function(HR, phi1, phi0, survivals, pmiss = 0, design = "MCAR") {
		miss 	<- any(pmiss!=0)  
   		beta  	<- log(HR)
   		J 		<- length(survivals)
   		Surv  	<- c(1, survivals)  	                   	
		if (length(pmiss)==1) pmiss <- rep(pmiss, J)
		sensitivity <- specificity <- 1	
		if (design=="MCAR" & !miss) {	
			## MCAR design with no missing tests
        	Cm <- diag(rep(1, J+1))
			prob 	<- 1
			n1 <- J:0
			n0 <- 0:J
		} else if (design=="MCAR" & miss) {
			## MCAR design with missing tests
			datagen <- function(n) {
				if (n==(J+1)) {
					return(as.matrix(do.call("expand.grid", rep(list(c(0, 2)), J)))[-2^J, ])
				} else {
					return(as.matrix(do.call("expand.grid", c(rep(list(c(0, 2)), n-1), 1, rep(list(1:2), J-n)))))
				}
			}
			data <- do.call("rbind", sapply(1:(J+1), datagen))
			o 		<- 3*outer(1:J, 1:(J+1), ">=") + 1 
			smat 	<- c(specificity, 1-specificity, 1, 1-sensitivity, sensitivity, 1)
			Cm	 	<- sapply(1:J, function(x) smat[outer(data[, x], o[x, ], "+")])
			Cm	 	<- apply(Cm, 1, prod)
			Cm	 	<- matrix(Cm, ncol=J+1)
			prob 	<- apply(pmiss^(t(data==2))*(1-pmiss)^(t(data!=2)), 2, prod)
			n1 <- rowSums(data==1)
			n0 <- rowSums(data==0)		
		} else if (design=="NTFP" & !miss) {
			## NTFP without missing tests
			Cm <- diag(rep(1, J+1))
			prob <- 1
			n1 <- c(rep(1, J), 0)
			n0 <- 0:J	
		} else if (design=="NTFP" & miss) {
			## NTFP with missing tests
			smat <- c(specificity, 1-specificity, 1, 1-sensitivity,  sensitivity, 1)
			crmat <- function(n) {
				if (n==0) {
					data <- matrix(1)		  
				} else if (n==J) {
					data <- as.matrix(do.call("expand.grid", rep(list(c(0, 2)), n)))[-2^J, ]
				} else {
					data <- cbind(as.matrix(do.call("expand.grid", rep(list(c(0, 2)), n))), 1)
				}
				n1 	<- min(n+1, J)
				o 	<- 3*outer(1:n1, 1:(J+1), ">=") + 1
				mats <- sapply(1:n1, function(x) smat[outer(data[, x], o[x, ], "+")])
				mats <- apply(mats, 1, prod)
				mats <- matrix(mats, ncol=J+1)
				p 	<- pmiss[1:n1]
				probs <-  apply(p^(t(data==2))*(1-p)^(t(data!=2)), 2, prod) 
				n1 <- rowSums(data==1)
				n0 <- rowSums(data==0)	
				mats <- cbind(mats, probs, n1, n0)
				return(mats)
			}
			Cm <- do.call("rbind", sapply(0:J, crmat))
			prob	<- Cm[, J+2]
			n1 <- Cm[, J+3]
			n0 <- Cm[, J+4]
			Cm		<- Cm[, 1:(J+1)]
		}
		## Create transformation matrix and obtain D matrix
   		mat     <- c(1, rep(0, J), rep(c(-1, 1, rep(0, J)), J-1), -1, 1)
   		S2theta <- matrix(mat, ncol=J+1)  
   		Dm      <- Cm%*%S2theta     ### no need to adjust baseline misclassification
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
   	 	IV1 <- IV2 <- matrix(0, nrow=J+3, ncol=J+3)	
    	IV1[1:(J+1), 1:(J+1)] <- I1
    	IV1[J+2, J+2] <- sum(c(Dm%*%Surv)*prob*n1/(phi1*(1-phi1)))
    	IV1[J+3, J+3] <- sum(c(Dm%*%Surv)*prob*n0/(phi0*(1-phi0)))
    	IV2[1:(J+1), 1:(J+1)] <- I2
    	IV2[J+2, J+2] <- sum(c(Dm%*%Surv2)*prob*n1/(phi1*(1-phi1)))
    	IV2[J+3, J+3] <- sum(c(Dm%*%Surv2)*prob*n0/(phi0*(1-phi0)))
    	return(list(IV1, IV2))
	}	
	
	IR <- RegFisher(HR, sensitivity, specificity, survivals, pmiss, design, negpred)
	IR1 <- IR[[1]]
	IR2 <- IR[[2]]

	IV <- VSetFisher(HR, sensitivity, specificity, survivals, pmiss, designval) 
	IV1 <- IV[[1]]
	IV2 <- IV[[2]]

	If       <- rhoval*rho*IV1 + rhoval*(1-rho)*IV2 + (1-rhoval)*rho*IR1 + (1-rhoval)*(1-rho)*IR2
   	inv.If   <- solve(If)
   	beta.var <- inv.If[1]

   	## Calculate Sample size or Power
   	if (!is.null(power)) {
    	N  <- ceiling((qnorm(1-alpha/2)+qnorm(power))^2*beta.var/beta^2)
        N1 <- round(rho*N)
        N2 <- N-N1
        return(list(result=data.frame(N,N1,N2,power), IR1=IR1, IR2=IR2, IV1=IV1, IV2=IV2))
   	} else {
        beta.varN <- beta.var/N
        za    <- qnorm(1-alpha/2)
        zb    <- beta/sqrt(beta.varN)
        power <- pnorm(-za,zb,1)+1-pnorm(za,zb,1)        
        N1 <- round(rho*N)
        N2 <- N-N1        
        return(list(result=data.frame(N,N1,N2,power), IR1=IR1, IR2=IR2, IV1=IV1, IV2=IV2))
  	}
}
