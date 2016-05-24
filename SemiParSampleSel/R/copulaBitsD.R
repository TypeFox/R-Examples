copulaBitsD <- function (VC, theta.star, F1, F2, F22, f2, precision) {

 if(VC$BivD=="FGM") {
    theta <- tanh(theta.star)
    FGM1 <- F1*F2*(1+theta*(1-F1)*(1-F2))
    FGM2 <- F1*F22*(1+theta*(1-F1)*(1-F22))
    
    FGM1 <- ifelse(FGM1>precision, FGM1, precision)
    FGM1 <- ifelse(FGM1<(1-precision), FGM1, 1-precision)
    FGM2 <- ifelse(FGM2>precision, FGM2, precision)
    FGM2 <- ifelse(FGM2<(1-precision), FGM2, 1-precision)
    lx <- f2 - FGM1 + FGM2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- F2*((1+theta*(1-F1)*(1-F2))-F1*theta*(1-F2))
    dcop11 <- F22*((1+theta*(1-F1)*(1-F22))-F1*theta*(1-F22))
    
    dcop2 <- F1*((1+theta*(1-F2)*(1-F1))-F2*theta*(1-F1))
    dcop22 <- F1*((1+theta*(1-F22)*(1-F1))-F22*theta*(1-F1))
    
    dcop.theta1 <- (1-F1)*(1-F2)*F1*F2
    dcop.theta2 <- (1-F1)*(1-F22)*F1*F22
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- F2*(-2*theta*(1-F2))
    d2Cdcop112 <- F22*(-2*theta*(1-F22))
    
    d2Cdcop22 <- F1*(-2*theta*(1-F1))
    d2Cdcop222 <- F1*(-2*theta*(1-F1))
    
    d2Cdcop1dcop2 <- (1+theta*(1-F1)*(1-F2)-F1*theta*(1-F2))+F2*(-theta*(1-F1)+F1*theta)
    d2Cdcop11dcop22 <- (1+theta*(1-F1)*(1-F22)-F1*theta*(1-F22))+F22*(-theta*(1-F1)+F1*theta)
    
    d2Cdcop.theta12 <- F2*((1-F1)*(1-F2)-F1*(1-F2))
    d2Cdcop.theta22 <- F22*((1-F1)*(1-F22)-F1*(1-F22))
    
    d2Cdcop1dcop.theta1 <- F1*((1-F2)*(1-F1)-F2*(1-F1))
    d2Cdcop11dcop.theta2 <- F1*((1-F22)*(1-F1)-F22*(1-F1))
    
    d2Cdcop2dcop.theta1 <- 0
    d2Cdcop22dcop.theta2 <- 0
    
   
    theta.append<-1-theta^2
    theta.append.der<- -2*theta*(1-theta^2)
  
    
  } else if (VC$BivD=="N") {
    theta <- tanh(theta.star)
    BN1 <- pbinorm(qnorm(F1), qnorm(F2), mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = theta) 
    BN2 <- pbinorm(qnorm(F1), qnorm(F22), mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = theta)
    
    BN1 <- ifelse(BN1>precision, BN1, precision)
    BN1 <- ifelse(BN1<(1-precision), BN1, 1-precision)
    BN2 <- ifelse(BN2>precision, BN2, precision)
    BN2 <- ifelse(BN2<(1-precision), BN2, 1-precision)
    
    lx <- f2 - BN1 + BN2
    lx <- ifelse(lx>precision, lx, precision)
    
    
    dcop1 <- pnorm((qnorm(F2)-theta*qnorm(F1))/sqrt(1-theta^2))
    dcop11 <- pnorm((qnorm(F22)-theta*qnorm(F1))/sqrt(1-theta^2))
    
    dcop2 <- pnorm((qnorm(F1)-theta*qnorm(F2))/sqrt(1-theta^2))
    dcop22 <- pnorm((qnorm(F1)-theta*qnorm(F22))/sqrt(1-theta^2))
    
    dcop.theta1 <- dbinorm(qnorm(F1), qnorm(F2), mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = theta, log = FALSE)
    dcop.theta2 <- dbinorm(qnorm(F1), qnorm(F22), mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = theta, log = FALSE)
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- dnorm((qnorm(F2)-theta*qnorm(F1))/sqrt(1-theta^2))*(-theta*(dnorm(qnorm(F1)))^(-1)*(1-theta^2)^(-1/2))
    d2Cdcop112 <- dnorm((qnorm(F22)-theta*qnorm(F1))/sqrt(1-theta^2))*(-theta*(dnorm(qnorm(F1)))^(-1)*(1-theta^2)^(-1/2))
    
    d2Cdcop22 <- dnorm((qnorm(F1)-theta*qnorm(F2))/sqrt(1-theta^2))*(-theta*(dnorm(qnorm(F2)))^(-1)*(1-theta^2)^(-1/2))
    d2Cdcop222 <- dnorm((qnorm(F1)-theta*qnorm(F22))/sqrt(1-theta^2))*(-theta*(dnorm(qnorm(F22)))^(-1)*(1-theta^2)^(-1/2))
    
    d2Cdcop1dcop2 <- dnorm((qnorm(F2)-theta*qnorm(F1))/sqrt(1-theta^2))*(dnorm(qnorm(F2)))^(-1)*(1-theta^2)^(-1/2)
    d2Cdcop11dcop22 <- dnorm((qnorm(F22)-theta*qnorm(F1))/sqrt(1-theta^2))*(dnorm(qnorm(F22)))^(-1)*(1-theta^2)^(-1/2)
    
    d2Cdcop.theta12 <- dnorm((qnorm(F2)-theta*qnorm(F1))/sqrt(1-theta^2))*(((-qnorm(F1))*(1-theta^2)^(1/2)+theta*(qnorm(F2)-theta*qnorm(F1))*(1-theta^2)^(-1/2))/(1-theta^2))
    d2Cdcop.theta22 <- dnorm((qnorm(F22)-theta*qnorm(F1))/sqrt(1-theta^2))*(((-qnorm(F1))*(1-theta^2)^(1/2)+theta*(qnorm(F22)-theta*qnorm(F1))*(1-theta^2)^(-1/2))/(1-theta^2))
    
    d2Cdcop1dcop.theta1 <- dnorm((qnorm(F1)-theta*qnorm(F2))/sqrt(1-theta^2))*(((-qnorm(F2))*(1-theta^2)^(1/2)+theta*(qnorm(F1)-theta*qnorm(F2))*(1-theta^2)^(-1/2))/(1-theta^2))
    d2Cdcop11dcop.theta2 <- dnorm((qnorm(F1)-theta*qnorm(F22))/sqrt(1-theta^2))*(((-qnorm(F22))*(1-theta^2)^(1/2)+theta*(qnorm(F1)-theta*qnorm(F22))*(1-theta^2)^(-1/2))/(1-theta^2))
    
    d2Cdcop2dcop.theta1  <- (1-theta^2)^(-1)*theta*dcop.theta1+dcop.theta1*(-(-2*qnorm(F1)*qnorm(F2)*(1-theta^2)+(qnorm(F1)^2-2*theta*qnorm(F1)*qnorm(F2)+qnorm(F2)^2)*(2*theta))/(2*(1-theta^2)^2))
    d2Cdcop22dcop.theta2  <- (1-theta^2)^(-1)*theta*dcop.theta2+dcop.theta2*(-(-2*qnorm(F1)*qnorm(F22)*(1-theta^2)+(qnorm(F1)^2-2*theta*qnorm(F1)*qnorm(F22)+qnorm(F22)^2)*(2*theta))/(2*(1-theta^2)^2))
    
    
    theta.append<-1-theta^2
    theta.append.der<- -2*theta*(1-theta^2)
  
    
    
  } else if (VC$BivD=="AMH") {
    theta <- tanh(theta.star)
    AMH1 <- (F1*F2)/(1-theta*(1-F1)*(1-F2))
    AMH2 <- (F1*F22)/(1-theta*(1-F1)*(1-F22))
    AMH1 <- ifelse(AMH1>precision, AMH1, precision)
    AMH1 <- ifelse(AMH1<(1-precision), AMH1, 1-precision)
    AMH2 <- ifelse(AMH2>precision, AMH2, precision)
    AMH2 <- ifelse(AMH2<(1-precision), AMH2, 1-precision)
    lx <- f2 - AMH1 + AMH2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- (F2*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F2)))/(1-theta*(1-F1)*(1-F2))^2
    dcop11 <- (F22*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F22)))/(1-theta*(1-F1)*(1-F22))^2
    
    dcop2 <- (F1*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F1)))/(1-theta*(1-F1)*(1-F2))^2
    dcop22 <- (F1*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F1)))/(1-theta*(1-F1)*(1-F22))^2
    
    dcop.theta1 <- ((F1*F2)/(1-theta*(1-F1)*(1-F2))^2)*(1-F1)*(1-F2)
    dcop.theta2 <- ((F1*F22)/(1-theta*(1-F1)*(1-F22))^2)*(1-F1)*(1-F22)
    
    
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- -((2 * (-1 + F2) * F2 * theta * (1 + (-1 + F2) * theta))/(-1 + (-1 + F1) * (-1 + F2) * theta)^3)
    d2Cdcop112 <- -((2 * (-1 + F22) * F22 * theta * (1 + (-1 + F22) * theta))/(-1 + (-1 + F1) * (-1 + F22) * theta)^3)
    
    d2Cdcop22 <- -((2*(-1 + F1)*F1*theta*(1 + (-1 + F1)*theta))/(-1 + (-1 + F1)*(-1 + F2)*theta)^3)
    d2Cdcop222 <- -((2*(-1 + F1)*F1*theta*(1 + (-1 + F1)*theta))/(-1 + (-1 + F1)*(-1 + F22)*theta)^3)
    
    d2Cdcop1dcop2 <- (-1 - (-2 + F1 + F2 + F1 * F2) * theta + (-1 + F1 + F2 - F1 * F2) * theta^2)/(-1 + (-1 + F1) * (-1 + F2) * theta)^3
    d2Cdcop11dcop22 <- (-1 - (-2 + F1 + F22 + F1 * F22) * theta + (-1 + F1 + F22 - F1 * F22) * theta^2)/(-1 + (-1 + F1) * (-1 + F22) * theta)^3
    
    d2Cdcop.theta12 <- -(((-1 + F2) * F2 * (-1 + theta - F2 * theta + F1 * (2 + (-1 + F2) * theta)))/(-1 + (-1 + F1) * (-1 + F2) * theta)^3)
    d2Cdcop.theta22 <- -(((-1 + F22) * F22 * (-1 + theta - F22 * theta + F1 * (2 + (-1 + F22) * theta)))/(-1 + (-1 + F1) * (-1 + F22) * theta)^3)
    
    d2Cdcop1dcop.theta1 <- -(((-1 + F1) * F1 * (-1 + theta - F1 * theta + F2 * (2 + (-1 + F1) * theta)))/(-1 + (-1 + F1) * (-1 + F2) * theta)^3)
    d2Cdcop11dcop.theta2 <- -(((-1 + F1) * F1 * (-1 + theta - F1 * theta + F22 * (2 + (-1 + F1) * theta)))/(-1 + (-1 + F1) * (-1 + F22) * theta)^3)
    
    d2Cdcop2dcop.theta1 <- -(2 * (-1 + F1)^2 * F1 * (-1 + F2)^2 * F2)/(-1 + (-1 + F1) * (-1 + F2) * theta)^3
    d2Cdcop22dcop.theta2 <- -(2 * (-1 + F1)^2 * F1 * (-1 + F22)^2 * F22)/(-1 + (-1 + F1) * (-1 + F22) * theta)^3
    
    
    theta.append <- 1 - theta^2
    theta.append.der <- -2 * theta * (1 - theta^2)
  
    
    
    
  } else if (VC$BivD %in% c("C0", "C90", "C180", "C270")) {
    
    if (VC$BivD == "C0") {
      F1.base <- F1
      F2.base <- F2
      F22.base <- F22
    } else if (VC$BivD == "C90") {
      F1.base <- 1 - F1
      F2.base <- F2
      F22.base <- F22
    } else if (VC$BivD == "C180") {
      F1.base <- 1 - F1
      F2.base <- 1 - F2
      F22.base <- 1 - F22
    } else if (VC$BivD == "C270") {
      F1.base <- F1
      F2.base <- 1 - F2
      F22.base <- 1 - F22
    } 
    
    
    theta <- exp(theta.star) + precision
    Clayton1 <- (F1.base^(-theta)+F2.base^(-theta)-1)^(-1/theta)
    Clayton2 <- (F1.base^(-theta)+F22.base^(-theta)-1)^(-1/theta)
    Clayton1 <- ifelse(Clayton1>precision, Clayton1, precision)
    Clayton1 <- ifelse(Clayton1<(1-precision), Clayton1, 1-precision)
    Clayton2 <- ifelse(Clayton2>precision, Clayton2, precision)
    Clayton2 <- ifelse(Clayton2<(1-precision), Clayton2, 1-precision)
    lx <- f2 - Clayton1 + Clayton2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- F1.base^(-theta-1)*(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-1)
    dcop11 <- F1.base^(-theta-1)*(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-1)
    
    dcop2 <- F2.base^(-theta-1)*(F2.base^(-theta)+F1.base^(-theta)-1)^(-(1/theta)-1)
    dcop22 <- F22.base^(-theta-1)*(F22.base^(-theta)+F1.base^(-theta)-1)^(-(1/theta)-1)
    
    dcop.theta1 <- ((-1+F1.base^(-theta)+F2.base^(-theta))^(-1/theta)*((theta*(F2.base^theta*log(F1.base)+F1.base^theta*log(F2.base)))/(F2.base^(theta)-F1.base^(theta)*(-1+F2.base^theta))+log(-1+F1.base^(-theta)+F2.base^(-theta))))/theta^2
    dcop.theta2 <- ((-1+F1.base^(-theta)+F22.base^(-theta))^(-1/theta)*((theta*(F22.base^theta*log(F1.base)+F1.base^theta*log(F22.base)))/(F22.base^(theta)-F1.base^(theta)*(-1+F22.base^theta))+log(-1+F1.base^(-theta)+F22.base^(-theta))))/theta^2
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- (F1.base^(-2+theta) * F2.base^theta*(-1 + F1.base^(-theta)+F2.base^-theta)^(-1/theta)*(-1 + F2.base^theta)*(1 + theta))/(F2.base^theta-F1.base^theta*(-1 + F2.base^theta))^2
    d2Cdcop112 <- (F1.base^(-2+theta) * F22.base^theta*(-1 + F1.base^(-theta)+F22.base^-theta)^(-1/theta)*(-1 + F22.base^theta)*(1 + theta))/(F22.base^theta-F1.base^theta*(-1 + F22.base^theta))^2
    
    d2Cdcop22 <- (-theta-1)*F2.base^(-theta-2)*(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-1)+F2.base^(-theta-1)*(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-2)*F2.base^(-theta-1)*(-theta)*(-1/theta-1)
    d2Cdcop222 <- (-theta-1)*F22.base^(-theta-2)*(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-1)+F22.base^(-theta-1)*(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-2)*F22.base^(-theta-1)*(-theta)*(-1/theta-1)
    
    d2Cdcop1dcop2 <- F1.base^(-theta-1)*(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*F2.base^(-theta-1)
    d2Cdcop11dcop22 <- F1.base^(-theta-1)*(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*F22.base^(-theta-1)
    
    d2Cdcop.theta12 <- F1.base^(-theta-1)*log(F1.base)*(-1)*(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-1)+F1.base^(-theta-1)*(dcop.theta1*(F1.base^(-theta)+F2.base^(-theta)-1)^(-1)+(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta))*(F1.base^(-theta)+F2.base^(-theta)-1)^(-2)*(F1.base^(-theta)*log(F1.base)+F2.base^(-theta)*log(F2.base)))
    d2Cdcop.theta22 <- F1.base^(-theta-1)*log(F1.base)*(-1)*(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-1)+F1.base^(-theta-1)*(dcop.theta2*(F1.base^(-theta)+F22.base^(-theta)-1)^(-1)+(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta))*(F1.base^(-theta)+F22.base^(-theta)-1)^(-2)*(F1.base^(-theta)*log(F1.base)+F22.base^(-theta)*log(F22.base)))
    
    d2Cdcop1dcop.theta1 <- F2.base^(-theta-1)*log(F2.base)*(-1)*(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-1)+F2.base^(-theta-1)*(dcop.theta1*(F1.base^(-theta)+F2.base^(-theta)-1)^(-1)+(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta))*(F1.base^(-theta)+F2.base^(-theta)-1)^(-2)*(F1.base^(-theta)*log(F1.base)+F2.base^(-theta)*log(F2.base)))
    d2Cdcop11dcop.theta2 <- F22.base^(-theta-1)*log(F22.base)*(-1)*(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-1)+F22.base^(-theta-1)*(dcop.theta2*(F1.base^(-theta)+F22.base^(-theta)-1)^(-1)+(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta))*(F1.base^(-theta)+F22.base^(-theta)-1)^(-2)*(F1.base^(-theta)*log(F1.base)+F22.base^(-theta)*log(F22.base)))
    
    d2Cdcop2dcop.theta1 <- -2/theta^3*log(F1.base^(-theta)+F2.base^(-theta)-1)*(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta))+1/theta^2*((F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-1)*(-F1.base^(-theta)*log(F1.base)-F2.base^(-theta)*log(F2.base))+log(F1.base^(-theta)+F2.base^(-theta)-1)*dcop.theta1)+(-1/theta^2*(F1.base^(-theta)*log(F1.base)+F2.base^(-theta)*log(F2.base))+1/theta*(-F1.base^(-theta)*(log(F1.base))^2-F2.base^(-theta)*(log(F2.base))^2))*(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-1)+1/theta*(F1.base^(-theta)*log(F1.base)+F2.base^(-theta)*log(F2.base))*(dcop.theta1*(F1.base^(-theta)+F2.base^(-theta)-1)^(-1)+(F1.base^(-theta)+F2.base^(-theta)-1)^(-(1/theta)-2)*(F1.base^(-theta)*log(F1.base)+F2.base^(-theta)*log(F2.base)))
    d2Cdcop22dcop.theta2 <- -2/theta^3*log(F1.base^(-theta)+F22.base^(-theta)-1)*(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta))+1/theta^2*((F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-1)*(-F1.base^(-theta)*log(F1.base)-F22.base^(-theta)*log(F22.base))+log(F1.base^(-theta)+F22.base^(-theta)-1)*dcop.theta2)+(-1/theta^2*(F1.base^(-theta)*log(F1.base)+F22.base^(-theta)*log(F22.base))+1/theta*(-F1.base^(-theta)*(log(F1.base))^2-F22.base^(-theta)*(log(F22.base))^2))*(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-1)+1/theta*(F1.base^(-theta)*log(F1.base)+F22.base^(-theta)*log(F22.base))*(dcop.theta2*(F1.base^(-theta)+F22.base^(-theta)-1)^(-1)+(F1.base^(-theta)+F22.base^(-theta)-1)^(-(1/theta)-2)*(F1.base^(-theta)*log(F1.base)+F22.base^(-theta)*log(F22.base)))
    
    
    
    theta.append <- theta-precision
    theta.append.der<- theta-precision
  

    
  } else if (VC$BivD=="F") {
    theta <- theta.star + precision
    
    Frank1 <- -theta^(-1)*log(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))
    Frank2 <- -theta^(-1)*log(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))
    
    Frank1 <- ifelse(Frank1>precision, Frank1, precision)
    Frank1 <- ifelse(Frank1<(1-precision), Frank1, 1-precision)
    Frank2 <- ifelse(Frank2>precision, Frank2, precision)
    Frank2 <- ifelse(Frank2<(1-precision), Frank2, 1-precision)
    lx <- f2 - Frank1 + Frank2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- -(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1)))*((exp(-theta*F2)-1)/(exp(-theta)-1))*exp(-theta*F1)*(-theta)
    dcop11 <- -(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1)))*((exp(-theta*F22)-1)/(exp(-theta)-1))*exp(-theta*F1)*(-theta) 
    
    dcop2 <- -(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1)))*((exp(-theta*F1)-1)/(exp(-theta)-1))*exp(-theta*F2)*(-theta)
    dcop22 <- -(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1)))*((exp(-theta*F1)-1)/(exp(-theta)-1))*exp(-theta*F22)*(-theta)
    
    dcop.theta1 <- (1/theta^2)*log(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1)))*((-exp(-theta*F1)*F1*(exp(-theta*F2)-1)-exp(-theta*F2)*F2*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)*exp(-theta))/(exp(-theta)-1)^2
    dcop.theta2 <- (1/theta^2)*log(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1)))*((-exp(-theta*F1)*F1*(exp(-theta*F22)-1)-exp(-theta*F22)*F22*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)*exp(-theta))/(exp(-theta)-1)^2
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- ((exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*(exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop112 <- ((exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*(exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop22 <- ((exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop222 <- ((exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop1dcop2 <- (exp(-theta*F1)/(exp(-theta)-1)*exp(-theta*F2)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop11dcop22 <- (exp(-theta*F1)/(exp(-theta)-1)*exp(-theta*F22)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop.theta12 <- (((exp(-theta)-1)*((exp(-theta*F2)-1)*exp(-theta*F1)*(-F1)+exp(-theta*F1)*exp(-theta*F2)*(-F2))+(exp(-theta*F2)-1)*exp(-theta*F1)*exp(-theta))/(exp(-theta)-1)^2*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*((exp(-theta)-1)*((exp(-theta*F2)-1)*exp(-theta*F1)*(-F1)+(exp(-theta*F1)-1)*exp(-theta*F2)*(-F2))+(exp(-theta*F2)-1)*(exp(-theta*F1)-1)*exp(-theta))/(exp(-theta)-1)^2 )/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop.theta22 <- (((exp(-theta)-1)*((exp(-theta*F22)-1)*exp(-theta*F1)*(-F1)+exp(-theta*F1)*exp(-theta*F22)*(-F22))+(exp(-theta*F22)-1)*exp(-theta*F1)*exp(-theta))/(exp(-theta)-1)^2*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*((exp(-theta)-1)*((exp(-theta*F22)-1)*exp(-theta*F1)*(-F1)+(exp(-theta*F1)-1)*exp(-theta*F22)*(-F22))+(exp(-theta*F22)-1)*(exp(-theta*F1)-1)*exp(-theta))/(exp(-theta)-1)^2 )/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop1dcop.theta1 <- (((exp(-theta)-1)*((exp(-theta*F1)-1)*exp(-theta*F2)*(-F2)+exp(-theta*F2)*exp(-theta*F1)*(-F1))+(exp(-theta*F1)-1)*exp(-theta*F2)*exp(-theta))/(exp(-theta)-1)^2*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*((exp(-theta)-1)*((exp(-theta*F2)-1)*exp(-theta*F1)*(-F1)+(exp(-theta*F1)-1)*exp(-theta*F2)*(-F2))+(exp(-theta*F2)-1)*(exp(-theta*F1)-1)*exp(-theta))/(exp(-theta)-1)^2 )/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop11dcop.theta2 <- (((exp(-theta)-1)*((exp(-theta*F1)-1)*exp(-theta*F22)*(-F22)+exp(-theta*F22)*exp(-theta*F1)*(-F1))+(exp(-theta*F1)-1)*exp(-theta*F22)*exp(-theta))/(exp(-theta)-1)^2*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*((exp(-theta)-1)*((exp(-theta*F22)-1)*exp(-theta*F1)*(-F1)+(exp(-theta*F1)-1)*exp(-theta*F22)*(-F22))+(exp(-theta*F22)-1)*(exp(-theta*F1)-1)*exp(-theta))/(exp(-theta)-1)^2 )/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop2dcop.theta1 <- -2/theta^3*log(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))+2/theta^2*1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))*((-exp(-theta*F1)*F1*(exp(-theta*F2)-1)-exp(-theta*F2)*F2*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)*exp(-theta))/(exp(-theta)-1)^2-1/theta*(-1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2*((exp(-theta)-1)*(exp(-theta*F1)*(-F1)*(exp(-theta*F2)-1)+(exp(-theta*F1)-1)*exp(-theta*F2)*(-F2))+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)*exp(-theta) )^2/(exp(-theta)-1)^4+1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))*( (exp(-theta)-1)^2*((exp(-theta*F1)*F1^2*(exp(-theta*F2)-1)+exp(-theta*F1)*F1*exp(-theta*F2)*F2+exp(-theta*F2)*F2^2*(exp(-theta*F1)-1)+exp(-theta*F2)*F2*exp(-theta*F1)*F1)*(exp(-theta)-1)-(-exp(-theta*F1)*F1*(exp(-theta*F2)-1)-exp(-theta*F2)*F2*(exp(-theta*F1)-1))*exp(-theta)+exp(-theta*F1-theta*F2-theta)*(-F1-F2-1)-exp(-theta*F1-theta)*(-F1-1)-exp(-theta*F2-theta)*(-F2-1)-exp(-theta))+2*(exp(-theta)-1)*exp(-theta)*((-exp(-theta*F1)*F1*(exp(-theta*F2)-1)-exp(-theta*F2)*F2*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)*exp(-theta) ) )/(exp(-theta)-1)^4         )  
    d2Cdcop22dcop.theta2 <- -2/theta^3*log(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))+2/theta^2*1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))*((-exp(-theta*F1)*F1*(exp(-theta*F22)-1)-exp(-theta*F22)*F22*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)*exp(-theta))/(exp(-theta)-1)^2-1/theta*(-1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2*((exp(-theta)-1)*(exp(-theta*F1)*(-F1)*(exp(-theta*F22)-1)+(exp(-theta*F1)-1)*exp(-theta*F22)*(-F22))+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)*exp(-theta) )^2/(exp(-theta)-1)^4+1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))*( (exp(-theta)-1)^2*((exp(-theta*F1)*F1^2*(exp(-theta*F22)-1)+exp(-theta*F1)*F1*exp(-theta*F22)*F22+exp(-theta*F22)*F22^2*(exp(-theta*F1)-1)+exp(-theta*F22)*F22*exp(-theta*F1)*F1)*(exp(-theta)-1)-(-exp(-theta*F1)*F1*(exp(-theta*F22)-1)-exp(-theta*F22)*F22*(exp(-theta*F1)-1))*exp(-theta)+exp(-theta*F1-theta*F22-theta)*(-F1-F22-1)-exp(-theta*F1-theta)*(-F1-1)-exp(-theta*F22-theta)*(-F22-1)-exp(-theta))+2*(exp(-theta)-1)*exp(-theta)*((-exp(-theta*F1)*F1*(exp(-theta*F22)-1)-exp(-theta*F22)*F22*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)*exp(-theta) ) )/(exp(-theta)-1)^4    )  
    
    
    theta.append<-1
    theta.append.der<- 0
  
    
  } else if (VC$BivD %in% c("G0", "G90", "G180", "G270")) {
    
    if (VC$BivD == "G0") {
      F1.base <- F1
      F2.base <- F2
      F22.base <- F22
    } else if (VC$BivD == "G90") {
      F1.base <- 1 - F1
      F2.base <- F2
      F22.base <- F22
    } else if (VC$BivD == "G180") {
      F1.base <- 1 - F1
      F2.base <- 1 - F2
      F22.base <- 1 - F22
    } else if (VC$BivD == "G270") {
      F1.base <- F1
      F2.base <- 1 - F2
      F22.base <- 1 - F22
    } 
    
    theta <- 1 + exp(theta.star)
    Gumbel1 <- exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))
    Gumbel2 <- exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))
    
    Gumbel1 <- ifelse(Gumbel1>precision, Gumbel1, precision)
    Gumbel1 <- ifelse(Gumbel1<(1-precision), Gumbel1, 1-precision)
    Gumbel2 <- ifelse(Gumbel2>precision, Gumbel2, precision)
    Gumbel2 <- ifelse(Gumbel2<(1-precision), Gumbel2, 1-precision)
    lx <- f2 - Gumbel1 + Gumbel2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*((-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(-log(F1.base))^(theta-1))*(-1/F1.base) 
    dcop11 <- exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*((-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(-log(F1.base))^(theta-1))*(-1/F1.base)  
    
    dcop2 <- exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*((-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(-log(F2.base))^(theta-1))*(-1/F2.base)
    dcop22 <- exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*((-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(-log(F22.base))^(theta-1))*(-1/F22.base) 
    
    dcop.theta1 <- -exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)*((-1/theta^2)*log((-log(F1.base))^(theta)+(-log(F2.base))^(theta)))+(1/theta)*(1/((-log(F1.base))^(theta)+(-log(F2.base))^(theta)))*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F2.base))^(theta)*log(-log(F2.base)))*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)*(-exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)))
    dcop.theta2 <- -exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)*((-1/theta^2)*log((-log(F1.base))^(theta)+(-log(F22.base))^(theta)))+(1/theta)*(1/((-log(F1.base))^(theta)+(-log(F22.base))^(theta)))*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F22.base))^(theta)*log(-log(F22.base)))*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)*(-exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)))
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- dcop1*((-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(-log(F1.base))^(theta-1))*(-1/F1.base)+exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-2)*theta^2*(-log(F1.base))^(2*theta-2)*(-1/F1.base)^2+((-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F1.base))^(theta-2)*(-1/F1.base)^2+theta*(-log(F1.base))^(theta-1)*(-1/F1.base)^2) )
    d2Cdcop112 <- dcop11*((-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(-log(F1.base))^(theta-1))*(-1/F1.base)+exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-2)*theta^2*(-log(F1.base))^(2*theta-2)*(-1/F1.base)^2+((-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F1.base))^(theta-2)*(-1/F1.base)^2+theta*(-log(F1.base))^(theta-1)*(-1/F1.base)^2) )
    
    d2Cdcop22 <- dcop2*((-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(-log(F2.base))^(theta-1))*(-1/F2.base)+exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-2)*theta^2*(-log(F2.base))^(2*theta-2)*(-1/F2.base)^2+((-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F2.base))^(theta-2)*(-1/F2.base)^2+theta*(-log(F2.base))^(theta-1)*(-1/F2.base)^2) )
    d2Cdcop222 <- dcop22*((-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(-log(F22.base))^(theta-1))*(-1/F22.base)+exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-2)*theta^2*(-log(F22.base))^(2*theta-2)*(-1/F22.base)^2+((-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F22.base))^(theta-2)*(-1/F22.base)^2+theta*(-log(F22.base))^(theta-1)*(-1/F22.base)^2) )
    
    d2Cdcop1dcop2 <- dcop2*((-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(-log(F1.base))^(theta-1))*(-1/F1.base)+exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-2)*theta*(-log(F2.base))^(theta-1)*(-1/F2.base))*theta*(-log(F1.base))^(theta-1)*(-1/F1.base)
    d2Cdcop11dcop22 <- dcop22*((-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(-log(F1.base))^(theta-1))*(-1/F1.base)+exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-2)*theta*(-log(F22.base))^(theta-1)*(-1/F22.base))*theta*(-log(F1.base))^(theta-1)*(-1/F1.base)
    
    d2Cdcop.theta12 <- dcop.theta1*(-1/theta*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(-log(F1.base))^(theta-1))*(-1/F1.base)+(-1/F1.base)*exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*((1/theta^2*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)))*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(-1)-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-2)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F2.base))^(theta)*log(-log(F2.base)))))*(theta*(-log(F1.base))^(theta-1))+(-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1)*((-log(F1.base))^(theta-1)+theta*(-log(F1.base))^(theta-1)*(log(-log(F1.base)))) )
    d2Cdcop.theta22 <- dcop.theta2*(-1/theta*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(-log(F1.base))^(theta-1))*(-1/F1.base)+(-1/F1.base)*exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*((1/theta^2*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)))*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(-1)-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-2)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F22.base))^(theta)*log(-log(F22.base)))))*(theta*(-log(F1.base))^(theta-1))+(-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1)*((-log(F1.base))^(theta-1)+theta*(-log(F1.base))^(theta-1)*(log(-log(F1.base)))) )
    
    d2Cdcop1dcop.theta1 <- dcop.theta1*(-1/theta*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1))*(theta*(-log(F2.base))^(theta-1))*(-1/F2.base)+(-1/F2.base)*exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*((1/theta^2*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)))*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(-1)-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-2)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F2.base))^(theta)*log(-log(F2.base)))))*(theta*(-log(F2.base))^(theta-1))+(-1/theta)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta-1)*((-log(F2.base))^(theta-1)+theta*(-log(F2.base))^(theta-1)*(log(-log(F2.base)))) )
    d2Cdcop11dcop.theta2 <- dcop.theta2*(-1/theta*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1))*(theta*(-log(F22.base))^(theta-1))*(-1/F22.base)+(-1/F22.base)*exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*((1/theta^2*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)))*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(-1)-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-2)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F22.base))^(theta)*log(-log(F22.base)))))*(theta*(-log(F22.base))^(theta-1))+(-1/theta)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta-1)*((-log(F22.base))^(theta-1)+theta*(-log(F22.base))^(theta-1)*(log(-log(F22.base)))) )
    
    d2Cdcop2dcop.theta1 <- -dcop.theta1*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)*(-1/theta^2*log(((-log(F1.base))^(theta)+(-log(F2.base))^(theta))))-exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta))*(dcop.theta1/(-exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log(F1.base))^(theta)+(-log(F2.base))^(theta)))+((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)*(2/theta^3*log(((-log(F1.base))^(theta)+(-log(F2.base))^(theta)))-1/theta^2*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(-1)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F2.base))^(theta)*log(-log(F2.base)))))+(-1/theta^2*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(-1)-1/theta*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(-2)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F2.base))^(theta)*log(-log(F2.base))))*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F2.base))^(theta)*log(-log(F2.base)))*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)*(-Gumbel1)+(-Gumbel1*(((-log(F1.base))^(theta)*log(-log(F1.base))^2+(-log(F2.base))^(theta)*log(-log(F2.base))^2)*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)+dcop.theta1/(-exp(-((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)))*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F2.base))^(theta)*log(-log(F2.base))))+((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(1/theta)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F2.base))^(theta)*log(-log(F2.base)))*(-dcop.theta1))*1/theta*((-log(F1.base))^(theta)+(-log(F2.base))^(theta))^(-1)
    d2Cdcop22dcop.theta2 <- -dcop.theta2*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)*(-1/theta^2*log(((-log(F1.base))^(theta)+(-log(F22.base))^(theta))))-exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta))*(dcop.theta2/(-exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log(F1.base))^(theta)+(-log(F22.base))^(theta)))+((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)*(2/theta^3*log(((-log(F1.base))^(theta)+(-log(F22.base))^(theta)))-1/theta^2*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(-1)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F22.base))^(theta)*log(-log(F22.base)))))+(-1/theta^2*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(-1)-1/theta*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(-2)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F22.base))^(theta)*log(-log(F22.base))))*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F22.base))^(theta)*log(-log(F22.base)))*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)*(-Gumbel2)+(-Gumbel2*(((-log(F1.base))^(theta)*log(-log(F1.base))^2+(-log(F22.base))^(theta)*log(-log(F22.base))^2)*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)+dcop.theta2/(-exp(-((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)))*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F22.base))^(theta)*log(-log(F22.base))))+((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(1/theta)*((-log(F1.base))^(theta)*log(-log(F1.base))+(-log(F22.base))^(theta)*log(-log(F22.base)))*(-dcop.theta2))*1/theta*((-log(F1.base))^(theta)+(-log(F22.base))^(theta))^(-1)
    
    
    theta.append <- theta - 1
    theta.append.der <- theta - 1
  
    
  } else if (VC$BivD %in% c("J0", "J90", "J180", "J270")) {
    
    
    if (VC$BivD == "J0") {
      F1.base <- F1
      F2.base <- F2
      F22.base <- F22
    } else if (VC$BivD == "J90") {
      F1.base <- 1 - F1
      F2.base <- F2
      F22.base <- F22
    } else if (VC$BivD == "J180") {
      F1.base <- 1 - F1
      F2.base <- 1 - F2
      F22.base <- 1 - F22
    } else if (VC$BivD == "J270") {
      F1.base <- F1
      F2.base <- 1 - F2
      F22.base <- 1 - F22
    } 
    
    theta <- 1 + exp(theta.star) + precision
    Joe1 <- 1-((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta)
    Joe2 <- 1-((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta)
    
    Joe1 <- ifelse(Joe1>precision, Joe1, precision)
    Joe1 <- ifelse(Joe1<(1-precision), Joe1, 1-precision)
    Joe2 <- ifelse(Joe2>precision, Joe2, precision)
    Joe2 <- ifelse(Joe2<(1-precision), Joe2, 1-precision)
    lx <- f2 - Joe1 + Joe2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1<- (-1/theta)*((1-F1.base)^theta+(1-F2.base)^theta-(1-F1.base)^theta*(1-F2.base)^theta)^(1/theta-1)*(-theta*(1-F1.base)^(theta-1)+(1-F2.base)^theta*theta*(1-F1.base)^(theta-1))
    dcop11<- (-1/theta)*((1-F1.base)^theta+(1-F22.base)^theta-(1-F1.base)^theta*(1-F22.base)^theta)^(1/theta-1)*(-theta*(1-F1.base)^(theta-1)+(1-F22.base)^theta*theta*(1-F1.base)^(theta-1))   
    
    dcop2<- (-1/theta)*((1-F1.base)^theta+(1-F2.base)^theta-(1-F1.base)^theta*(1-F2.base)^theta)^(1/theta-1)*(-theta*(1-F2.base)^(theta-1)+(1-F1.base)^theta*theta*(1-F2.base)^(theta-1))
    dcop22<- (-1/theta)*((1-F1.base)^theta+(1-F22.base)^theta-(1-F1.base)^theta*(1-F22.base)^theta)^(1/theta-1)*(-theta*(1-F22.base)^(theta-1)+(1-F1.base)^theta*theta*(1-F22.base)^(theta-1)) 
    
    dcop.theta1<- -(-1/theta^2)*log((1-F1.base)^theta+(1-F2.base)^theta-(1-F1.base)^theta*(1-F2.base)^theta)*((1-F1.base)^theta+(1-F2.base)^theta-(1-F1.base)^theta*(1-F2.base)^theta)^(1/theta)-(1/theta)*((1-F1.base)^theta+(1-F2.base)^theta-(1-F1.base)^theta*(1-F2.base)^theta)^(1/theta-1)*((1-F1.base)^theta*log(1-F1.base)+(1-F2.base)^theta*log(1-F2.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F2.base)^theta-(1-F2.base)^theta*log(1-F2.base)*(1-F1.base)^theta) 
    dcop.theta2<- -(-1/theta^2)*log((1-F1.base)^theta+(1-F22.base)^theta-(1-F1.base)^theta*(1-F22.base)^theta)*((1-F1.base)^theta+(1-F22.base)^theta-(1-F1.base)^theta*(1-F22.base)^theta)^(1/theta)-(1/theta)*((1-F1.base)^theta+(1-F22.base)^theta-(1-F1.base)^theta*(1-F22.base)^theta)^(1/theta-1)*((1-F1.base)^theta*log(1-F1.base)+(1-F22.base)^theta*log(1-F22.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F22.base)^theta-(1-F22.base)^theta*log(1-F22.base)*(1-F1.base)^theta)
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- -1/theta*(1/theta-1)*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-2)*(-theta*(1-F1.base)^(theta-1)+(1-F2.base)^(theta)*theta*(1-F1.base)^(theta-1))^2-1/theta*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F1.base)^(theta-2)-(1-F2.base)^(theta)*theta*(theta-1)*(1-F1.base)^(theta-2))
    d2Cdcop112 <- -1/theta*(1/theta-1)*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-2)*(-theta*(1-F1.base)^(theta-1)+(1-F22.base)^(theta)*theta*(1-F1.base)^(theta-1))^2-1/theta*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F1.base)^(theta-2)-(1-F22.base)^(theta)*theta*(theta-1)*(1-F1.base)^(theta-2))
    
    d2Cdcop22 <- -1/theta*(1/theta-1)*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-2)*(-theta*(1-F2.base)^(theta-1)+(1-F1.base)^(theta)*theta*(1-F2.base)^(theta-1))^2-1/theta*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F2.base)^(theta-2)-(1-F1.base)^(theta)*theta*(theta-1)*(1-F2.base)^(theta-2))
    d2Cdcop222 <- -1/theta*(1/theta-1)*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-2)*(-theta*(1-F22.base)^(theta-1)+(1-F1.base)^(theta)*theta*(1-F22.base)^(theta-1))^2-1/theta*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F22.base)^(theta-2)-(1-F1.base)^(theta)*theta*(theta-1)*(1-F22.base)^(theta-2))
    
    d2Cdcop1dcop2 <- -1/theta*(1/theta-1)*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-2)*(-theta*(1-F2.base)^(theta-1)+(1-F1.base)^(theta)*theta*(1-F2.base)^(theta-1))*(-theta*(1-F1.base)^(theta-1)+(1-F2.base)^(theta)*theta*(1-F1.base)^(theta-1))-1/theta*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*(-theta^2*(1-F2.base)^(theta-1)*(1-F1.base)^(theta-1))
    d2Cdcop11dcop22 <- -1/theta*(1/theta-1)*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-2)*(-theta*(1-F22.base)^(theta-1)+(1-F1.base)^(theta)*theta*(1-F22.base)^(theta-1))*(-theta*(1-F1.base)^(theta-1)+(1-F22.base)^(theta)*theta*(1-F1.base)^(theta-1))-1/theta*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*(-theta^2*(1-F22.base)^(theta-1)*(1-F1.base)^(theta-1))
    
    d2Cdcop.theta12 <- 1/theta^2*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*(-theta*(1-F1.base)^(theta-1)+(1-F2.base)^(theta)*theta*(1-F1.base)^(theta-1))-1/theta*((-theta*(1-F1.base)^(theta-1)+(1-F2.base)^(theta)*theta*(1-F1.base)^(theta-1))*(-dcop.theta1*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(-1)-((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-2)*((1-F1.base)^theta*log(1-F1.base)+(1-F2.base)^theta*log(1-F2.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F2.base)^theta-(1-F2.base)^theta*log(1-F2.base)*(1-F1.base)^theta))+((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*(-(1-F1.base)^(theta-1)-theta*(1-F1.base)^(theta-1)*log(1-F1.base)+theta*(1-F2.base)^(theta)*log(1-F2.base)*(1-F1.base)^(theta-1)+((1-F1.base)^(theta-1)+theta*(1-F1.base)^(theta-1)*log(1-F1.base))*(1-F2.base)^theta))
    d2Cdcop.theta22 <- 1/theta^2*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*(-theta*(1-F1.base)^(theta-1)+(1-F22.base)^(theta)*theta*(1-F1.base)^(theta-1))-1/theta*((-theta*(1-F1.base)^(theta-1)+(1-F22.base)^(theta)*theta*(1-F1.base)^(theta-1))*(-dcop.theta2*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(-1)-((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-2)*((1-F1.base)^theta*log(1-F1.base)+(1-F22.base)^theta*log(1-F22.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F22.base)^theta-(1-F22.base)^theta*log(1-F22.base)*(1-F1.base)^theta))+((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*(-(1-F1.base)^(theta-1)-theta*(1-F1.base)^(theta-1)*log(1-F1.base)+theta*(1-F22.base)^(theta)*log(1-F22.base)*(1-F1.base)^(theta-1)+((1-F1.base)^(theta-1)+theta*(1-F1.base)^(theta-1)*log(1-F1.base))*(1-F22.base)^theta))
    
    d2Cdcop1dcop.theta1 <- 1/theta^2*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*(-theta*(1-F2.base)^(theta-1)+(1-F1.base)^(theta)*theta*(1-F2.base)^(theta-1))-1/theta*((-theta*(1-F2.base)^(theta-1)+(1-F1.base)^(theta)*theta*(1-F2.base)^(theta-1))*(-dcop.theta1*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(-1)-((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-2)*((1-F1.base)^theta*log(1-F1.base)+(1-F2.base)^theta*log(1-F2.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F2.base)^theta-(1-F2.base)^theta*log(1-F2.base)*(1-F1.base)^theta))+((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*(-(1-F2.base)^(theta-1)-theta*(1-F2.base)^(theta-1)*log(1-F2.base)+theta*(1-F1.base)^(theta)*log(1-F1.base)*(1-F2.base)^(theta-1)+((1-F2.base)^(theta-1)+theta*(1-F2.base)^(theta-1)*log(1-F2.base))*(1-F1.base)^theta))
    d2Cdcop11dcop.theta2 <- 1/theta^2*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*(-theta*(1-F22.base)^(theta-1)+(1-F1.base)^(theta)*theta*(1-F22.base)^(theta-1))-1/theta*((-theta*(1-F22.base)^(theta-1)+(1-F1.base)^(theta)*theta*(1-F22.base)^(theta-1))*(-dcop.theta2*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(-1)-((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-2)*((1-F1.base)^theta*log(1-F1.base)+(1-F22.base)^theta*log(1-F22.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F22.base)^theta-(1-F22.base)^theta*log(1-F22.base)*(1-F1.base)^theta))+((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*(-(1-F22.base)^(theta-1)-theta*(1-F22.base)^(theta-1)*log(1-F22.base)+theta*(1-F1.base)^(theta)*log(1-F1.base)*(1-F22.base)^(theta-1)+((1-F22.base)^(theta-1)+theta*(1-F22.base)^(theta-1)*log(1-F22.base))*(1-F1.base)^theta))
    
    d2Cdcop2dcop.theta1 <- -2/theta^3*log(((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta)))*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta)+1/theta^2*(1/((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))*((1-F1.base)^theta*log(1-F1.base)+(1-F2.base)^theta*log(1-F2.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F2.base)^theta-(1-F2.base)^theta*log(1-F2.base)*(1-F1.base)^theta)*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta)+log(((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta)))*(-dcop.theta1))+1/theta^2*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*((1-F1.base)^theta*log(1-F1.base)+(1-F2.base)^theta*log(1-F2.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F2.base)^theta-(1-F2.base)^theta*log(1-F2.base)*(1-F1.base)^theta)-1/theta*(((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-1)*((1-F1.base)^theta*log(1-F1.base)^2+(1-F2.base)^theta*log(1-F2.base)^2-(1-F1.base)^theta*log(1-F1.base)^2*(1-F2.base)^theta-(1-F2.base)^theta*log(1-F2.base)^2*(1-F1.base)^theta-2*(1-F2.base)^theta*log(1-F2.base)*(1-F1.base)^theta*log(1-F1.base))+((1-F1.base)^theta*log(1-F1.base)+(1-F2.base)^theta*log(1-F2.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F2.base)^theta-(1-F2.base)^theta*log(1-F2.base)*(1-F1.base)^theta)*(-dcop.theta1*((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(-1)-((1-F1.base)^(theta)+(1-F2.base)^(theta)-(1-F1.base)^(theta)*(1-F2.base)^(theta))^(1/theta-2)*((1-F1.base)^theta*log(1-F1.base)+(1-F2.base)^theta*log(1-F2.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F2.base)^theta-(1-F2.base)^theta*log(1-F2.base)*(1-F1.base)^theta)))
    d2Cdcop22dcop.theta2 <- -2/theta^3*log(((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta)))*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta)+1/theta^2*(1/((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))*((1-F1.base)^theta*log(1-F1.base)+(1-F22.base)^theta*log(1-F22.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F22.base)^theta-(1-F22.base)^theta*log(1-F22.base)*(1-F1.base)^theta)*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta)+log(((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta)))*(-dcop.theta2))+1/theta^2*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*((1-F1.base)^theta*log(1-F1.base)+(1-F22.base)^theta*log(1-F22.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F22.base)^theta-(1-F22.base)^theta*log(1-F22.base)*(1-F1.base)^theta)-1/theta*(((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-1)*((1-F1.base)^theta*log(1-F1.base)^2+(1-F22.base)^theta*log(1-F22.base)^2-(1-F1.base)^theta*log(1-F1.base)^2*(1-F22.base)^theta-(1-F22.base)^theta*log(1-F22.base)^2*(1-F1.base)^theta-2*(1-F22.base)^theta*log(1-F22.base)*(1-F1.base)^theta*log(1-F1.base))+((1-F1.base)^theta*log(1-F1.base)+(1-F22.base)^theta*log(1-F22.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F22.base)^theta-(1-F22.base)^theta*log(1-F22.base)*(1-F1.base)^theta)*(-dcop.theta2*((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(-1)-((1-F1.base)^(theta)+(1-F22.base)^(theta)-(1-F1.base)^(theta)*(1-F22.base)^(theta))^(1/theta-2)*((1-F1.base)^theta*log(1-F1.base)+(1-F22.base)^theta*log(1-F22.base)-(1-F1.base)^theta*log(1-F1.base)*(1-F22.base)^theta-(1-F22.base)^theta*log(1-F22.base)*(1-F1.base)^theta)))
    

    theta.append <- theta-1-precision
    theta.append.der <- theta-1-precision
  
    
  } 
  
  
  
  if (VC$BivD %in% c("C90", "J90", "G90")) {
    
    if (VC$BivD == "C90") {
      Copula1 <- Clayton1
      Copula2 <- Clayton2
    } else if (VC$BivD == "J90") {
      Copula1 <- Joe1
      Copula2 <- Joe2
    } else if (VC$BivD == "G90") {
      Copula1 <- Gumbel1
      Copula2 <- Gumbel2
    } 
    

    Copula1 <- F2 - Copula1
    Copula2 <- F22 - Copula2
    Copula1 <- ifelse(Copula1>precision, Copula1, precision)
    Copula1 <- ifelse(Copula1<(1-precision), Copula1, 1-precision)
    Copula2 <- ifelse(Copula2>precision, Copula2, precision)
    Copula2 <- ifelse(Copula2<(1-precision), Copula2, 1-precision)
    lx <- f2 - Copula1 + Copula2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- dcop1
    dcop11 <- dcop11    
    
    dcop2 <- 1 - dcop2
    dcop22 <- 1 - dcop22
    
    dcop.theta1 <- - ( - dcop.theta1)
    dcop.theta2 <- - ( - dcop.theta2)
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- - d2Cdcop12
    d2Cdcop112 <- - d2Cdcop112
    
    d2Cdcop22 <- - d2Cdcop22
    d2Cdcop222 <- - d2Cdcop222
    
    d2Cdcop1dcop2 <- d2Cdcop1dcop2 
    d2Cdcop11dcop22 <- d2Cdcop11dcop22 
    
    d2Cdcop.theta12 <- - d2Cdcop.theta12
    d2Cdcop.theta22 <- - d2Cdcop.theta22
    
    d2Cdcop1dcop.theta1 <- - (- d2Cdcop1dcop.theta1)
    d2Cdcop11dcop.theta2 <- - (- d2Cdcop11dcop.theta2)
    
    d2Cdcop2dcop.theta1 <- - d2Cdcop2dcop.theta1
    d2Cdcop22dcop.theta2 <- - d2Cdcop22dcop.theta2
    

    
    theta.append <- -theta.append
    theta.append.der <- -theta.append.der
  
    
  } 
  
  if (VC$BivD %in% c("C180", "J180", "G180")) {
    
    if (VC$BivD == "C180") {
      Copula1 <- Clayton1
      Copula2 <- Clayton2
    } else if (VC$BivD == "J180") {
      Copula1 <- Joe1
      Copula2 <- Joe2
    } else if (VC$BivD == "G180") {
      Copula1 <- Gumbel1
      Copula2 <- Gumbel2
    } 
    
    
    Copula1 <- F1 + F2 - 1 + Copula1
    Copula2 <- F1 + F22 - 1 + Copula2
    Copula1 <- ifelse(Copula1>precision, Copula1, precision)
    Copula1 <- ifelse(Copula1<(1-precision), Copula1, 1-precision)
    Copula2 <- ifelse(Copula2>precision, Copula2, precision)
    Copula2 <- ifelse(Copula2<(1-precision), Copula2, 1-precision)
    lx <- f2 - Copula1 + Copula2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- 1 - dcop1
    dcop11 <- 1 - dcop11    
    
    dcop2 <- 1 - dcop2
    dcop22 <- 1 - dcop22
    
    dcop.theta1 <- ( dcop.theta1)
    dcop.theta2 <- ( dcop.theta2)
    
    
    # Hessian derivative components
    
    d2Cdcop12 <-  d2Cdcop12
    d2Cdcop112 <-  d2Cdcop112
    
    d2Cdcop22 <-  d2Cdcop22
    d2Cdcop222 <-  d2Cdcop222
    
    d2Cdcop1dcop2 <- d2Cdcop1dcop2 
    d2Cdcop11dcop22 <- d2Cdcop11dcop22 
    
    d2Cdcop.theta12 <- - d2Cdcop.theta12
    d2Cdcop.theta22 <- - d2Cdcop.theta22
    
    d2Cdcop1dcop.theta1 <- - d2Cdcop1dcop.theta1
    d2Cdcop11dcop.theta2 <- - d2Cdcop11dcop.theta2
    
    d2Cdcop2dcop.theta1 <-  d2Cdcop2dcop.theta1
    d2Cdcop22dcop.theta2 <-  d2Cdcop22dcop.theta2
    

  
    theta.append <- theta.append
    theta.append.der <- theta.append.der
  
    
  } 
  
  if (VC$BivD %in% c("C270", "J270", "G270")) {
    
    if (VC$BivD == "C270") {
      Copula1 <- Clayton1
      Copula2 <- Clayton2
    } else if (VC$BivD == "J270") {
      Copula1 <- Joe1
      Copula2 <- Joe2
    } else if (VC$BivD == "G270") {
      Copula1 <- Gumbel1
      Copula2 <- Gumbel2
    } 
    
    
    Copula1 <- F1 - Copula1
    Copula2 <- F1 - Copula2
    Copula1 <- ifelse(Copula1>precision, Copula1, precision)
    Copula1 <- ifelse(Copula1<(1-precision), Copula1, 1-precision)
    Copula2 <- ifelse(Copula2>precision, Copula2, precision)
    Copula2 <- ifelse(Copula2<(1-precision), Copula2, 1-precision)
    lx <- f2 - Copula1 + Copula2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- 1 - dcop1
    dcop11 <- 1 - dcop11    
    
    dcop2 <- dcop2
    dcop22 <- dcop22
    
    dcop.theta1 <- - ( - dcop.theta1)
    dcop.theta2 <- - ( - dcop.theta2)
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- - d2Cdcop12
    d2Cdcop112 <- - d2Cdcop112
    
    d2Cdcop22 <- - d2Cdcop22
    d2Cdcop222 <- - d2Cdcop222
    
    d2Cdcop1dcop2 <- d2Cdcop1dcop2 
    d2Cdcop11dcop22 <- d2Cdcop11dcop22 
    
    d2Cdcop.theta12 <- - ( - d2Cdcop.theta12 )
    d2Cdcop.theta22 <- - ( - d2Cdcop.theta22 )
    
    d2Cdcop1dcop.theta1 <- - ( d2Cdcop1dcop.theta1 )
    d2Cdcop11dcop.theta2 <- - ( d2Cdcop11dcop.theta2 )
    
    d2Cdcop2dcop.theta1 <- - d2Cdcop2dcop.theta1
    d2Cdcop22dcop.theta2 <- - d2Cdcop22dcop.theta2
    
    
    theta.append <- -theta.append
    theta.append.der <- -theta.append.der
  
    
  } 
  
  
  list(theta=theta, lx=lx,
       dcop1=dcop1, dcop11=dcop11, 
       dcop2=dcop2, dcop22=dcop22, 
       dcop.theta1=dcop.theta1, dcop.theta2=dcop.theta2, 
       d2Cdcop12=d2Cdcop12, d2Cdcop112=d2Cdcop112, 
       d2Cdcop22=d2Cdcop22, d2Cdcop222=d2Cdcop222,
       d2Cdcop1dcop2=d2Cdcop1dcop2, d2Cdcop11dcop22=d2Cdcop11dcop22, 
       d2Cdcop.theta12=d2Cdcop.theta12, d2Cdcop.theta22=d2Cdcop.theta22, 
       d2Cdcop1dcop.theta1=d2Cdcop1dcop.theta1, d2Cdcop11dcop.theta2=d2Cdcop11dcop.theta2,
       d2Cdcop2dcop.theta1=d2Cdcop2dcop.theta1, d2Cdcop22dcop.theta2=d2Cdcop22dcop.theta2,
       theta.append=theta.append, theta.append.der=theta.append.der
    )
  
  
  
}




