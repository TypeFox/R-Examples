NPHMC<-function(n=NULL, power=0.8, alpha=0.05, accrualtime=NULL, followuptime=NULL, p=0.5, 
            accrualdist=c("uniform","increasing","decreasing"),
		hazardratio=NULL, oddsratio=NULL, pi0=NULL, survdist=c("exp","weib"), k=1, lambda0=NULL, data=NULL){
	N<-list()
	class(N) <- c("NPHMC") 
 	if (is.null(n) & is.null(power)) 
             stop("'n' and 'power' cannot be missing simultaneously! User must input 'n' to calculate power
                   input 'power' to calculate sample size 'n'.")
      if (is.null(data)){
      if (hazardratio<=0) stop("Hazardratio must be greater than 0.")
	if (hazardratio==1) stop("Hazardratio cannot be 1 since beta is the denominator of the sample size formula and cannot be 0.")
	if (oddsratio<0) stop("Oddsratio cannot be less than 0.")
        if (pi0==0 | oddsratio==0) {
             i1 <- integrate(f1,0,followuptime,survdist,k,lambda0)$value
	       i2 <- integrate(f2,followuptime,(accrualtime+followuptime),accrualtime,followuptime,accrualdist,survdist,k,lambda0)$value	  
	       beta0 <- log(hazardratio)
         	 pdeath <- i1+i2
            if (is.null(n)) {
	       nsizeph <- ceiling((qnorm(power)-qnorm(alpha/2))^2/(p*(1-p)*beta0^2*pdeath))
             cat("====================================================================== \n")
	       cat("SAMPLE SIZE CALCULATION BASED ON STANDARD PH MODEL (NO CURE FRACTION) \n")
	       cat("====================================================================== \n")
		 cat("At alpha =",alpha,"and power =",power,":","\n")
             cat("Standard PH Model: n =",nsizeph,"\n")}
           else {
             pw <- round(pnorm(sqrt(n*p*(1-p)*beta0^2*pdeath)+qnorm(alpha/2)), digits = 2) 
             cat("====================================================================== \n")
	       cat("POWER CALCULATION BASED ON STANDARD PH MODEL (NO CURE FRACTION) \n")
	       cat("====================================================================== \n") 
             cat("Standard PH Model: Power =",pw,"\n")
               }
               }
        else {
	  i1 <- integrate(f1,0,followuptime,survdist,k,lambda0)$value
	  i2 <- integrate(f2,followuptime,(accrualtime+followuptime),accrualtime,followuptime,accrualdist,survdist,k,lambda0)$value	  
	  beta0 <- log(hazardratio)
	  gamma0 <- log(oddsratio)
	  i3 <- integrate(f3,0,followuptime,beta0,gamma0,pi0,survdist,k,lambda0)$value
	  i4 <- integrate(f4,followuptime,(accrualtime+followuptime),accrualtime,followuptime,accrualdist,beta0,gamma0,pi0,survdist,k,lambda0)$value
      pdeath <- i1+i2
     if (is.null(n)) {
      nsize <- ceiling((qnorm(power)-qnorm(alpha/2))^2*(i1+i2)/((i3+i4)^2*p*(1-p)*(1-pi0)*beta0^2))
	nsizeph <- ceiling((qnorm(power)-qnorm(alpha/2))^2/(p*(1-p)*beta0^2*pdeath))
	cat("\n")
	cat("======================================================================== \n")
	cat("SAMPLE SIZE CALCULATION FOR PH MIXTURE CURE MODEL AND STANDARD PH MODEL \n")
	cat("======================================================================== \n")
		 cat("At alpha =",alpha,"and power =",power,":","\n") 
	cat("PH Mixture Cure Model: n =",nsize,"\n")
	cat("Standard PH Model: n =",nsizeph,"\n")
      N$nsize <- nsize
	N$nsizeph <- nsizeph }
      else {
             pw <- round(pnorm(sqrt(n*((i3+i4)^2*p*(1-p)*(1-pi0)*beta0^2)/(i1+i2))+qnorm(alpha/2)), digits = 2)
             pwph <- round(pnorm(sqrt(n*p*(1-p)*beta0^2*pdeath)+qnorm(alpha/2)), digits = 2)
             cat("====================================================================== \n")
	       cat("POWER CALCULATION FOR PH MIXTURE CURE MODEL AND STANDARD PH MODEL \n")
	       cat("====================================================================== \n")
             cat("PH Mixture Cure Model: Power =",pw,"\n")
             cat("Standard PH Model: Power =",pwph,"\n")
		if (length(n)>1) {
		 plot(n, pw, type="o", main="Power Analysis \n PH Cure Model vs. Standard PH Model",
			sub="Solid line - PH Cure Model; Dash Line - Standard PH Model",
 			xlab="sample size", ylab="Power", lty=1, col="blue", xlim=c(0, max(n)), ylim=c(0, 1))
		 lines(n, pwph, type="o", pch=22, lty=2, col="red") }
		 N$pw <- pw
		 N$pwph <- pwph
               }
	}
     }
 
  if (!is.null(data)){  
  if(try(!is.null(hazardratio) | !is.null(oddsratio))){
      stop("The 'hazardratio' and 'oddsratio' are not needed when data is specified.")}     	 
	ta=accrualtime
	tf=followuptime
	ttot=ta+tf
	t<-data[,1]
      colnames(data)<-c("Time","Status","X")
	Time=data[,1]
	Status=data[,2]
	X=data[,3]
      f=smcure(Surv(Time, Status)~X,~X,data=data,model="ph",Var=FALSE)
      N$f <- f
	time<-sort(t[Status==1])
	beta0nocure <- coxph(Surv(Time, Status)~X,method="breslow", data=data)$coef
	death_point <- sort(unique(subset(Time, Status==1)))
	coxexp <- exp(beta0nocure*X)
	lambda <- numeric()
    	event <- numeric()
      for(i in 1: length(death_point)){
       event[i] <- sum(Status*as.numeric(Time==death_point[i]))
        temp <- sum(as.numeric(Time>=death_point[i])*Status*drop(coxexp))
       	temp1 <- event[i]
       lambda[i] <- temp1/temp
        }
    HHazard <- numeric()
    for(i in 1:length(Time)){
        HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
        if(Time[i]>max(death_point))HHazard[i] <- Inf
        if(Time[i]<min(death_point))HHazard[i] <- 0
        }
 	 	snocure <- exp(-HHazard)
		beta0 <- f$beta
		gamma0 <- -f$b[2]
		pi0=1-exp(f$b[1])/(1 + exp(f$b[1]))
    		s=sort(f$s[Status==1],decreasing = TRUE)
		snocure <- sort(snocure[Status==1],decreasing = TRUE)
    		 f0<-diff(s)
     		 f0nocure <- diff(snocure)
    		 s1 <- sum(-f0*as.numeric(time<=tf)[-1])
       	 s1nocure <- sum(-f0nocure*as.numeric(time<=tf)[-1])
		sc=(ta+tf-time)/ta
	 	 s2 <- sum(-diff(s)*sc[-length(sc)]*as.numeric(time>tf)[-1])
       	 s2nocure <- sum(-diff(snocure)*sc[-length(sc)]*as.numeric(time>tf)[-1])
		Spop=pi0+(1-pi0)*s
		m=(gamma0/beta0-log(s))*pi0/Spop-1
		s3 <- sum(-diff(s)*m[-length(m)]*as.numeric(time<=tf)[-1])
		Spop4=pi0+(1-pi0)*s
		m4=(gamma0/beta0-log(s))*pi0/Spop4-1
	  	s4 <- sum(-diff(s)*m4[-length(m4)]*sc[-length(sc)]*as.numeric((time>tf) & (time<=ttot))[-1])
		pdeathNonpar <- s1+s2
		pdeathNonpar <- s1nocure+s2nocure

      if (is.null(n)) {
		nonpar=ceiling((qnorm(power)-qnorm(alpha/2))^2*(s1+s2)/((s3+s4)^2*p*(1-p)*(1-pi0)*beta0^2))
		N$nonpar<- nonpar
		N$HR <- exp(beta0) 
		N$OR <- exp(gamma0) 
		N$pi0<- pi0
		cat("\n")
		cat("======================================================================== \n")
		cat("SAMPLE SIZE CALCULATION FOR PH MIXTURE CURE MODEL AND STANDARD PH MODEL \n")
		cat("======================================================================== \n")
 		cat("At alpha =",alpha,"and power =",power,":","\n") 
		cat("PH Mixture Cure Model with KM estimators: n =",nonpar,"\n")
		nonparPH<- ceiling((qnorm(power)-qnorm(alpha/2))^2/(p*(1-p)*beta0nocure^2*pdeathNonpar))
		cat("Standard PH Model with KM estimators: n =",nonparPH,"\n")
		N$nonpar<- nonpar
		N$nonparPH<- nonparPH	}
       else {
             nonparpw <- round(pnorm(sqrt(n*((s3+s4)^2*p*(1-p)*(1-pi0)*beta0^2)/(s1+s2))+qnorm(alpha/2)), digits = 2) 
             nonparpwph <- round(pnorm(sqrt(n*p*(1-p)*beta0nocure^2*pdeathNonpar)+qnorm(alpha/2)), digits = 2)
             cat("====================================================================== \n")
	       cat("POWER CALCULATION FOR PH MIXTURE CURE MODEL AND STANDARD PH MODEL \n")
	       cat("====================================================================== \n")
             cat("PH Mixture Cure Model: Power =",nonparpw,"\n")
             cat("Standard PH Model: Power =",nonparpwph,"\n")
		if (length(n)>1) {
		 plot(n, nonparpw, type="o", main="Power Analysis (by existing historical data) \n PH Cure Model vs. Standard PH Model",
			sub="Solid line - PH Cure Model; Dash Line - Standard PH Model",
 			xlab="sample size", ylab="Power", lty=1, col="blue", xlim=c(0, max(n)), ylim=c(0, 1))
		 lines(n, nonparpwph, type="o", pch=22, lty=2, col="red") }
		 N$nonparpw <- nonparpw 
		 N$nonparpwph <- nonparpwph 
         }    
	
  	}

invisible(N)
 }
