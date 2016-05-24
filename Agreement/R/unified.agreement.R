unified.agreement <-
function(dataset, var=NA, k, m, CCC_a=0.9775, CCC_a_intra=0.995, CCC_a_inter=0.9775, CCC_a_total=0.9775, CP_a=0.9, TDI_a=150, TDI_a_intra=75, TDI_a_inter=150, TDI_a_total=150, tran=0, error="prop", dec=4, alpha=0.025,iter=35,toe=1e-10){

if (sum(is.na(var))==0) dataset <- dataset[,var];

n <- dim(dataset)[1];
Inputs <- cbind(k,m,n,tran,error,dec);

y <- as.numeric(as.matrix(t(dataset)));
rep <- rep(1:m,length(y)/m);
method <- rep(1:k,each=m,len=length(y));
id <- rep(1:(length(y)/(k*m)),each=k*m);
rating <- (method-1)*m+rep;

ccc_intra <- NA;
ccc_inter <- NA;
ccc_total <- NA;
ccc_intra_1 <- NA;
ccc_inter_1 <- NA;
ccc_total_1 <- NA;
rho_intra <- NA;
rho_inter <- NA;
rho_total <- NA;
accuracy_intra <- NA;
accuracy_inter <- NA;
accuracy_total <- NA;
TDI_intra <- NA;
TDI_inter <- NA;
TDI_total <- NA;
cp_intra <- NA;
cp_inter <- NA;
cp_total <- NA;
ccc_intra_tran_lower <- NA;
ccc_inter_tran_lower <- NA;
ccc_total_tran_lower <- NA;
rho_intra_tran_lower <- NA;
rho_inter_tran_lower <- NA;
rho_total_tran_lower <- NA;
accuracy_intra_tran_lower <- NA;
accuracy_inter_tran_lower <- NA;
accuracy_total_tran_lower <- NA;
TDI_intra_tran_upper <- NA;
TDI_inter_tran_upper <- NA;
TDI_total_tran_upper <- NA;
cp_intra_tran_lower <- NA;
cp_inter_tran_lower <- NA;
cp_total_tran_lower <- NA;
RBS_intra <- NA;
RBS_inter <- NA;
RBS_total <- NA;

stat1 <- "Estimate";
ccca <- (1-alpha)*100;
stat2 <- paste(ccca,"% Conf. Limit",sep="");
stat3 <- "Allowance";

if (tran==1 & tolower(error)=="prop"){
	y <- log(y); 
	dataset <- log(dataset);
}

k1 <- k-1;
x <- array(0,c(length(y),k-1));
for (i in 1:(k-1)) x[method==i,i] <- 1;

if (m>1){
	var <- array(0,k);
	for (i in 1:k) var[i] <- mean(apply(dataset[((i-1)*m+1):(i*m)],1,var));
	cellvar <- mean(var);
}

nx <- dim(x)[2];
beta <- lm(y~x)$coefficients;

if (m==1){
	nbeta <- length(beta);
	n <- max(id); # total number of persons
	nc <- k*m;
	sigmaalpha <- 0;
	sigmaerror <- 0;
	sigmabeta <- 0;
	U_i <- 0;
	nsigma <- k+3;
	sigma <- rep(0,nsigma);
	Nst <- array(n,c(nc,nc));
	muk <- apply(dataset,2,sum)/diag(Nst);
	V <- var(dataset)*(Nst-1)/(Nst);
	for (s in 1:(k-1)){
		for (t in (s+1):k){
			sigmabeta <- sigmabeta+(muk[s]-muk[t])^2/(k*(k-1));
			sigmaalpha <- sigmaalpha+2*V[s,t]/(k*(k-1));
		}
	}
	sigmaerror <- sum(diag(V))/k-sigmaalpha;
	sigma <- muk;
	sigma[k+1] <- sigmabeta;
	sigma[k+2] <- sigmaalpha;
	sigma[k+3] <- sigmaerror;
	
	crit1 <- 1;
	it <- 1;
	while ((crit1>toe)&(it<iter)){
          U2 <- rep(0,nsigma);
          FHF <- array(0,c(nsigma,nsigma));
          u2sq <- FHF;
          for (i in 1:n){
		Y_i <- y[id==i];
		R_i <- rep[id==i];
		YY_i <- rep(0,nsigma);
		for (s in 1:k) YY_i[s] <- YY_i[s]+Y_i[s];
		for (s in 1:(k-1)){
			for (t in 2:k) YY_i[k+1] <- YY_i[k+1]+(Y_i[s]-Y_i[t])^2/(k*(k-1));
		}
		for (s in 1:k) YY_i[k+2]  <-  YY_i[k+2] + (Y_i[s]- muk[s])^2/k;
		for (s in 1:(k-1)){
			for (t in (s+1):k) YY_i[k+3]  <-  YY_i[k+3] + 2*((Y_i[s] - muk[s])*(Y_i[t] - muk[t]))/(k*(k-1));
		}
		delta_i <- muk;
		delta_i[k+1] <- sigmabeta+sigmaerror;
        	delta_i[k+2] <- sigmaalpha+sigmaerror;
        	delta_i[k+3] <- sigmaalpha;
		
		H_i <- diag(sigmaalpha+sigmaerror,nsigma);
		H_i[k+1,k+1] <- (4*sigmaerror^2+8*sigmabeta*sigmaerror)/(k*(k-1));
                H_i[k+2,k+2] <- (2*k*sigmaalpha^2+2*sigmaerror^2+4*sigmaalpha*sigmaerror)/k;
        	H_i[k+3,k+3] <- (2*k*(k-1)*sigmaalpha^2+2*k*sigmaerror*sigmaalpha+2*sigmaerror^2)/(k*(k-1));
		F_i <- diag(0,nsigma);
		for (s in (1:k)) F_i[s,s] <- 1;
		F_i[k+1,k+1] <- 1;
        	F_i[k+1,k+3] <- 1;
        	F_i[k+2,k+2] <- 1;
        	F_i[k+2,k+3] <- 1;
        	F_i[k+3,k+2] <- 1;

		u2_i <- t(F_i) %*% solve(H_i) %*% (YY_i-delta_i);
		U2 <- U2+u2_i;
        	FHF <- FHF+t(F_i) %*% solve(H_i) %*% F_i;
		u2sq <- u2sq + u2_i%*%t(u2_i);
          }
          DELTA2 <- solve(FHF,U2);
          crit1 <- max(abs(DELTA2))/max(sigma);
          sigma <- sigma+DELTA2;
          it <- it+1;
	}
	if ((it>=iter)&(crit1>toe)) warning("Iteration does not converge, soluations may be inaccurate!");
	var <- solve(FHF) %*% u2sq %*% t(solve(FHF));
	var_sigmabeta <- var[k+1,k+1];
	var_sigmaalpha <- var[k+2,k+2];
	var_sigmaerror <- var[k+3,k+3];
	cov_betaalpha <- var[k+1,k+2];
	cov_betaerror <- var[k+1,k+3];
	cov_alphaerror <- var[k+2,k+3];

	rho <- sigmaalpha/(sigmaalpha + sigmaerror);
	if (rho<=0.999999){
		rho_z <- log((1+rho)/(1-rho))/2;
	}else rho_z <- NA;
	accuracy <- (sigmaalpha + sigmaerror)/(sigmaalpha + sigmabeta + sigmaerror);
	if (accuracy<=0.999999){
		accuracy_l <- log(accuracy/(1-accuracy));
	}else accuracy_l <- NA;
	ccc <- sigmaalpha/(sigmaalpha + sigmabeta + sigmaerror);
	if (ccc<=0.999999){
		ccc_z <- log((1+ccc)/(1-ccc))/2;
	}else ccc_z <- NA
	esquare <- (2*sigmaerror+2*sigmabeta);

	if (tran==1){
		TDI <- qnorm(1-(1-CP_a)/2)*sqrt(esquare);
		TDI_w <- log(esquare);
		kp0_ <- TDI_a;
       	        if (error=="prop") kp0_ <- log(TDI_a/100+1);
		cp <- 1-2*(1-pnorm(kp0_/sqrt(esquare)));
		if (cp<=0.999999){
			cp_l <- log(cp/(1-cp));
		}else cp_l <- NA;
	}

	var_rho <- ((1-rho)^2*var_sigmaalpha+rho^2*var_sigmaerror-2*(1-rho)*rho*cov_alphaerror)/(sigmaalpha + sigmaerror)^2;
	if (var_rho<=0) var_rho <- NA;
	se_rho <- sqrt(var_rho);
	rho_notran_lower <- rho-qnorm(1-alpha)*se_rho;
	if (rho<=0.999999){
		se_rho_z <- se_rho/(1-rho^2);
	}else se_rho_z <- NA;
	rho_z_low <- rho_z-qnorm(1-alpha)*se_rho_z;
	rho_tran_lower <- (exp(2*rho_z_low)-1)/(exp(2*rho_z_low)+1);

	var_ccc <- ((1-ccc)^2*var_sigmaalpha+ccc^2*(var_sigmabeta+var_sigmaerror+2*cov_betaerror)-(2*(1-ccc)*ccc*(cov_betaalpha+cov_alphaerror)))/(sigmaalpha + sigmabeta + sigmaerror)^2;
	if (var_ccc<=0) var_ccc <- NA;
	se_ccc <- sqrt(var_ccc);
	ccc_notran_lower <- ccc-qnorm(1-alpha)*se_ccc;
	if (ccc<=0.999999){
		se_ccc_z <- se_ccc/(1-ccc^2);
	}else se_ccc_z <- NA;
	ccc_z_low <- ccc_z-qnorm(1-alpha)*se_ccc_z;
	ccc_tran_lower <- (exp(2*ccc_z_low)-1)/(exp(2*ccc_z_low)+1);

	var_accuracy <- ((1-accuracy)^2*(var_sigmaalpha+var_sigmaerror+2*cov_alphaerror)+accuracy^2*var_sigmabeta-2*(1-accuracy)*accuracy*(cov_betaalpha+cov_betaerror))/(sigmaalpha + sigmabeta + sigmaerror)^2;
	if (var_accuracy<=0) var_accuracy <- NA;
	se_accuracy <- sqrt(var_accuracy);
	accuracy_notran_lower <- accuracy-qnorm(1-alpha)*se_accuracy;
	if (accuracy<=0.999999){
		se_accuracy_l <- se_accuracy/(accuracy*(1-accuracy));
	}else se_accuracy_l <- NA;
	accuracy_l_low <- accuracy_l-qnorm(1-alpha)*se_accuracy_l;
	accuracy_tran_lower <- exp(accuracy_l_low)/(1+exp(accuracy_l_low));

	var_esquare <- 4*(var_sigmabeta+var_sigmaerror+2*cov_betaerror);
	if (var_esquare<=0) var_esquare <- NA;
	se_esquare <- sqrt(var_esquare);

	if (tran==1){
		se_TDI_w <- sqrt(var_esquare)/esquare;
		TDI_w_up <- TDI_w+qnorm(1-alpha)*se_TDI_w;
		TDI_e_up <- sqrt(exp(TDI_w_up));
		TDI_tran_upper <- qnorm(1-(1-CP_a)/2)*TDI_e_up;
		var_cp <- exp(-(kp0_)^2/(esquare))*(1+(kp0_)^2/(esquare))^2*var_esquare/(8*pi*(kp0_)^2*esquare);
		if (var_cp<=0) var_cp <- NA;
		se_cp <- sqrt(var_cp);
		cp_notran_lower <- cp-qnorm(1-alpha)*se_cp;
		if (cp<=0.999999){
			se_cp_l <- se_cp/(cp*(1-cp));
		}else se_cp_l <- NA;
		cp_l_low <- cp_l-qnorm(1-alpha)*se_cp_l;
		cp_tran_lower <- exp(cp_l_low)/(1+exp(cp_l_low));
		RBS <- sigmabeta/sigmaerror;
		
		if (error=="prop"){
			TDI <- 100*(exp(TDI)-1);
			TDI_tran_upper <- 100*(exp(TDI_tran_upper)-1);
		}

		TDI <- round(TDI,dec);
		TDI_tran_upper <- round(TDI_tran_upper,dec);

		Stat <- cbind(rbind(stat1,stat2,stat3),round(rbind(cbind(ccc,rho,accuracy,TDI,cp,RBS),cbind(ccc_tran_lower, rho_tran_lower, accuracy_tran_lower, TDI_tran_upper, cp_tran_lower, NA),cbind(CCC_a,NA,NA,TDI_a, CP_a,NA)),4));
	}
	if (tran==0){
		Stat <- cbind(rbind(stat1,stat2,stat3),round(rbind(cbind(ccc,rho,accuracy,NA,NA,NA),cbind(ccc_notran_lower, rho_notran_lower, accuracy_notran_lower, NA, NA, NA),cbind(CCC_a,NA,NA,NA, NA,NA)),4));
	}
	estimates <- list(Stat=Stat, Inputs=Inputs);
	class(estimates) <- "unified_agreement";

}else{
	nbeta <- length(beta);
	n <- max(id); # total number of persons
	nc <- k*m;
	sigmaalpha <- 0;
	sigmaerror <- cellvar;
	sigmabeta <- 0;
	sigmagamma <- 0;
	U_i <- 0;
	nsigma <- k+4;
	sigma <- rep(0,nsigma);

	Nst <- array(n,c(nc,nc));
	Ytotal <- apply(dataset,2,sum);
	Ybar <- Ytotal/diag(Nst);
	VC <- var(dataset)*(Nst-1)/(Nst); # covariance structure for columns
	V <- 0;
	if (k>2){
		for (i in 1:n){
			mu <- cbind(rep(1,k*m),x[id==i,]) %*% beta;
			V <- V+(y[id==i]-mu)%*%t(y[id==i]-mu);
		}
	}else{
		for (i in 1:n){
			mu <- cbind(rep(1,k*m),x[id==i]) %*% beta;
			V <- V+(y[id==i]-mu)%*%t(y[id==i]-mu);
		}
	}
	muk <- rep(0,k); # mean of method
	VM <- array(0,c(k,k)); # variance of method
	for (i in 1:k){
		muk[i] <- mean(y[method==i]);
		VM[i,i] <- var(y[method==i])*(length(y[method==i])-1)/length(y[method==i]);
	}
	sigmabeta <- 0;
	for (s in 1:(k-1)){
		for (t in (s+1):k) sigmabeta <- sigmabeta + (muk[s]-muk[t])^2/(k*(k-1));
	}
	sigmaalpha <- 0;
	for (s in 1:(k-1)){
		for (t in (s+1):k){
			for (l in 1:m){
				for (r in 1:m) sigmaalpha = sigmaalpha + 2*VC[(s-1)*m+l, (t-1)*m+r]/(m*m*k*(k-1));
			}
		}
	}
	VC1 <- sum(diag(VC))/(m^2*k);
	VC2 <- 0;
	for (s in 1:k){
		for (l in 1:(m-1)){
			for (r in (l+1):m) VC2 <- VC2 + 2*VC[(s-1)*m+l,(s-1)*m+r]/(m*m*k);
		}
	}
	sigmagamma  <-  VC1 + VC2 - sigmaerror/m - sigmaalpha;
	sigma <- muk;
	sigma[k+1] <- sigmabeta;
    	sigma[k+2] <- sigmaalpha;
    	sigma[k+3] <- sigmaerror;
    	sigma[k+4] <- sigmagamma;
	
#######################################

	crit1 <- 1;
	it <- 1;
	while ((crit1>toe)&(it<iter)){
	  U2 <- rep(0,nsigma);
	  FHF <- array(0,c(nsigma,nsigma));
	  u2sq <- FHF;
	  for (i in 1:n){
		Y_i <- y[id==i];
		R_i <- rep[id==i];
		YY_i <- rep(0,nsigma);
		for (s in 1:k){
			for (l in 1:m)	YY_i[s] <- YY_i[s] + Y_i[(s-1)*m+l];
			YY_i[s] <- YY_i[s]/m;
		}
		for (s in 1:(k-1)){
			for (t in 2:k) YY_i[k+1] <- YY_i[k+1] + (YY_i[s]-YY_i[t])^2/(k*(k-1));
		}
		for (s in 1:(k-1)){
			for (t in (s+1):k){
				for (l in 1:m){
					for (r in 1:m)  YY_i[k+2] <- YY_i[k+2] + 2*((Y_i[(s-1)*m+l] - Ybar[(s-1)*m+l])*(Y_i[(t-1)*m+r] - Ybar[(t-1)*m+r]))/(m*m*k*(k-1));
				}
			}
		}
		A_i <- 0;
		for (s in 1:k){
			for (l in 1:m){
				YY_i[k+3] <- YY_i[k+3] + (Y_i[(s-1)*m+l]- YY_i[s])^2/(k*(m-1));
				A_i <- A_i + (Y_i[(s-1)*m+l] - Ybar[(s-1)*m+l])^2/(m*m*k);
			}
		}
		B_i <- 0;
		for (s in 1:k){
			for (l in 1:(m-1)){
				for (r in (l+1):m) B_i <- B_i + 2*(Y_i[(s-1)*m+l] - Ybar[(s-1)*m+l])*(Y_i[(s-1)*m+r] - Ybar[(s-1)*m+r])/(m*m*k);
			}
		}
		YY_i[k+4] <- YY_i[k+4] + A_i + B_i;
		delta_i <- muk;
		delta_i[k+1] <- sigmabeta+sigmaerror/m+sigmagamma;
		delta_i[k+2] <- sigmaalpha;
		delta_i[k+3] <- sigmaerror;
		delta_i[k+4] <- sigmagamma+sigmaalpha+sigmaerror/m;		
		H_i <- diag(sigmagamma+sigmaalpha+sigmaerror/m,nsigma);
		H_i[k+1, k+1] <- (k*(k-1)*(3*k-2)/(m*m))*sigmaerror^2+k*(k-1)*(3*k-2)*sigmagamma^2+(4*k*k*(k-1)/m)*sigmagamma*sigmaerror+(8*k*(k-1)/m)*sigmabeta*sigmaerror+8*k*(k-1)*sigmabeta*sigmagamma;
		H_i[k+2, k+2] <- m^4*k*(k-1)*(2*k-3)*sigmaalpha^2+m^4*k*(k-1)*(2*k-3)/2*sigmagamma^2+(k*(k-1)*m*m/2)*sigmaerror^2+m^4*k*(k-1)*(2*k-3)*sigmaalpha*sigmagamma+((2+(m-1)^2+2*m*(m-1)*(k-2))*k*(k-1)*m^2/2)*sigmaalpha*sigmaerror+((2+(m-1)^2+2*m*(m-1)*(k-2))*k*(k-1)*m^2/2)*sigmagamma*sigmaerror;
		H_i[k+3, k+3] <- (2*k/(m-1))*sigmaerror^2;
		H_i[k+4, k+4] <- (2*m^2*k^2+k*m*(m-1)*(2*m-3+(k-1)*m*(m-1)/2)+m^2*k^2*(3*m-5)/2)*sigmaalpha^2+(m^2*k*((1-3*k)/2+m+m*k/2)+k*m*(m-1)*(2*m-3))*sigmagamma^2+(m*k*(m-3)/2)*sigmaerror^2+(2*m^2*k*(2-k)+m*k*(m-1)*(m*k+5*m-4))*sigmaalpha*sigmagamma+(m*k*(4-m*k+(m-1)^2*m*k*(m-1)/2))*sigmaalpha*sigmaerror+(4+(m-1)^2-m*k+(m-1)*(m*k+4)/2)*m*k*sigmagamma*sigmaerror;
		F_i <- diag(0,nsigma);
		for (s in 1:k) F_i[s,s] <- 1;
		F_i[k+1,k+1] <- 1;
		F_i[k+1,k+3] <- 1/m;
   		F_i[k+1,k+4] <- 1;
   		F_i[k+2,k+2] <- 1;
   		F_i[k+3,k+3] <- 1;
   		F_i[k+4,k+2] <- 1;
   		F_i[k+4,k+3] <- 1/m;
   		F_i[k+4,k+4] <- 1;
		u2_i <- t(F_i) %*% solve(H_i) %*% (YY_i-delta_i);
		U2 <- U2+u2_i;
		FHF <- FHF+t(F_i) %*% solve(H_i) %*% F_i;
		u2sq <- u2sq + u2_i%*%t(u2_i);
	  }
	  DELTA2 <- solve(FHF,U2);
        crit1 <- max(abs(DELTA2))/max(sigma);
	  sigma <- sigma+DELTA2;	  
	  it <- it+1;
	}
        if ((it>=iter)&(crit1>toe)) warning("Iteration does not converge, soluations may be inaccurate!");        
	var <- solve(FHF) %*% u2sq %*% t(solve(FHF));
	var_sigmabeta <- var[k+1,k+1];
	var_sigmaalpha <- var[k+2,k+2];
	var_sigmaerror <- var[k+3,k+3];
	var_sigmagamma <- var[k+4,k+4];
	cov_betaalpha <- var[k+1,k+2];
	cov_betaerror <- var[k+1,k+3];
	cov_betagamma <- var[k+1,k+4];
	cov_alphaerror <- var[k+2,k+3];
	cov_alphagamma <- var[k+2,k+4];
	cov_errorgamma <- var[k+3,k+4];
####################################################
	rho_intra <- (sigmaalpha + sigmagamma)/(sigmaalpha + sigmagamma + sigmaerror);
	if (rho_intra<=0.999999){
		rho_intra_z <- log((1+rho_intra)/(1-rho_intra))/2;
	}else rho_intra_z <- NA;
	ccc_intra <- (sigmaalpha + sigmagamma)/(sigmaalpha + sigmagamma + sigmaerror);
	if (ccc_intra<=0.999999){
		ccc_intra_z <- log((1+ccc_intra)/(1-ccc_intra))/2;
	}else ccc_intra_z <- NA;
	esquare_intra <- 2*sigmaerror;

	if (tran==1){
		TDI_intra <- qnorm(1-(1-CP_a)/2)*sqrt(esquare_intra);
		TDI_intra_w <- log(2*sigmaerror);
		kp0_intra <- TDI_a_intra;
		if (error=="prop") kp0_intra <- log(TDI_a_intra/100+1);
		cp_intra <- 1-2*(1-pnorm(kp0_intra/sqrt(esquare_intra)));
		if (cp_intra<=0.999999){
			cp_intra_l <- log(cp_intra/(1-cp_intra));
		}else cp_intra_l <- NA;
	}

	#inter
	rho_inter <- sigmaalpha/(sigmaalpha + sigmagamma + sigmaerror/m);
	if (rho_inter<=0.999999){
		rho_inter_z <- log((1+rho_inter)/(1-rho_inter))/2;

	}else rho_inter_z <- NA;
	ccc_inter_1 <- sigmaalpha/(sigmaalpha + sigmagamma + sigmabeta + sigmaerror/m);
	if (ccc_inter_1<=0.999999){
		ccc_inter_1_z <- log((1+ccc_inter_1)/(1-ccc_inter_1))/2;
	}else ccc_inter_1_z <- NA;
	accuracy_inter <- (sigmaalpha + sigmagamma + sigmaerror/m)/(sigmaalpha + sigmagamma + sigmaerror/m + sigmabeta);
	if (accuracy_inter<=0.999999){
		accuracy_inter_l <- log(accuracy_inter/(1-accuracy_inter));
	}else accuracy_inter_l <- NA;

	esquare_inter <- 2*(sigmagamma+sigmabeta+sigmaerror/m);
	if (tran==1){
      	        TDI_inter <- qnorm(1-(1-CP_a)/2)*sqrt(esquare_inter);
		TDI_inter_w <- log(esquare_inter);
		kp0_inter <- TDI_a_inter;
       	if (tolower(error)=="prop") kp0_inter <- log(TDI_a_inter/100+1);
		cp_inter <- 1-2*(1-pnorm(kp0_inter/sqrt(esquare_inter)));
		if (cp_inter<=0.999999){
			cp_inter_l <- log(cp_inter/(1-cp_inter));
		}else cp_inter_l <- NA;
	}

	#total	
	rho_total <- sigmaalpha/(sigmaalpha + sigmagamma + sigmaerror);
	if (rho_total<=0.999999){
		rho_total_z <- log((1+rho_total)/(1-rho_total))/2;
	}else rho_total_z <- NA;
	ccc_total <- sigmaalpha/(sigmaalpha + sigmagamma + sigmaerror + sigmabeta);
	if (ccc_total<=0.999999){
		ccc_total_z <- log((1+ccc_total)/(1-ccc_total))/2;
	}else ccc_total_z <- NA;
	accuracy_total <- (sigmaalpha + sigmagamma + sigmaerror)/(sigmaalpha + sigmagamma + sigmaerror + sigmabeta);
	if (accuracy_total<=0.999999){
		accuracy_total_l <- log(accuracy_total/(1-accuracy_total));
	}else accuracy_total_l <- NA;
	esquare_total <- 2*(sigmaerror+sigmagamma+sigmabeta);
	if (tran==1){
		TDI_total <- qnorm(1-(1-CP_a)/2)*sqrt(esquare_total);
		TDI_total_w <- log(esquare_total);
		kp0_total <- TDI_a_total;
		if (tolower(error)=="prop") kp0_total <- log(TDI_a_total/100+1);
		cp_total <- 1-2*(1-pnorm(kp0_total/sqrt(esquare_total)));
		if (cp_total<=0.999999){
			cp_total_l <- log(cp_total/(1-cp_total));
		}else cp_total_l <- NA;
	}
	
	#se of intra
	rhoc0 <- ccc_intra;
        var_ccc_intra <- ((1-rhoc0)^2*(var_sigmaalpha+var_sigmagamma+2*cov_alphagamma)+rhoc0^2*var_sigmaerror-2*(1-rhoc0)*rhoc0*(cov_alphaerror+cov_errorgamma))/(sigmaalpha + sigmagamma + sigmaerror)^2;
	if (var_ccc_intra<=0) var_ccc_intra <- NA;
	se_ccc_intra <- sqrt(var_ccc_intra);
	ccc_intra_notran_lower <- ccc_intra-qnorm(1-alpha)*se_ccc_intra;
	if (ccc_intra<=0.999999){
		se_ccc_intra_z <- se_ccc_intra/(1-ccc_intra^2);
	}else se_ccc_intra_z <- NA;
        ccc_intra_z_low <- ccc_intra_z-qnorm(1-alpha)*se_ccc_intra_z;
        ccc_intra_tran_lower <- (exp(2*ccc_intra_z_low)-1)/(exp(2*ccc_intra_z_low)+1);

	se_rho_intra <- se_ccc_intra;
	se_rho_intra_z <- se_ccc_intra_z;
	rho_intra_z_low <- ccc_intra_z_low;
	rho_intra_notran_lower <- ccc_intra_notran_lower;
	rho_intra_tran_lower <- ccc_intra_tran_lower;

	var_esqua_intra <- 4*(var_sigmaerror);
	if (var_esqua_intra<=0) var_esqua_intra <- NA;
	se_esqua_intra <- sqrt(var_esqua_intra);

	if (tran==1){
		se_TDI_intra_w <- sqrt(var_esqua_intra)/(2*sigmaerror);
		TDI_intra_w_up <- TDI_intra_w+qnorm(1-alpha)*se_TDI_intra_w;
		TDI_intra_e_up <- sqrt(exp(TDI_intra_w_up));
		TDI_intra_tran_upper <- qnorm(1-(1-CP_a)/2)*TDI_intra_e_up;
		var_notran_cpintra <- (1/(exp(kp0_intra^2/esquare_intra))*(1/kp0_intra+kp0_intra/esquare_intra)^2*var_esqua_intra)/(8*pi*esquare_intra);
		if (var_notran_cpintra<=0) var_notran_cpintra <- NA;
		se_notran_cpintra <- sqrt(var_notran_cpintra);
		if (cp_intra<=0.999999){
			var_tran_cpintra <- var_notran_cpintra/((cp_intra*(1-cp_intra))^2);
		}else var_tran_cpintra <- NA;
		if ((var_tran_cpintra<=0) & (!is.na(var_tran_cpintra))) var_tran_cpintra <- NA;
		se_tran_cpintra <- sqrt(var_tran_cpintra);
		cp_intra_notran_lower <- cp_intra - qnorm(1-alpha)*se_notran_cpintra;
		cp_intra_l_low <- cp_intra_l-qnorm(1-alpha)*se_tran_cpintra;
		cp_intra_tran_lower <- exp(cp_intra_l_low)/(1+exp(cp_intra_l_low));
	}

	#se of inter
	rhoc1 <- ccc_inter_1;
	var_ccc_inter <- ((1-rhoc1)^2*var_sigmaalpha+rhoc1^2*(var_sigmabeta+var_sigmagamma+var_sigmaerror/(m^2)+2*cov_betagamma+2*cov_betaerror/m+2*cov_errorgamma/m)-(2*(1-rhoc1)*rhoc1*(cov_betaalpha+cov_alphagamma+cov_alphaerror/m)))/(sigmaalpha + sigmabeta + sigmagamma + sigmaerror/m)^2;
	if (var_ccc_inter<=0) var_ccc_inter <- NA;
	se_ccc_inter <- sqrt(var_ccc_inter);
	ccc_inter_notran_lower <- ccc_inter_1 - qnorm(1-alpha)*se_ccc_inter;
	if (ccc_inter_1<=0.999999){
		se_ccc_inter_1_z <- se_ccc_inter/(1-(ccc_inter_1)^2);
	}else se_ccc_inter_1_z <- NA;
	ccc_inter_z_low <- ccc_inter_1_z-qnorm(1-alpha)*se_ccc_inter_1_z;
	ccc_inter_tran_lower <- (exp(2*ccc_inter_z_low)-1)/(exp(2*ccc_inter_z_low)+1);

	rhoc2b <- rho_inter;
	var_rho_inter <- ((1-rhoc2b)^2*var_sigmaalpha+rhoc2b^2*(var_sigmagamma+var_sigmaerror/(m^2)+2*cov_errorgamma/m)-(2*(1-rhoc2b)*rhoc2b*(cov_alphagamma+cov_alphaerror/m)))/(sigmaalpha + sigmagamma + sigmaerror/m)^2;
	if (var_rho_inter<=0) var_rho_inter <- NA;
	se_rho_inter <- sqrt(var_rho_inter);
	rho_inter_notran_lower <- rho_inter - qnorm(1-alpha)*se_rho_inter;
	if (rho_inter<=0.999999){
		se_rho_inter_z <- se_rho_inter/(1-rho_inter^2);
	}else se_rho_inter_z <- NA;
	rho_inter_z_low <- rho_inter_z-qnorm(1-alpha)*se_rho_inter_z;
	rho_inter_tran_lower <- (exp(2*rho_inter_z_low)-1)/(exp(2*rho_inter_z_low)+1);

	var_accuracy_inter <- ((1-accuracy_inter)^2*(var_sigmaalpha+var_sigmagamma+var_sigmaerror/(m^2)+2*cov_alphaerror/m+2*cov_alphagamma+2*cov_errorgamma/m)+accuracy_inter^2*var_sigmabeta-2*(1-accuracy_inter)*accuracy_inter*(cov_betaalpha+cov_betaerror/m+cov_betagamma))/(sigmaalpha + sigmabeta + sigmaerror/m+sigmagamma)^2;
        if (var_accuracy_inter<=0) var_accuracy_inter <- NA;
	se_accuracy_inter <- sqrt(var_accuracy_inter);
	accuracy_inter_notran_lower <- accuracy_inter - qnorm(1-alpha)*se_accuracy_inter;
	if (accuracy_inter<=0.999999){
		se_accuracy_inter_l <- se_accuracy_inter/(accuracy_inter*(1-accuracy_inter));
	}else se_accuracy_inter_l <- NA;
	accuracy_inter_l_low <- accuracy_inter_l-qnorm(1-alpha)*se_accuracy_inter_l;
	accuracy_inter_tran_lower <- exp(accuracy_inter_l_low)/(1+exp(accuracy_inter_l_low));

	var_esqua_inter <- 4*(var_sigmabeta+var_sigmagamma+var_sigmaerror/(m^2)+2*cov_betagamma+2*cov_betaerror/m+2*cov_errorgamma/m);
	if (var_esqua_inter<=0) var_esqua_inter <- NA;
	se_esqua_inter <- sqrt(var_esqua_inter);

	if (tran==1){
		se_TDI_inter_w <- sqrt(var_esqua_inter)/(2*sigmaerror/m+2*sigmabeta+2*sigmagamma);
		TDI_inter_w_low <- TDI_inter_w-qnorm(1-alpha)*se_TDI_inter_w;
		TDI_inter_w_up <- TDI_inter_w+qnorm(1-alpha)*se_TDI_inter_w;
		TDI_inter_e_low <- sqrt(exp(TDI_inter_w_low));
		TDI_inter_e_up <- sqrt(exp(TDI_inter_w_up));
		TDI_inter_tran_lower <- qnorm(1-(1-CP_a)/2)*TDI_inter_e_low;
		TDI_inter_tran_upper <- qnorm(1-(1-CP_a)/2)*TDI_inter_e_up;
		var_notran_cpinter <- (1/(exp(kp0_inter^2/esquare_inter))*(1/kp0_inter+kp0_inter/esquare_inter)^2*var_esqua_inter)/(8*pi*esquare_inter);
	  	if (var_notran_cpinter<=0) var_notran_cpinter <- NA;
	  	se_notran_cpinter <- sqrt(var_notran_cpinter);
		cp_inter_notran_lower <- cp_inter- qnorm(1-alpha)*se_notran_cpinter;
		if (cp_inter<=0.999999){
			var_tran_cpinter <- var_notran_cpinter/((cp_inter*(1-cp_inter))^2);
		}else var_tran_cpinter <- NA
		if ((var_tran_cpinter<=0) & (!is.na(var_tran_cpinter))) var_tran_cpinter <- NA;
		se_tran_cpinter <- sqrt(var_tran_cpinter);
		cp_inter_l_low <- cp_inter_l-qnorm(1-alpha)*se_tran_cpinter;
		cp_inter_tran_lower <- exp(cp_inter_l_low)/(1+exp(cp_inter_l_low));
	        RBS_inter <- sigmabeta/(sigmaerror/m+sigmagamma);
	}
	
	#se of total
	rhoc3 <- ccc_total;
	var_ccc_total <- ((1-rhoc3)^2*var_sigmaalpha+rhoc3^2*(var_sigmabeta+var_sigmagamma+var_sigmaerror)-2*(1-rhoc3)*rhoc3*(cov_betaalpha+cov_alphagamma+cov_alphaerror)+2*rhoc3^2*(cov_betaerror+cov_betagamma+cov_errorgamma))/(sigmaalpha + sigmabeta + sigmagamma + sigmaerror)^2;
	if (var_ccc_total<=0) var_ccc_total <- NA;
	se_ccc_total <- sqrt(var_ccc_total);
	ccc_total_notran_lower <- ccc_total - qnorm(1-alpha)*se_ccc_total;
	if (ccc_total<=0.999999){
		se_ccc_total_z <- se_ccc_total/(1-(ccc_total)^2);
	}else se_ccc_total_z <- NA;
        ccc_total_z_low <- ccc_total_z-qnorm(1-alpha)*se_ccc_total_z;
	ccc_total_tran_lower <- (exp(2*ccc_total_z_low)-1)/(exp(2*ccc_total_z_low)+1);
	
	rhoc3b <- rho_total;
	var_rho_total <- ((1-rhoc3b)^2*var_sigmaalpha+rhoc3b^2*(var_sigmagamma+var_sigmaerror+2*cov_errorgamma)-(2*(1-rhoc3b)*rhoc3b*(cov_alphagamma+cov_alphaerror)))/(sigmaalpha + sigmagamma + sigmaerror)^2;
	if (var_rho_total<=0) var_rho_total <- NA;
	se_rho_total <- sqrt(var_rho_total);
	rho_total_notran_lower <- rho_total - qnorm(1-alpha)*se_rho_total;
	if (rho_total<=0.999999){
		se_rho_total_z <- se_rho_total/(1-rho_total^2);
	}else se_rho_total_z <- NA;
	rho_total_z_low <- rho_total_z-qnorm(1-alpha)*se_rho_total_z;
	rho_total_tran_lower <- (exp(2*rho_total_z_low)-1)/(exp(2*rho_total_z_low)+1);

	var_accuracy_total <- ((1-accuracy_total)^2*(var_sigmaalpha+var_sigmagamma+var_sigmaerror+2*cov_alphaerror+2*cov_alphagamma+2*cov_errorgamma)+accuracy_total^2*var_sigmabeta-2*(1-accuracy_total)*accuracy_total*(cov_betaalpha+cov_betaerror+cov_betagamma))/(sigmaalpha + sigmabeta + sigmaerror+sigmagamma)^2;
	if (var_accuracy_total<=0) var_accuracy_total <- NA;
	se_accuracy_total <- sqrt(var_accuracy_total);
	accuracy_total_notran_lower <- accuracy_total - qnorm(1-alpha)*se_accuracy_total;
	if (accuracy_total<=0.999999){
		se_accuracy_total_l <- se_accuracy_total/(accuracy_total*(1-accuracy_total));
	}else se_accuracy_total_l <- NA;
	accuracy_total_l_low <- accuracy_total_l-qnorm(1-alpha)*se_accuracy_total_l;
	accuracy_total_tran_lower <- exp(accuracy_total_l_low)/(1+exp(accuracy_total_l_low));

	var_esqua_total <- 4*(var_sigmabeta+var_sigmagamma+var_sigmaerror+2*cov_betagamma+2*cov_betaerror+2*cov_errorgamma);
	if (var_esqua_total<=0) var_esqua_total <- NA;
	se_esqua_total <- sqrt(var_esqua_total);

	if (tran==1){
		se_TDI_total_w <- sqrt(var_esqua_total)/(2*sigmaerror+2*sigmabeta+2*sigmagamma);
		TDI_total_w_up <- TDI_total_w+qnorm(1-alpha)*se_TDI_total_w;
		TDI_total_e_up <- sqrt(exp(TDI_total_w_up));
		TDI_total_tran_upper <- qnorm(1-(1-CP_a)/2)*TDI_total_e_up;
		var_notran_cptotal <- (1/(exp(kp0_total^2/esquare_total))*(1/kp0_total+kp0_total/esquare_total)^2*var_esqua_total)/(8*pi*esquare_total);
		if (var_notran_cptotal<=0) var_notran_cptotal <- NA;
		se_notran_cptotal <- sqrt(var_notran_cptotal);
		if (cp_total<=0.999999){
			var_tran_cptotal <- var_notran_cptotal/((cp_total*(1-cp_total))^2);
		}else var_tran_cptotal <- NA;
		if ((var_tran_cptotal<=0) & (!is.na(var_tran_cptotal))) var_tran_cptotal <- NA;
		se_tran_cptotal <- sqrt(var_tran_cptotal);
		cp_total_notran_lower <- cp_total - qnorm(1-alpha)*se_notran_cptotal;
		cp_total_l_low <- cp_total_l-qnorm(1-alpha)*se_tran_cptotal;
		cp_total_tran_lower <- exp(cp_total_l_low)/(1+exp(cp_total_l_low));
		RBS_total <- sigmabeta/(sigmaerror+sigmagamma);
	}
	if (tran==0){
		ccc_intra_tran_lower <- ccc_intra_notran_lower;
		ccc_inter_tran_lower <- ccc_inter_notran_lower;
		ccc_total_tran_lower <- ccc_total_notran_lower;
		rho_intra_tran_lower <- rho_intra_notran_lower;
		rho_inter_tran_lower <- rho_inter_notran_lower;
		rho_total_tran_lower <- rho_total_notran_lower;
		accuracy_inter_tran_lower <- accuracy_inter_notran_lower;
		accuracy_total_tran_lower <- accuracy_total_notran_lower;
	}
	if (tran==1){
		if (tolower(error)=="prop"){
			TDI_intra <- 100*(exp(TDI_intra)-1);
			TDI_inter <- 100*(exp(TDI_inter)-1);
			TDI_total <- 100*(exp(TDI_total)-1);
			TDI_intra_tran_upper <- 100*(exp(TDI_intra_tran_upper)-1);
			TDI_inter_tran_upper <- 100*(exp(TDI_inter_tran_upper)-1);
			TDI_total_tran_upper <- 100*(exp(TDI_total_tran_upper)-1);
		}
	}
	TDI_intra <- round(TDI_intra,dec);
	TDI_inter <- round(TDI_inter,dec);
	TDI_total <- round(TDI_total,dec);
	TDI_intra_tran_upper <- round(TDI_intra_tran_upper,dec);
	TDI_inter_tran_upper <- round(TDI_inter_tran_upper,dec);
	TDI_total_tran_upper <- round(TDI_total_tran_upper,dec);

	Intra <- cbind(rbind(stat1,stat2,stat3),round(rbind(cbind(ccc_intra,rho_intra,accuracy_intra,TDI_intra,cp_intra,RBS_intra),cbind(ccc_intra_tran_lower, rho_intra_tran_lower, accuracy_intra_tran_lower, TDI_intra_tran_upper, cp_intra_tran_lower, NA),	cbind(CCC_a_intra,CCC_a_intra,NA,TDI_a_intra, CP_a,NA)),4));
	Inter <- cbind(rbind(stat1,stat2,stat3),round(rbind(cbind(ccc_inter_1,rho_inter,accuracy_inter,TDI_inter,cp_inter,RBS_intra),cbind(ccc_inter_tran_lower, rho_inter_tran_lower, accuracy_inter_tran_lower, TDI_inter_tran_upper, cp_inter_tran_lower, NA),cbind(CCC_a_inter,NA,NA,TDI_a_inter, CP_a,NA)),4));
	Total <- cbind(rbind(stat1,stat2,stat3),round(rbind(cbind(ccc_total,rho_total,accuracy_total,TDI_total,cp_total,RBS_total),cbind(ccc_total_tran_lower, rho_total_tran_lower, accuracy_total_tran_lower, TDI_total_tran_upper, cp_total_tran_lower, NA),	cbind(CCC_a_total,NA,NA,TDI_a_total, CP_a,NA)),4));
	estimates <- list(Intra=Intra, Inter=Inter, Total=Total, Inputs=Inputs);
	class(estimates) <- "unified_agreement";
}
return(estimates)
}

