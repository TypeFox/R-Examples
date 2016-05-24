IIR <-
function(data,k,m,testindex,refindex,error,alpha=0.05){
Y <- data;
n <- dim(data)[1];
obs <- k*m*n;

if (error=="prop") Y <- log(Y);
Ydot <- matrix(0,n,k);

U <- matrix(0,n,k);

for (i in 1:n){
	for (s in 1:k){
		Ydot[i,s] <- sum(Y[i,((s-1)*m+1):(s*m)])/m;
		U[i,s] <- sum((Y[i,((s-1)*m+1):(s*m)]-Ydot[i,s])^2)/(m-1);
	}
}

Ydot2 <- Ydot^2;
Z <- matrix(0,n,k*(k-1)/2);
p <- 1;
for (s in 1:(k-1)){
	for (t in (s+1):k){
		Z[,p] <- Ydot[,s]*Ydot[,t];
		p <- p+1;
	}
}

#calculate estimates of mu_star, sigma2_star, delta2 and rho_mu
mustar <- apply(Ydot,2,sum)/n;
sigma2star <- apply(U,2,sum)/n;
delta2 <- apply(Ydot2,2,sum)/n-mustar^2-sigma2star/m;
rhomu <- matrix(0,k*(k-1)/2,1);

tempmu1 <- matrix(0,k,k);
tempmu2 <- matrix(0,k,k);
tempdelta <- matrix(0,k,k);
tempmudelta1 <- matrix(0,k,k);
tempmudelta2 <- matrix(0,k,k);
tempsigmadelta <- matrix(0,k,k);

for (s in 1:(k-1)){
	for (t in (s+1):k){
		tempmu1[s,t] <- mustar[s]*mustar[t];
		tempmu2[s,t] <- (mustar[s]-mustar[t])^2;
		tempdelta[s,t] <- sqrt(delta2[s]*delta2[t]);
		tempmudelta1[s,t] <- mustar[s]^2*(delta2[t]+sigma2star[t]/m);
		tempmudelta2[s,t] <- mustar[t]^2*(delta2[s]+sigma2star[s]/m);
		tempsigmadelta[s,t] <- (delta2[s]+sigma2star[s])*(delta2[t]+sigma2star[t]);
	}
}

muproc <- matrix(0,k*(k-1)/2,1);
deltaproc <- matrix(0,k*(k-1)/2,1);
musum <- matrix(0,k*(k-1)/2,1);
mudeltaproc1 <- matrix(0,k*(k-1)/2,1);
mudeltaproc2 <- matrix(0,k*(k-1)/2,1);
sigmadeltaproc <- matrix(0,k*(k-1)/2,1);

p <- 1;
for (s in 1:(k-1)){
	for (t in (s+1):k){
		muproc[p,1] <- mustar[s]*mustar[t];
		deltaproc[p,1] <- sqrt(delta2[s]*delta2[t]);
		musum[p,1] <- (mustar[s]-mustar[t])^2;
		mudeltaproc1[p,1] <- mustar[s]^2*(delta2[t]+sigma2star[t]/m);
		mudeltaproc2[p,1] <- mustar[t]^2*(delta2[s]+sigma2star[s]/m);
		sigmadeltaproc[p,1] <- (delta2[s]+sigma2star[s])*(delta2[t]+sigma2star[t]);
		p <- p+1;
	}
}

for (s in 1:(k*(k-1)/2)) rhomu[s,1] <- (sum(Z[,s])/n-muproc[s,1])/deltaproc[s,1];


#The following calculates the variance matrix H_i for YY_i
Di1 <- diag(k);
vi1 <- diag(k);
diag(vi1) <- apply((Ydot-matrix(colMeans(Ydot),n,k,byrow=TRUE))^2,2,sum)/(n-1);

Di2 <- diag(k);
vi2 <- diag(k);
diag(vi2) <- apply((U-matrix(colMeans(U),n,k,byrow=TRUE))^2,2,sum)/(n-1);

Di3 <- diag(k);
vi3 <- diag(k);
diag(vi3) <- 2*(delta2+sigma2star/m)^2+4*mustar^2*(delta2+sigma2star/m);

Di4 <- diag(k*(k-1)/2);
vi4 <- diag(k*(k-1)/2);
diag(Di4) <- deltaproc;
diag(vi4) <- (rhomu*deltaproc)^2+2*muproc*rhomu*deltaproc+sigmadeltaproc+mudeltaproc1+mudeltaproc2;

#The following calculates the derivative of the function w.r.t each parameter
Gi1 <- matrix(k);
Gi2 <- diag(2*mustar);
Gi3 <- diag(k)/m;

index <- matrix(0,2,k*(k-1)/2);
p <- 1;
for (s in 1:(k-1)){
	for (t in (s+1):k){
		index[1,p] <- s;
		index[2,p] <- t;
		p <- p+1;
	}
}

Gi4 <- matrix(0,k,k*(k-1)/2);
for (s in 1:(k*(k-1)/2)){
	Gi4[index[1,s],s] <- mustar[index[2,s]];
	Gi4[index[2,s],s] <- mustar[index[1,s]];
}

Gi6 <- matrix(0,k,k*(k-1)/2);
for (s in 1:(k*(k-1)/2)){
	Gi6[index[1,s],s] <- rhomu[s]*sqrt(delta2[index[2,s]])/(2*sqrt(delta2[index[1,s]]));
	Gi6[index[2,s],s] <- rhomu[s]*sqrt(delta2[index[1,s]])/(2*sqrt(delta2[index[2,s]]));
}

#the following calculates the psi matrix
Psi <- matrix(0,3*k+k*(k-1)/2,3*k+k*(k-1)/2);
Psi[1:k,1:k] <- n * Di1 %*% solve(vi1) %*% t(Di1);
Psi[(k+1):(2*k),(k+1):(2*k)] <- n*Di2 %*% solve(vi2) %*% t(Di2);
Psi[(2*k+1):(3*k),1:k] <- n*Di3 %*% solve(vi3) %*% t(Gi2);
Psi[(2*k+1):(3*k),(k+1):(2*k)] <- n*Di3 %*% solve(vi3) %*% t(Gi3);
Psi[(2*k+1):(3*k),(2*k+1):(3*k)] <- n*Di3 %*% solve(vi3) %*% t(Di3);
Psi[(3*k+1):(3*k+k*(k-1)/2),1:k] <- n*Di4 %*% solve(vi4) %*% t(Gi4);
Psi[(3*k+1):(3*k+k*(k-1)/2),(2*k+1):(3*k)] <- n*Di4 %*% solve(vi4) %*% t(Gi6);
Psi[(3*k+1):(3*k+k*(k-1)/2),(3*k+1):(3*k+k*(k-1)/2)] <- n*Di4 %*% solve(vi4) %*% t(Di4);

#the following calculates the Lambda matrix
Lambda <- matrix(0,3*k+k*(k-1)/2,3*k+k*(k-1)/2);
varY <- matrix(0,k,k);
varY2 <- matrix(0,k,k);
varU <- matrix(0,k,k);
covYU <- matrix(0,k,k);
covYY2 <- matrix(0,k,k);
covUY2 <- matrix(0,k,k);
varZ <- matrix(0,k*(k-1)/2,k*(k-1)/2);
covYZ <- matrix(0,k,k*(k-1)/2);
covUZ <- matrix(0,k,k*(k-1)/2);
covY2Z <- matrix(0,k,k*(k-1)/2);

for (i in 1:n){
	varY <- varY+(Ydot[i,]-mustar)%*%t(Ydot[i,]-mustar);
	varU <- varU+(U[i,]-sigma2star)%*%t(U[i,]-sigma2star);
	varY2 <- varY2+(Ydot2[i,]-delta2-sigma2star/m-mustar^2)%*%t(Ydot2[i,]-delta2-sigma2star/m-mustar^2);
	covYU <- covYU+(Ydot[i,]-mustar)%*%t(U[i,]-sigma2star);
	covYY2 <- covYY2+(Ydot[i,]-mustar)%*%t(Ydot2[i,]-delta2-sigma2star/m-mustar^2);
	covUY2 <- covUY2+(U[i,]-sigma2star)%*%t(Ydot2[i,]-delta2-sigma2star/m-mustar^2);
	varZ <- varZ+(Z[i,]-rhomu*deltaproc-muproc)%*%t(Z[i,]-rhomu*deltaproc-muproc);
	covYZ <- covYZ+(Ydot[i,]-mustar)%*%t(Z[i,]-rhomu*deltaproc-muproc);
	covUZ <- covUZ+(U[i,]-sigma2star)%*%t(Z[i,]-rhomu*deltaproc-muproc);
	covY2Z <- covY2Z+(Ydot2[i,]-delta2-sigma2star/m-mustar^2)%*%t(Z[i,]-rhomu*deltaproc-muproc);
}

Lambda[1:k,1:k] <- Di1 %*% solve(vi1) %*% varY %*% solve(vi1) %*% t(Di1);
Lambda[1:k,(k+1):(2*k)] <- Di1 %*% solve(vi1) %*% covYU %*% solve(vi2) %*% t(Di2);
Lambda[(k+1):(2*k),1:k] <- t(Lambda[1:k,(k+1):(2*k)]);
Lambda[1:k,(2*k+1):(3*k)] <- Di1 %*% solve(vi1) %*% covYY2 %*% solve(vi3) %*% t(Di3);
Lambda[(2*k+1):(3*k),1:k] <- t(Lambda[1:k,(2*k+1):(3*k)]);
Lambda[1:k,(3*k+1):(3*k+k*(k-1)/2)] <- Di1 %*% solve(vi1) %*% covYZ %*% solve(vi4) %*% t(Di4);
Lambda[(3*k+1):(3*k+k*(k-1)/2),1:k] <- t(Lambda[1:k,(3*k+1):(3*k+k*(k-1)/2)]);

Lambda[(k+1):(2*k),(k+1):(2*k)] <- Di2 %*% solve(vi2) %*% varU %*% solve(vi2) %*% t(Di2);
Lambda[(k+1):(2*k),(2*k+1):(3*k)] <- Di2 %*% solve(vi2) %*% covUY2 %*% solve(vi3) %*% t(Di3);
Lambda[(2*k+1):(3*k),(k+1):(2*k)] <- t(Lambda[(k+1):(2*k),(2*k+1):(3*k)]);
Lambda[(k+1):(2*k),(3*k+1):(3*k+k*(k-1)/2)] <- Di2 %*% solve(vi2) %*% covUZ %*% solve(vi4) %*% t(Di4);
Lambda[(3*k+1):(3*k+k*(k-1)/2),(k+1):(2*k)] <- t(Lambda[(k+1):(2*k),(3*k+1):(3*k+k*(k-1)/2)]);

Lambda[(2*k+1):(3*k),(2*k+1):(3*k)] <- Di3 %*% solve(vi3) %*% varY2 %*% solve(vi3) %*% t(Di3);
Lambda[(2*k+1):(3*k),(3*k+1):(3*k+k*(k-1)/2)] <- Di3 %*% solve(vi3) %*% covY2Z %*% solve(vi4) %*% t(Di4);
Lambda[(3*k+1):(3*k+k*(k-1)/2),(2*k+1):(3*k)] <- t(Lambda[(2*k+1):(3*k),(3*k+1):(3*k+k*(k-1)/2)]);
Lambda[(3*k+1):(3*k+k*(k-1)/2),(3*k+1):(3*k+k*(k-1)/2)] <- Di4 %*% solve(vi4) %*% varZ %*% solve(vi4) %*% t(Di4);

COV <- solve(Psi) %*% Lambda %*% t(solve(Psi));#estimates of variance covariance matrix for all parameters;

ntest <- length(testindex);
nref <- length(refindex);
msdw <- 2*sigma2star;
withintest <- msdw[testindex];
withinref <- msdw[refindex];
msdtest <- log(sum(withintest)/(ntest));
msdref <- log(sum(withinref)/(nref));
IIR <- exp(msdtest-msdref);

Dvector <- rep(0,k);

sumtestsigma2 <- sum(withintest)/2;
sumrefsigma2 <- sum(withinref)/2;

#the following calculates the derivative with respect to sigma2star
Dvector[testindex] <- nref/(ntest*sumrefsigma2*IIR);
Dvector[refindex] <- -nref*sumtestsigma2/(ntest*sumrefsigma2^2*IIR);

COV_sigma2 <- COV[(k+1):(2*k),(k+1):(2*k)];

var_logIIR <- t(Dvector) %*% COV_sigma2 %*% Dvector; #Variance for log transformed IIR;
se_var_logIIR <- sqrt(var_logIIR); #standard deviation for log transformed TIR;
logIIR_lower_limit <- log(IIR)-qnorm(1-alpha/2)*se_var_logIIR; #lower limit for log transformed IIR;
logIIR_upper_limit <- log(IIR)+qnorm(1-alpha/2)*se_var_logIIR; #upper limit for log transformed IIR;
IIR_lower_limit <- exp(logIIR_lower_limit); #lower limit for IIR;
IIR_upper_limit <- exp(logIIR_upper_limit); #upper limit for IIR;

return(list(IIR=IIR,IIR_upper=IIR_upper_limit,IIR_lower=IIR_lower_limit));
}

