## Itterative Prediction for gaussian process models
##
## May 18th, 2012

predict.GP <- function(object,xnew = object$X,M=1,...){
if (is.GP(object) == FALSE){
	stop("The object in question is not of class \"GP\" \n")
}
if (is.matrix(xnew) == FALSE){
	xnew = as.matrix(xnew)
}
if (M <= 0){
	M = 1;
	warning("M was assigned non-positive, changed to 1. \n")
}

X = object$X
Y = object$Y
n = nrow(X)
d = ncol(X)
corr = object$correlation_param
power = corr$power
nu = corr$nu

if (d != ncol(xnew)){
	stop("The training and prediction data sets are of 
	different dimensions. \n")
}

beta = object$beta;
sig2 = object$sig2;
delta = object$delta;

dim(beta) = c(d,1);
R = corr_matrix(X,beta,corr);

One = rep(1,n);
LO = diag(n);
Sig = R + delta*LO;
L = chol(Sig);

## adding in the check on delta to see about the itterative approach
if (delta == 0){
	Sig_invOne = solve(L,solve(t(L),One));
	Sig_invY = solve(L,solve(t(L),Y));
	Sig_invLp = solve(L,solve(t(L),LO));
} else {
## Adding in the itterative approach section
	s_Onei = One;
	s_Yi = Y;
	s_Li = LO;
	Sig_invOne = matrix(0, ncol = 1, nrow = n);
	Sig_invY = matrix(0, ncol = 1, nrow = n);
	Sig_invLp = matrix(0, ncol = n, nrow = n)
	for (it in 1:M)
	{
		s_Onei = solve(L,solve(t(L),delta*s_Onei));
		Sig_invOne = Sig_invOne + s_Onei/delta;

		s_Yi = solve(L,solve(t(L),delta*s_Yi));
		Sig_invY = Sig_invY + s_Yi/delta;

		s_Li = solve(L,solve(t(L),delta*s_Li));
		Sig_invLp = Sig_invLp + s_Li/delta;
	}
}

nnew = nrow(xnew);
Y_hat = rep(0,nnew);
MSE = rep(0,nnew);
## prediction is different between exponential and matern correlation functions
if (corr$type == "exponential"){
	for(kk in 1:nnew)
	{
		## Changing to accomadate the itterative approach
		xn = matrix(xnew[kk,],nrow=1);
		r=exp(-(abs(X-as.matrix(rep(1,n))%*%(xn))^power)%*%(10^beta));
		yhat = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invY;
		Y_hat[kk] = yhat;
	
		## Adding itterative steps
		if (delta == 0)
		{
			Sig_invr = solve(L,solve(t(L),r));
		} else {
			## if delta != 0, start itterations
			s_ri = r;
			Sig_invr = matrix(0, ncol = 1, nrow = n);
			for (it in 1:M)
			{
				s_ri = solve(L,solve(t(L),delta*s_ri));
				Sig_invr = Sig_invr + s_ri/delta;
			}
		}
		cp_delta_r = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invr;
	
		cp_delta_Lp = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invLp;
		mse = sig2*(1-2*cp_delta_r+cp_delta_Lp%*%R%*%t(cp_delta_Lp));
		MSE[kk] = mse*(mse>0);
	}
}
if (corr$type == "matern"){
	for(kk in 1:nnew)
	{
		## Changing to accomadate the itterative approach
		xn = matrix(xnew[kk,],nrow=1);
		temp = 10^beta
		temp = matrix(temp,ncol=d,nrow=(length(X)/d),byrow=TRUE)
		temp = 2*sqrt(nu)*abs(X-as.matrix(rep(1,n))%*%(xn))*(temp)
		ID = which(temp==0)

		rd=(1/(gamma(nu)*2^(nu-1)))*(temp^nu)*besselK(temp,nu)	
		rd[ID]=1;

		r = matrix(apply(rd,1,prod),ncol=1)		
		yhat = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invY;
		Y_hat[kk] = yhat;
	
		## Adding itterative steps
		if (delta == 0)
		{
			Sig_invr = solve(L,solve(t(L),r));
		} else {
			## if delta != 0, start itterations
			s_ri = r;
			Sig_invr = matrix(0, ncol = 1, nrow = n);
			for (it in 1:M)
			{
				s_ri = solve(L,solve(t(L),delta*s_ri));
				Sig_invr = Sig_invr + s_ri/delta;
			}
		}
		cp_delta_r = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invr;
	
		cp_delta_Lp = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invLp;
		mse = sig2*(1-2*cp_delta_r+cp_delta_Lp%*%R%*%t(cp_delta_Lp));
		MSE[kk] = mse*(mse>0);
	}
}

prediction = NULL

names = c()
for (i in 1:d)
{
	names[i] = paste("xnew.",i, sep="")
}
names[d+1] = "Y_hat"
names[d+2] = "MSE"

full_pred = cbind(xnew, Y_hat, MSE)
colnames(full_pred) = names
prediction$Y_hat = Y_hat
prediction$MSE = MSE
prediction$complete_data = full_pred

return(prediction)
}
