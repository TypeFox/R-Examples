## Optimizes the deviance using multiple optim starts
## 
## May 8th, 2012

GP_fit <- function(X,Y,control=c(200*d,80*d,2*d),nug_thres=20,trace=FALSE,maxit = 100,
	corr = list(type="exponential",power=1.95), optim_start = NULL){
if (is.matrix(X) == FALSE){
	X = as.matrix(X)
}

## Checking to make sure that the X design is scaled between 0 and 1
if ((min(X) <0) | (max(X)>1) | ((max(X) - min(X)) <=0.5)){
	cat("The optimization routine assumes that inputs are in [0,1].\n")
}

n = nrow(X)
d = ncol(X)
## Checking the dimensions between the different inputs
if (n != length(Y)){
	stop("The dimensions of X and Y do not match. \n")
}
if (n != length(Y)){
	stop("The dimension of the training points and simulator values are different. \n")
}
if (nug_thres < 10){
	warning("nug_thres is outside of the normal range of [10, 25].")
}
if (nug_thres > 25){
	warning("nug_thres is outside of the normal range of [10, 25].")
}
## Need to make sure that the control values are relatively higher / lower and need 
## to print a warning if they are not
if (length(control) != 3){
	stop("control is defined incorrectly. Wrong number of arguments. \n")
}
if (control[1] < control[2]){ 
	stop("control is defined incorrectly. Need control[1] >= control[2]. \n")
}
if (control[1] < control[3]){
	stop("control is defined incorrectly. Need control[1] >= control[3]. \n")
}
if (control[2] < control[3]){
	stop("control is defined incorrectly. Need control[2] >= control[3]. \n")
}

## Checking to see if the vigen starting values for the optim runs are in the correct 
## format
if (!is.null(optim_start)){
	if(!is.matrix(optim_start)){
		if(length(optim_start)/d != floor(length(optim_start)/d)){
			stop("The dimension of optim_start does not match the dimension
					of the problem \n")
		}
		optim_start = matrix(optim_start, byrow = TRUE, ncol = d)
	} else if(ncol(optim_start)!= d){
		stop("The dimension of optim_start does not match the dimension
				of the problem \n")
	}
}
 
param_search = control[1];
param_percent = control[2];
param_clust = control[3];
###########################################################
## Using a grid of control[1] points as default, need to find the 
## control[2] corresponding values with the lowest deviance
## starting out based on the defined range for beta
if (corr$type == "exponential"){
	beta_range = c((corr$power-4)-log10(d),log10(5)+corr$power-log10(d));
	param_init_200d = maximinLHS(param_search,d)*(beta_range[2]-beta_range[1])+beta_range[1];
	param_lower = rep((-10),d);
	param_upper = rep(10,d);
}
if (corr$type == "matern"){
	beta_range = c((2)-log10(d),log10(5)+2-log10(d));
	param_init_200d = maximinLHS(param_search,d)*(beta_range[2]-beta_range[1])+beta_range[1];
	param_lower = rep((-10),d);
	param_upper = rep(10,d);
}

##-------------------------------------------------------##
## Need to evaluate the deviance for the control[1] points
deviance1 = NULL;
for (i in 1:nrow(param_init_200d)){
	temp = GP_deviance(param_init_200d[i,],X,Y,nug_thres,corr = corr);
	deviance1 = rbind(deviance1,c(temp,param_init_200d[i,]));
}
## Need to order the initial values based on their deviance
deviance2 = deviance1[order(deviance1[,1]),];

##-------------------------------------------------------##
## Taking the control[2] smallest deviance values
deviance3 = deviance2[1:param_percent,];

##-------------------------------------------------------##
## Going to cluster these control[2] "best" observations into control[3]
## groups using K-means, over 5 iterations and the smallest WSS as 
## the criterion
points_percen = deviance3[,2:(d+1)];
## Taking the best of 5 different runs of K-means
k = list();
for (i in 1:5){
	k[[i]] = kmeans(points_percen,param_clust);
}
wss_sum = c(sum(k[[1]]$withinss),sum(k[[2]]$withinss),sum(k[[3]]$withinss),sum(k[[4]]$withinss),
	sum(k[[5]]$withinss));
ord = order(wss_sum,decreasing = TRUE);
k[[1]]=k[[ord[1]]];

##---------------------------------------------------------##
## Need to set the control[3] cluster centers as starting values

#param_init = k[[1]]$centers;


##-------------------------------------------------------------------##
## Need to set the control[3] best points in cluster as starting values

param_init = matrix(0,param_clust,d)
for(i in 1:param_clust){
  ID = which(k[[1]]$cluster==i)
  fID = which.min(deviance3[ID,1])
  param_init[i,] = deviance3[ID[fID],2:(d+1)]
}






#############################################################
## For the diagonal search: use 3 points along the diagonal, 
## 1 near each of the ends and one near the middle of the range, 
## this is only necessary above 1 dimension
if (d >= 2) {
	param_wrap = matrix(c(0.20*(beta_range[2]-beta_range[1])+beta_range[1],
		0.5*(beta_range[2]-beta_range[1])+beta_range[1],0.8*(beta_range[2]-beta_range[1])+
		beta_range[1]),byrow=TRUE);

	## Need to run optim() on the wrapped function the 3 times in order to find
	## the lowest starting ponit of the 3 values
	dev = NULL;
	for (i in 1:nrow(param_wrap)){
		temp = optim(param_wrap[i],dev_wrapper,X=X,Y=Y,nug_thres=nug_thres,
			corr=corr,method="L-BFGS-B",lower=param_lower,upper=param_upper, 
			control = c(maxit = maxit)); 
		dev = rbind(dev,c(temp$par,temp$value,param_wrap[i]));
	}
	## Take the best of the 3 based on the likelohood value
	dev = dev[order(dev[,2]),];

##---------------------------------------------------------##
	## Combining the 2*d centers from the clusters with the single point from 
	## the diagonal search
	param_init = rbind(rep(dev[1,1],d),param_init);
}

#############################################################
## We now have the control[3]+1 initializing points for the optim search, as well
## as any points that are specified by the user
dev_val = NULL;
param_init = rbind(param_init, optim_start)
for (i in 1:nrow(param_init)){
	temp = optim(param_init[i,],GP_deviance,X=X,Y=Y,nug_thres=nug_thres,
		corr=corr,method="L-BFGS-B",lower=param_lower,upper=param_upper, 
		control = c(maxit = maxit));
	dev_val = rbind(dev_val,c(temp$par,temp$value));
}
#############################################################
## Making a print statement of what the progress of the optimizer
## only if trace == TRUE
if (trace == TRUE){
	optim_result = cbind(param_init,dev_val)
	col_name = NULL
	if(d==1){
		row_name = NULL
	} else {
		row_name = c("Diagonal Search")
	}
	for (i in 1:d){
		col_name = cbind(col_name, paste("Beta", as.character(i), "Start"))
	}
	for (i in 1:d){
		col_name = cbind(col_name, paste("Beta", as.character(i), "Final"))
	}
	col_name = cbind(col_name, "Deviance Value")
	for (i in 1:(param_clust)){
		row_name = cbind(row_name, paste("Start", as.character(i)))
	}
	colnames(optim_result) = col_name;
	rownames(optim_result) = row_name;
	print(optim_result)
}
dev_val = dev_val[order(dev_val[,(d+1)]),];
beta = (dev_val[1,1:d]);

dim(beta) = c(d,1);
R = corr_matrix(X,beta,corr);
temp = eigen(R,symmetric = TRUE, only.values = TRUE);
eig_val = temp$values;
condnum = kappa(R,triangular = TRUE,exact=TRUE);
max_eigval = eig_val[1];
delta = max(c(0,abs(max_eigval)*(condnum-exp(nug_thres))/(condnum*(exp(nug_thres)-1))));

One = rep(1,n);
dim(Y) = c(n,1);

LO = diag(n);
Sig = R + delta*LO;

L = chol(Sig);

Sig_invOne = solve(L,solve(t(L),One));	
Sig_invY = solve(L,solve(t(L),Y));	

mu_hat = solve(t(One)%*%Sig_invOne,t(One)%*%Sig_invY);
Sig_invb = solve(L,solve(t(L),(Y-One%*%mu_hat)));
sig2 = t(Y-One%*%mu_hat)%*%Sig_invb/n;

GP = NULL;
GP$X = X;
GP$Y = Y;
GP$sig2 = as.vector(sig2);
GP$beta = beta;
GP$delta = delta
GP$nugget_threshold_parameter = nug_thres
GP$correlation_param = corr

class(GP) = "GP"
return(GP)
}