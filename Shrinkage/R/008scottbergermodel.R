
# Corey M. Yanofsky, February 17, 2009

# implements a slightly altered version of the model of 
# James G. Scott and James O. Berger, "An exploration of aspects of Bayesian multiple testing", 
# J. Stat. Plan. Inf., 136, 2133-2162, 2006

get_xbar_and_sample_size_factor<-function(x,y=NULL){
	z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
	x<-z$x;y<-z$y
	
	if(length(y)==0){
		xbar <- rowMeans(x,na.rm = TRUE)
		n1 <- ncol(x) - rowSums(is.na(x))	
		sample_size_factor <- 1/n1
	}
	else{	
		xbar <- rowMeans(x,na.rm = TRUE)
		xbar <- xbar - rowMeans(y,na.rm = TRUE)
		nx <- ncol(x) - rowSums(is.na(x))
		ny <- ncol(y) - rowSums(is.na(y))
		sample_size_factor <- (nx + nx)/(nx*ny)
	}
	
	list(xbar = xbar, sample_size_factor = sample_size_factor)
}


#----
scottbergermodel <- function(xbar, sample_size_factor, logsigma2, logV, logodds)
{
	xbar <- na.omit(xbar)
	sigma2 <- exp(logsigma2)
	V <- exp(logV)
	odds <- exp(logodds)
	p <- odds/(1 + odds)
	logprior <- -2*log(V + sigma2)
	loglike <- sum(log(p*dnorm(xbar, mean = 0, sd = sqrt(sample_size_factor*sigma2)) + (1-p)*dnorm(xbar, mean = 0, sd = sqrt(sample_size_factor*sigma2 + V))))
	logJacobian <- logsigma2 + logV + logodds - 2*log(1 + odds)
		
	logprior + loglike + logJacobian
}

extract_original_parameters <- function(states)
{
	n_chains <- length(states[[1]])
	n_param <- length(states[[1]][[1]])
	output <- vector("list")
	for(i in 1:n_param)
	{
		output[[i]] <- sapply(states, function(x) x[[1]][i])
		if(n_chains > 1)
		{
			for(j in 2:n_chains)
			{
			output[[i]] <- c(output[[i]], sapply(states, function(x) x[[j]][i]))
			}
		}
	}
	output
}

scottbergermodel_posterior_sample <- function(obj, n_iter, burnin, nchains, opt_init_val)
{
	xbar <- obj$xbar
	sample_size_factor <- obj$sample_size_factor
	logtarget <- function(x) scottbergermodel(xbar, sample_size_factor, x[1], x[2], x[3])
	if(missing(opt_init_val)) opt_init_val <- c(log(var(xbar)/20),log(var(xbar)/20),log(9))
	ctr <- 20
	while(is.infinite(logtarget(opt_init_val)) && ctr > 0)
	{
		ctr <- ctr - 1
		opt_init_val <- c(log(var(xbar)/ctr),log(var(xbar)/ctr),log(9))
	}
	optpoint <- optim(par=opt_init_val, logtarget, hessian = TRUE, control = list(fnscale = -1))
	init_dist_chol <- chol(-solve(optpoint$hessian))
	n_param <- 3
	init_state <- vector("list")
	for(i in 1:nchains) 
	{
		temp <- as.vector(optpoint$par + t(t(init_dist_chol)%*%matrix(rnorm(n_param),nrow = n_param)))
		names(temp) <- c("logsigma2", "logV", "logodds")
		init_state[[i]] <- temp
	}
	DEMC_sample(logtarget,burnin + n_iter,init_state)[(burnin + (1:n_iter))]
}
			

#-------------

k.scottbergermodel.sev<-function(xbar, sample_size_factor, states,estimator_value=NULL,scale_factor=NULL,
			       opt="posterior_mean"){opt<-opt[1]
	stopifnot(length(sample_size_factor)==length(xbar))
	if(length(estimator_value)>0){stopifnot(length(estimator_value) == length(xbar))}
	if(length(scale_factor)>0){stopifnot(length(scale_factor) == length(xbar))}	
	states <- extract_original_parameters(states)
	states[[1]] <- exp(states[[1]])
	states[[2]] <- exp(states[[2]])
	states[[3]] <- exp(-states[[3]])
	#-
	k.fun_posterior_mean<-function(i){
		sig2 <- sample_size_factor[i]*states[[1]]
		p0 <- mean((1 + states[[3]]*sqrt(sig2/(sig2+states[[2]]))*exp((xbar[i]^2)*states[[2]]/(2*sig2*(sig2+states[[2]]))))^-1)
		mu_given_DE <- mean(states[[2]]*xbar[i]/(sig2 + states[[2]]))
		c(p0,mu_given_DE)
	}
	k.fun_posterior_expected_squared_error_loss<-function(i){
		sig2 <- sample_size_factor[i]*states[[1]]
		p0 <- (1 + states[[3]]*sqrt(sig2/(sig2+states[[2]]))*exp((xbar[i]^2)*states[[2]]/(2*sig2*(sig2+states[[2]]))))^-1
		mu_given_DE <- states[[2]]*xbar[i]/(sig2 + states[[2]])
		tau2 <- states[[2]]*sig2/(sig2 + states[[2]])
					
		nexpected_loss<- mean(p0*(estimator_value[i]^2) + (1-p0)*(tau2 + (mu_given_DE - estimator_value[i])^2))
		nexpected_loss
	}
	k.fun_scaled_data_posterior_expected_squared_error_loss<-function(i){
		sig2 <- sample_size_factor[i]*states[[1]]
		p0 <- (1 + states[[3]]*sqrt(sig2/(sig2+states[[2]]))*exp((xbar[i]^2)*states[[2]]/(2*sig2*(sig2+states[[2]]))))^-1
		mu_given_DE <- states[[2]]*xbar[i]/(sig2 + states[[2]])
		tau2 <- states[[2]]*sig2/(sig2 + states[[2]])
					
		nnexpected_loss<- mean(p0*(estimator_value[i]^2) + (1-p0)*(tau2*scale_factor[i]^2 + (mu_given_DE*scale_factor[i] - estimator_value[i])^2))

		nnexpected_loss
	}
	k.fun_scaled_data_posterior_mean<-function(i){
		sig2 <- sample_size_factor[i]*states[[1]]
		p0 <- ((1 + states[[3]]*sqrt(sig2/(sig2+states[[2]]))*exp((xbar[i]^2)*states[[2]]/(2*sig2*(sig2+states[[2]]))))^-1)
		mu_given_DE <- (states[[2]]*xbar[i]/(sig2 + states[[2]]))
		pos_mean <- mean(scale_factor[i]*(1 - p0)*mu_given_DE)
		pos_mean
	}
	k.out<-function(fun,nn=1){ # ,...
		zo<-vapply(1:length(xbar),FUN.VALUE=numeric(nn),FUN=function(i){
			zaux<-try(fun(i=i))
			if(is_err(zaux)){zaux<-rep(as.numeric(NA),nn)}
			zaux})
		zo
	}
	
	#---------------
	if(opt%in% c("posterior_mean","posterior.mean","post_mean","post.mean") ){
		zo<-k.out(fun=k.fun_posterior_mean,nn=2)
		p0<-zo[1,];mu_given_DE<-zo[2,]
		zout<-list(p0 = p0, mu_given_DE = mu_given_DE);return(zout)
	}
	else if(opt%in% c("posterior_expected_squared_error_loss","post.exp.sqerr")){
			fun<-k.fun_posterior_expected_squared_error_loss}
	else if(opt%in% c("scaled_data_posterior_mean","sc.post.mean")){
			fun<-k.fun_scaled_data_posterior_mean}
	else if(opt%in% c("scaled_data_posterior_expected_squared_error_loss","sc.post.exp.sqerr")){
			fun<-k.fun_scaled_data_posterior_expected_squared_error_loss}
	else {stop("bad opt")}
	
	zo<-k.out(fun=fun,nn=1)
	zo
}
scottbergermodel_posterior_mean <- function(xbar, sample_size_factor, states){
	k.scottbergermodel.sev(xbar=xbar, sample_size_factor=sample_size_factor, states=states,estimator_value=NULL,scale_factor=NULL,
			       opt="posterior_mean")	
}
scottbergermodel_posterior_expected_squared_error_loss <- function(xbar, sample_size_factor, states,estimator_value){
	k.scottbergermodel.sev(xbar=xbar, sample_size_factor=sample_size_factor, states=states,estimator_value=estimator_value,scale_factor=NULL,
			       opt="posterior_expected_squared_error_loss")
}
scottbergermodel_scaled_data_posterior_mean <- function(xbar, sample_size_factor, scale_factor, states){
	k.scottbergermodel.sev(xbar=xbar, sample_size_factor=sample_size_factor, states=states,estimator_value=NULL,scale_factor=scale_factor,
			       opt="scaled_data_posterior_mean")
}
scottbergermodel_scaled_data_posterior_expected_squared_error_loss <- function(xbar, sample_size_factor,
					scale_factor, states, estimator_value){
	k.scottbergermodel.sev(xbar=xbar, sample_size_factor=sample_size_factor, states=states,estimator_value=estimator_value,scale_factor=scale_factor,
			       opt="scaled_data_posterior_mean")
}
##----------------
