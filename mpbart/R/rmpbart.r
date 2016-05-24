#' Multinomial Probit Bayesian Additive Regression Trees
#'
#' A function to implement multinomial probit regression via Bayesian Addition Regression Trees using partial marginal data augmentation.
#'
#' @param x.train Training data predictors.
#' @param y.train Training data observed classes.
#' @param x.test Test data predictors.
#' @param Prior List of Priors for MPBART: e.g., Prior = list(nu=p+2,  V= diag(p - 1), ntrees=200,  kfac=2.0,  pbd=1.0, pb=0.5 , beta = 2.0, alpha = 0.95, nc = 100, priorindep = 0,  minobsnode = 10)
#' @param Mcmc List of MCMC starting values, burn-in ...: e.g.,     list(sigma0 = diag(p - 1), keep = 1, burn = 100, ndraws = 1000, keep_sigma_draws=FALSE)
#' @param seedvalue random seed value: e.g., seedvalue = 99
#' @useDynLib mpbart
#' @examples
#' 
#' set.seed(64)
#' library(mpbart)
#' p=3
#' train_wave = mlbench.waveform(50)
#' test_wave = mlbench.waveform(100)
#' traindata = data.frame(train_wave$x, y = train_wave$classes) 
#' testdata = data.frame(test_wave$x, y = test_wave$classes)
#' 
#' x.train = data.frame(train_wave$x)
#' x.test = data.frame(test_wave$x)
#' 
#' y.train = train_wave$classes
#' 
#' sigma0 = diag(p-1)
#' burn = 100
#' ndraws = 200 # a higher number >=1000 is more appropriate.
#' 
#' Mcmc1=list(sigma0=sigma0, burn = burn, ndraws = ndraws)
#' Prior1 = list(nu=p+2,
#'               V=(p+2)*diag(p-1),
#'               ntrees = 5, #typically 200 trees is good 
#'               kfac = 2.0, 
#'               pbd = 1.0, 
#'               pb = 0.5, 
#'               alpha = 0.99,  
#'               beta =  2.0, 
#'               nc = 200, 
#'               priorindep = FALSE)
#' 
#' 
#' 
#' out = rmpbart(x.train = x.train, y.train = y.train, x.test = x.test, 
#'             Prior = Prior1, Mcmc=Mcmc1, seedvalue = 99)
#' 
#' #confusion matrix train
#' table(y.train, out$predicted_class_train)
#' table(y.train==out$predicted_class_train)/sum(table(y.train==out$predicted_class_train))
#'   
#' 
#' #confusion matrix test
#' table(test_wave$classes, out$predicted_class_test)
#' 
#' test_err <- sum(test_wave$classes != out$predicted_class_test)/
#'     sum(table(test_wave$classes == out$predicted_class_test))
#' 
#' cat("test error :", test_err )
#' @export
rmpbart =
  function(x.train, y.train, x.test = NULL, Prior = NULL, Mcmc = NULL, seedvalue = NULL) 
  {
    
if(is.null(seedvalue)){
	seedvalue = 99
} else {
	set.seed(seedvalue)
}

XEx = NULL;
for(i in 1:nrow(x.train)){
XEx = rbind(XEx, matrix(rep(x.train[i,], p-1), byrow = TRUE, ncol = ncol(x.train) ) )
}


if(!is.na(x.test)[1]){
	testXEx = NULL;
	for(i in 1:nrow(x.test)){
	testXEx = rbind(testXEx, matrix(rep(x.test[i,], p-1), byrow = TRUE, ncol = ncol(x.test) ) )
	}
} else {
	testXEx = 0
}
p = length(unique(y.train))

Data = list(p=p,y=y.train,X= XEx)

testData = list(p=p,X= testXEx)

p=Data$p
y=Data$y
X=Data$X
testX = testData$X
   
    
    levely=as.numeric(levels(as.factor(y)))
    
    cat("Table of y values",fill=TRUE)
    print(table(y))
    
    n=length(y)
    k=ncol(X)
    pm1=p-1
    testn = nrow(testX)/(p-1)
    
    if(missing(Prior)) 
    {nu=pm1+3; V=nu*diag(pm1);
      ntrees=200; kfac=2.0;pbd=1.0;pb=0.5;beta = 2.0;alpha = 0.95; nc = 100; priorindep = 0; minobsnode = 10;
    }
    else 
    {if(is.null(Prior$nu)) {nu=pm1+3} else {nu=Prior$nu}
     if(is.null(Prior$V)) {V=nu*diag(pm1)} else {V=Prior$V}
     if(is.null(Prior$ntrees)) {ntrees=200} else {ntrees=Prior$ntrees}
     if(is.null(Prior$kfac)) {kfac=2.0} else {kfac=Prior$kfac}
     if(is.null(Prior$pbd)) {pbd=1.0} else {pbd=Prior$pbd}
     if(is.null(Prior$pb)) {pb=0.5} else {pb=Prior$pb}
     if(is.null(Prior$beta)) {beta = 2.0} else {beta=Prior$beta}
     if(is.null(Prior$alpha)) {alpha = 0.95} else {alpha=Prior$alpha}
     if(is.null(Prior$nc)) {nc=100} else {nc=Prior$nc}
     if(is.null(Prior$priorindep)) {priorindep= FALSE} else {priorindep=Prior$priorindep}
     if(is.null(Prior$minobsnode)) {minobsnode= 10} else {minobsnode=Prior$minobsnode}
     
     
    }

        if(is.null(Mcmc$sigma0)) {sigma0=diag(pm1)} else {sigma0=Mcmc$sigma0}
    
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    if(is.null(Mcmc$burn)) {burn=100} else {burn=Mcmc$burn}
    if(is.null(Mcmc$ndraws)) {ndraws=1000} else {ndraws=Mcmc$ndraws}
    if(is.null(Mcmc$keep_sigma_draws)) {keep_sigma_draws=FALSE} else {keep_sigma_draws=Mcmc$keep_sigma_draws}
    
    
    
    

	
    C=chol(solve(sigma0))
    #
    #  C is upper triangular root of sigma^-1 (G) = C'C
    #
    sigmai=crossprod(C)
    
   
	if( (priorindep ==TRUE) || (keep_sigma_draws==FALSE)){
	sigmasample = as.double(0);
	savesigma = 0;
	} else {
	sigmasample = as.double(rep(sigma0, ndraws+burn));
	savesigma = 1;
	}
  
	
    res =   .C("rmnpMDA",w=as.double(rep(0,nrow(X))),
               trainx= as.double(t(X)), 
               testx= as.double(t(testX)),
               mu = as.double(rep(0,nrow(X))),
               sigmai = as.double(sigmai),
               V = as.double(V),
               n = as.integer(length(y)),
               n_dim = as.integer(ncol(sigmai)),
               y = as.integer(y), 
               n_cov = as.integer(k), 
               nu = as.integer(nu), 
               trainpred = as.double(rep(0,p*n)) , 
               testn = as.integer(testn), 
               testpred = as.double(rep(0,p*testn)), 
               ndraws = as.integer(ndraws), 
               burn = as.integer(burn),
               ntrees = as.integer(ntrees),
               kfac = as.double(kfac), 
               pbd = as.double(pbd), 
               pb = as.double(pb), 
               alpha = as.double(alpha),  
               beta =  as.double(beta),
               nc = as.integer(nc),
				 savesigma = as.integer(savesigma),
				 minobsnode = as.integer(minobsnode),
               sigmasample = sigmasample,
			   PACKAGE="mpbart")      

class_prob_train = matrix(res$trainpred,ncol = p, byrow = TRUE)
predicted_class_train = apply(class_prob_train,1,which.max)

class_prob_test = matrix(res$testpred,ncol = p, byrow = TRUE)
predicted_class_test = apply(class_prob_test,1,which.max)

ret = list(class_prob_train = class_prob_train, 
			predicted_class_train = predicted_class_train,
			class_prob_test = class_prob_test, 
			predicted_class_test = predicted_class_test);
			
			
class(ret) = "mpbart"

return(ret)
}