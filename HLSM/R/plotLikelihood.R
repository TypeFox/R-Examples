plotLikelihood = function(object,burnin = 0, thin = 1){
    if(class(object) != 'HLSM')(stop("object must be of class HLSM"))
    xx = getLikelihood(object, burnin = burnin, thin = thin)
    plot(xx,type='l', main = 'Log-likelihood from MCMC draws',ylab = 'estimated log-likelihood')
}

