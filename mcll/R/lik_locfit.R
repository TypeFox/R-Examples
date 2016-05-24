# Parameter estimation 

lik_locfit <- function(x,  data.t, ch, bb, prior.func, alp) {

    # data.t: transformed data (of the posterior samples of the parameters)
    # ch: lower triangular matrix of the cholesky for the original posterior samples 
    # bb: means of the original posterior samples 
    # prior.func: prior function (user specified) that provies a log prior density 
    # alp: alpha smoothing parameter 
          
    vec.t <- as.numeric(x)  # transformed scale 
    npar <- length(vec.t) # dimension 

    # transfomation back to the original scale 
    vec <- t((ch%*%vec.t) +bb)

    # posterior 
    fit <- locfit(~. , data= as.data.frame(data.t), deg=2, alpha=alp, ev=vec.t, scale=0, kt="prod" )   
    pred <- predict(fit)
    log.post <- log(pred)
        
    # prior 
    log.prior <- prior.func(vec)
        
    # log likelihood 
    log.lik <- log.post - log.prior -log(det(ch))

    return(-log.lik)  
}
