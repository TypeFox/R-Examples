test.equality.mixed=function (y, x = NULL, w = NULL, arb.R = TRUE, arb.sigma = FALSE, 
    lambda = NULL, mu = NULL, sigma = NULL, R = NULL, alpha = NULL, 
    ...) 
{
    if (arb.R == arb.sigma) 
        stop("Change either 'arb.R' or 'arb.sigma'!")
    v = 1
    while (v == 1) {
        if (arb.R) {
            H0 = regmixEM.mixed(y = y, x = x, w = w, arb.R = TRUE, 
                arb.sigma = FALSE, lambda = lambda, sigma = sigma, 
                mu = mu, R = R, alpha = alpha, ...)
            p = nrow(H0$posterior.beta[[1]])
            k = length(H0$lambda)
            if (is.null(w)) 
                alpha = NULL
            else alpha = H0$alpha

#		n.i=sapply(y,length)
#		N=length(y)
#		tmp=apply(H1$posterior.z*n.i,2,sum)
#		common.sig=sum(tmp/N*H1$sigma)



#            H1 = regmixEM.mixed(y = y, x = x, w = w, arb.R = TRUE, 
#                arb.sigma = TRUE, lambda = H0$lambda, sigma = rep(H0$sigma,k), 
#                mu = H0$mu, R = H0$R, alpha = alpha, ...)

            H1 = regmixEM.mixed(y = y, x = x, w = w, arb.R = TRUE, 
                arb.sigma = TRUE, lambda = H0$lambda, sigma = NULL, 
                mu = H0$mu, R = H0$R, alpha = alpha, ...)

            D = 2 * (H1$loglik - H0$loglik)
            df = k - 1
            alpha = 1 - pchisq(D, df = df)
        }
        else {
            H0 = regmixEM.mixed(y = y, x = x, w = w, arb.R = FALSE, 
                arb.sigma = TRUE, lambda = lambda, sigma = sigma, 
                mu = mu, R = R, alpha = alpha, ...)
            p = nrow(H0$posterior.beta[[1]])
            k = length(H0$lambda)
            if (is.null(w)) 
                alpha = NULL
            else alpha = H0$alpha

#		N=length(y)
#		tmp=(apply(H1$posterior.z,2,sum))/N
#		common.R=matrix(0,ncol=p,nrow=p)
#		for(i in 1:length(H1$lambda)){
#			common.R=common.R+tmp[i]*H1$R[[i]]
#		}


#            H1 = regmixEM.mixed(y = y, x = x, w = w, arb.R = TRUE, 
#                arb.sigma = TRUE, lambda = H0$lambda, sigma = H0$sigma, 
#                mu = H0$mu, R = lapply(1:k,function(i) H0$R), alpha = alpha, ...)
            H1 = regmixEM.mixed(y = y, x = x, w = w, arb.R = TRUE, 
                arb.sigma = TRUE, lambda = H0$lambda, sigma = H0$sigma, 
                mu = H0$mu, R = NULL, alpha = alpha, ...)


            D = 2 * (H1$loglik - H0$loglik)
            df = p * (p + 1) * (k - 1)/2
            alpha = 1 - pchisq(D, df = df)
        }
        if (D < 0) {
            v = 1
	    lambda=NULL
	    sigma=NULL
	    mu=NULL
	    R=NULL
	    alpha=NULL}
        else v = 2
    }
    a = list(chi.sq = D, df = df, p.value = alpha)
    a
}
