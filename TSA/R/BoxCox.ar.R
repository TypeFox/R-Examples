`BoxCox.ar` <-
function (y, order, lambda=seq(-2,2,.01), plotit = TRUE,
method=c("mle","yule-walker", "burg", "ols", "yw"), ...) 
{
# Find the power transformation assuming the model 
# y^lambda=(y**lambda-1)/lambda, for lambda not =0
# y^lambda=log(y) for lambda=0
# The model is y^lambda(t)=phi_0+phi_1*y^lambda(t-1)+...+phi_k*y^lambda(t-p)+e(t)
# That is we assume an AR(k) model for y^lambda, with known or unknown order.
# All models are fitted by Yule-Walker method with the noise variance estimated
# by the formula sigma**2=data variance/(1-phi_1*rho_1-...-phi_k*rho_k)
# We approximate the likelihood method by approximating the loglikelihood by
# -n/2*log(sigma**2)+Jacobian due to the power transformation.
# 
# Input: y=time series (must be positive)
#        order= AR order 
#        if order is missing, it will be estimated by minimizing AIC for each
#        power lambda.
#        lambda=the vector of powers over which the transformation is searched.
#        plotit=TRUE (default) whether or not to plot the "log-likelihood" for
#               choosing the power.
#
# Output: a list containing the following components
#         lambda
#         loglike=the vector of loglikelihood of the power parameter
#         mle="maximum likelihood estimate" of lambda
#         ci= approximate 95% C.I. of lambda (note that the algorithm 
#             assumes that the set over which the likelihood is within 1.92
#             of the maximum log likelihood is an interval. If this assumption
#             is incorrect, then the interval is incorrect. So check the diagram.
#
#         It is useful to try different lambda vector to zoom in the log-likelihood
#         to get a good picture.
#
#  Example
#
#         > hare.transf=BoxCox.ar(y=hare2,lambda=seq(-0.1,0.7,.01))
#         > hare.transf$mle
#         [1] 0.36
#         > hare.transf$ci
#         [1] 0.14 0.65
#
#  So, the mle of lambda is .36 and the approximate 95% confidence interval
#  (0.14, 0.65) which includes 0.5. Hence, we may apply the square root
#  transformation to the hare data.
#
#
                    
                    if(missing(method)) method='mle'  
                    y=as.vector(y/(max(abs(y))+1)) 
                    if (any(y <= 0)) 
                stop("Data values must be positive")
                order=ar(log(y),method='mle')$order
                nlngmy <- sum(log(y))
        if (!missing(lambda)) 
                xl <- lambda
        else xl <- seq(-2, 2, 0.1)
        loglik <- as.vector(xl)
        for (i in 1:length(xl)) if (abs(xl[i]) > 0.00) {
                if( missing(order)) ar.result=ar((y^xl[i]-1)/xl[i],method=method) else 
                                   ar.result=ar((y^xl[i]-1)/xl[i],method=method,order.max=order)
                 n=length(y)-ar.result$order
                 ar.res=ar.result$resid
                 n=length(y)
                loglik[i] <-  -n/2 * log(ar.result$var.pred) + 
                (xl[i] - 1) * nlngmy}
        else {
                if( missing(order)) ar.result=ar(log(y),method=method) else 
                                   ar.result=ar(log(y),method=method,order.max=order)
                 n=length(y)-ar.result$order
                 ar.res=ar.result$resid
                 n=length(y)

                 loglik[i] <- -n/2 * log(ar.result$var.pred) - 
                nlngmy}
        if (plotit) {
                plot(xl, loglik, xlab = expression(lambda), ylab = "Log Likelihood", 
                        type = "l", ylim=c(min(loglik),max(loglik)))
                lambdahat <- loglik[loglik == max(loglik)]
                limit <- lambdahat - 0.5 * qchisq(0.95, 1)
                in.interval=xl[loglik>=limit]
                lower=in.interval[1]
                upper=rev(in.interval)[1]
                mle=(xl[loglik==max(loglik)])[1]
               lines(x=c(lower,lower),y=c(min(loglik),limit),lty=2)
               lines(x=c(upper,upper),y=c(min(loglik),limit),lty=2)
               lines(x=c(mle,mle),y=c(min(loglik),max(loglik)),lty=2)

                abline(limit, 0,lty=2)
                scal <- (par("usr")[4] - par("usr")[3])/par("pin")[2]
                text(c(xl[1])+0.1, limit + 0.08 * scal, " 95%")
        }
        invisible(list(lambda = xl, loglike = loglik, mle=mle, ci=c(lower,upper)))
}

