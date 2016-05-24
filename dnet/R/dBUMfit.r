#' Function to fit a p-value distribution under beta-uniform mixture model
#'
#' \code{dBUMfit} is supposed to take as input a vector of p-values for deriving their distribution under beta-uniform mixture model (see Note below). The density distribution of input p-values is expressed as a mixture of two components: one for the null hypothesis (the noise component) and the other for the alternative hypothesis (the signal component). The noise component is the uniform density, while the signal component is the remainder of the mixture distribution. It returns an object of class "BUM".
#'
#' @param x a vector containing input p-values
#' @param ntry an integeter specifying how many trys are used to find the optimised parameters by maximum likelihood estimation
#' @param hist.bum logical to indicate whether the histogram graph should be drawn
#' @param contour.bum logical to indicate whether a contour plot should be drawn to show the log likelihood as a function of two parameters (a and lambda) in the beta-uniform mixture model
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' an object of class "BUM", a list with following elements:
#' \itemize{
#'  \item{\code{lambda}: estimated mixture parameter}
#'  \item{\code{a}: estimated shape parameter}
#'  \item{\code{NLL}: Negative log-likelihood}
#'  \item{\code{pvalues}: the input pvalues}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The probability density function of p-values under the Beta-Uniform Mixture model is formulated as: \eqn{f(x|\lambda,a) = \lambda + (1-\lambda)*a*x^{a-1}}. The model names after mixing two distributions: 
#' \itemize{
#' \item{the uniform distribution with the density function as \eqn{\frac{1}{b-a}|_{a=0}^{b=1}=1}}
#' \item{the beta distribution with the density function as \eqn{\frac{\Gamma(a+b)}{\Gamma(a)+\Gamma(b)}*x^{a-1}*(1-x)^{b-1}|_{b=1}=a*x^{a-1}}}
#' }
#' Both are mixed via \eqn{\lambda}. The mixture parameter \eqn{\lambda} measures the contribution from the uniform distribution. Accordingly, \eqn{1-\lambda} measures the contribution from the beta distribution. Notably, the probability density function of the beta distribution can be splitted into two parts (rather than the exclusitive signal): 
#' \itemize{
#' \item{the constant part as noise: \eqn{a*x^{a-1}|_{x=1}=a}}
#' \item{the rest part as signal: \eqn{a*(x^{a-1}-1)}}
#' }
#' In other words, there is no signal at \eqn{x=1} but all being noise. It is a conservative, upper bound estimation of the noise. Therefore, the probability density function in the model can be decomposed into signal-noise components:
#' \itemize{
#' \item{the signal component: \eqn{(1-\lambda)*a*(x^{a-1}-1)}}
#' \item{the noise component: \eqn{\lambda + (1-\lambda)*a}}
#' }
#' It is misleading to simply view \eqn{\lambda} as the noise component and \eqn{(1-\lambda)*a*x^{a-1}} as the signal component, just as wrongly do in the literatures (e.g. \url{http://www.ncbi.nlm.nih.gov/pubmed/18586718})
#' @export
#' @seealso \code{\link{dBUMscore}}
#' @include dBUMfit.r
#' @examples
#' # 1) generate an vector consisting of random values from beta distribution
#' x <- rbeta(1000, shape1=0.5, shape2=1)
#'
#' # 2) fit a p-value distribution under beta-uniform mixture model
#' fit <- dBUMfit(x)
#' fit$lambda
#' fit$a

dBUMfit <- function(x, ntry=1, hist.bum=T, contour.bum=T, verbose=T)
{

    ## Initial values for the parameters (a and lambda) to be optimized over
    a <- stats::runif(ntry, 0.1, 0.9)
    lambda <- stats::runif(ntry, 0.1, 0.9)
  
    ## beta-uniform mixture model with shape2=1
    fbum <- function(x, lambda, a){
        lambda+(1-lambda)*a*x^(a-1)
    }
    ## A function to be minimized: negative log likelihood of BUM model
    fn <- function(parms, x){
        -1*sum(log(fbum(x, parms[1], parms[2])))
    }
    
    # A function to return the gradient for "L-BFGS-B" methods
    gr <- function(parms, x){
        lambda <- parms[1]
        a <- parms[2]
        
        deno <- lambda+(1-lambda)*a*x^(a-1)
        
        ## partial erivative in terms of lambda
        d_lambda <- -1*sum( (1-a*x^(a-1)) / deno )
        ## partial erivative in terms of a
        d_a <- -1*sum( ((1-lambda)*x^(a-1) + a*(1-lambda)*x^(a-1)*log(x)) / deno )

        return(c(d_lambda,d_a))
    }
    
    value <- Inf
    best <- list()
    for(i in 1:ntry){
        test.optim <- try(opt <- stats::optim(c(lambda[i],a[i]), fn=fn, gr=gr, x=x, lower=rep(1e-5,3), method="L-BFGS-B", upper=rep(1-1e-5,3)))
        
        if ((!class(test.optim)=="try-error") && all(opt$par >= 1e-5) && all(opt$par <= 1-1e-5)){
            
            if(opt$value < value){
                value <- opt$value
                best <- opt
            }
            
            if(0){
                message(sprintf("Try: %d\t with Log-likelihood: %.1f, Mixture parameter (lambda): %1.3f, Shape parameter (a): %1.3f", i,-1*opt$value,opt$par[1],opt$par[2]))
            }
        }
    }
    
    if(length(best)==0){
        return(warning("BUM model could not be fitted to data"))
    }else{
        if (any(opt$par == 1e-5) || any(opt$par == 1-1e-5)){
            #warning("One or both parameters are on the limit of the defined parameter space")
        }
        
        fit <- list( lambda = best$par[1], 
                     a= best$par[2],
                     NLL = best$value, 
                     pvalues = x,
                     call = match.call()
                     )
        class(fit) <- "BUM"
        
        if(verbose){
            message(sprintf("\tA total of p-values: %d", length(fit$pvalues)), appendLF=T)
            message(sprintf("\tMaximum Log-Likelihood: %.1f", -1*fit$NLL), appendLF=T)
            message(sprintf("\tMixture parameter (lambda): %1.3f", fit$lambda), appendLF=T)
            message(sprintf("\tShape parameter (a): %1.3f", fit$a), appendLF=T)
        }
        
        if(hist.bum){
            
            grDevices::dev.new()
            
            ## A function to return the upper bound of pi (when pvalue equels 1)
            piUbound <- function (x){
                return(x$lambda + (1 - x$lambda) * x$a)
            }
    
            graphics::hist(fit$pvalues, breaks=100, probability=T, main="Histogram of p-values and Beta-Uniform-Mixture (BUM) model", xlab="P-values", ylab="Density (%)")
            px <- seq(from=0, to=1, 1/100)
            graphics::lines(px, fit$lambda+(1-fit$lambda)*fit$a*px^(fit$a-1), lwd=3, col="darkblue")
            graphics::lines(px, piUbound(fit)*(px<=1), col="darkgreen", lwd=2)
            graphics::lines(px, fit$lambda+(1-fit$lambda)*fit$a*px^(fit$a-1)-piUbound(fit), lwd=2, col="darkred")
            
            #abline(h=piUbound(fit), col="darkgreen", lwd=2)
            #axis(side=2, labels=expression(pi), at=piUbound(fit))
            
            leg.txt <- c(
            expression(paste("Density fitted under BUM model: ", f(x), "=", lambda+(1-lambda)*a*x^(a-1), sep="")),
            expression(paste("With density for the noise component: ", pi, "=", lambda+(1-lambda)*a, sep="")),
            expression(paste("With density for the signal component: ", f(x)-pi, "=", (1-lambda)*a*(x^(a-1)-1), sep=""))
            )
            graphics::legend("top", legend=leg.txt, pch=15, col=c("darkblue", "darkgreen", "darkred"), border="transparent", box.col="transparent", cex=0.8)
            
        }
        
        if(contour.bum){
        
            grDevices::dev.new()
        
            v <- seq(0.05, 0.95, 0.05)
            z <- matrix(0, nrow=length(v), ncol=length(v))
            for(i in 1:length(v)){
                for(j in 1:length(v)){
                    z[i,j] <- -1*fn(c(v[i],v[j]), fit$pvalues)
                }
            }
  
            Lines <- list(bquote(lambda==.(round(fit$lambda,3))), bquote(a==.(round(fit$a,3))))
  
            graphics::filled.contour(v, v, z, 
                nlevels=24, color.palette=grDevices::colorRampPalette(unlist(strsplit("darkblue-lightblue-white-lightyellow-darkorange","-"))),
                main="Log-likelihood as a function of parameters", xlab=expression(lambda), ylab="a",
                plot.axes={
                    graphics::axis(1, seq(0,1,0.1))
                    graphics::axis(2, seq(0,1,0.1))       
                    graphics::abline(v=fit$lambda, lty=2, col="black")
                    graphics::abline(h=fit$a, lty=2, col="black")
                    graphics::points(fit$lambda, fit$a, cex=2)
                    graphics::text(fit$lambda, fit$a+(graphics::strheight("X")*1.5*seq(length(Lines))), do.call(expression, Lines), adj=c(-0.2,0))
                }
            )
        }
        
        return(fit)
    }
}