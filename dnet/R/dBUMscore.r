#' Function to transform p-values into scores according to the fitted beta-uniform mixture model and/or after controlling false discovery rate
#'
#' \code{dBUMscore} is supposed to take as input a vector of p-values, which are transformed into scores according to the fitted beta-uniform mixture model. Also if the FDR threshold is given, it is used to make sure that p-values below this are considered significant and thus scored positively. Instead, those p-values above the given FDR are considered insigificant and thus scored negatively. 
#'
#' @param fit an object of class "BUM"
#' @param method the method used for the transformation. It can be either "pdf" for the method based on the probability density function of the fitted model, or "cdf" for the method based on the cumulative distribution function of the fitted model
#' @param fdr the given FDR threshold. By default, it is set to NULL, meaning there is no constraint. If given, those p-values with the FDR below this are considered significant and thus scored positively. Instead, those p-values with the FDR above this given FDR are considered insigificant and thus scored negatively
#' @param scatter.bum logical to indicate whether the scatter graph of scores against p-values should be drawn. Also indicated is the p-value (called tau) corresponding to the given FDR threshold (if any)
#' @return
#' \itemize{
#'  \item{\code{scores}: a vector of scores}
#' }
#' @note The transformation from the input p-value \eqn{x} to the score \eqn{S(x)} is based on the fitted beta-uniform mixture model with two parameters \eqn{\lambda} and \eqn{a}: \eqn{f(x|\lambda,a) = \lambda + (1-\lambda)*a*x^{a-1}}. Specifically, it considers the log-likelyhood ratio between the signal and noise compoment of the model. The probability density function (pdf) of the signal component and the noise component are \eqn{(1-\lambda)*a*(x^{a-1}-1)} and \eqn{\lambda + (1-\lambda)*a}, respectively. Accordingly, the cumulative distribution function (cdf) of the signal component and the noise component are \eqn{\int_0^x (1-\lambda)*a*(x^{a-1}-1) \, \mathrm{d}x} and \eqn{\int_0^x \lambda+(1-\lambda)*a \, \mathrm{d}x}. In order to take into account the significance of the p-value, the \eqn{fdr} threshold is also used for down-weighting the score. According to how to measure both components, there are two methods implemented for deriving the score \eqn{S(x)}:
#' \itemize{
#' \item{The method "pdf": \eqn{S(x) = log_2\frac{(1-\lambda)*a*(x^{a-1}-1)}{\lambda+(1-\lambda)*a} - log_2\frac{(1-\lambda)*a*(\tau^{a-1}-1)}{\lambda+(1-\lambda)*a} = log_2\big(\frac{x^{a-1}-1}{\tau^{a-1}-1}\big)}. For the purpose of down-weighting scores, it must ensure \eqn{log_2\frac{(1-\lambda)*a*(\tau^{a-1}-1)}{\lambda+(1-\lambda)*a} \geq 0 }, that is, the constraint via \eqn{\tau \leq \big(\frac{\lambda+2*a*(1-\lambda)}{a*(1-\lambda)}\big)^{\frac{1}{a-1}}}}
#' \item{The method "cdf": \eqn{S(x) = log_2\frac{\int_0^x (1-\lambda)*a*(x^{a-1}-1) \, \mathrm{d}x}{\int_0^x \lambda+(1-\lambda)*a \, \mathrm{d}x} - log_2\frac{\int_0^\tau (1-\lambda)*a*(\tau^{a-1}-1) \, \mathrm{d}x}{\int_0^\tau \lambda+(1-\lambda)*a \, \mathrm{d}x} = log_2\frac{(1-\lambda)*(x^{a-1}-a)}{\lambda+(1-\lambda)*a} - log_2\frac{(1-\lambda)*(\tau^{a-1}-a)}{\lambda+(1-\lambda)*a} = log_2\big(\frac{x^{a-1}-a}{\tau^{a-1}-a}\big)}. For the purpose of down-weighting scores, it must ensure \eqn{log_2\frac{(1-\lambda)*(\tau^{a-1}-a)}{\lambda+(1-\lambda)*a} \geq 0 }, that is, the constraint via \eqn{\tau \leq \big(\frac{\lambda+2*a*(1-\lambda)}{1-\lambda}\big)^{\frac{1}{a-1}}}}
#' \item{Where \eqn{\tau =\big[\frac{\lambda+(1-\lambda)*a-fdr*\lambda}{fdr*(1-\lambda)}\big]^{\frac{1}{a-1}}}, i.e. the p-value corresponding to the exact \eqn{fdr} threshold. It can be deduced from the definition of the false discovery rate: \eqn{fdr \doteq \frac{\int_0^\tau \lambda+(1-\lambda)*a \, \mathrm{d}x}{\int_0^\tau \lambda+(1-\lambda)*a*x^{a-1} \, \mathrm{d}x}}. Notably, if the calculated \eqn{\tau} exceeds the contraint, it will be reset to the maximum end of that constraint}
#' }
#' @export
#' @seealso \code{\link{dBUMfit}}
#' @include dBUMscore.r
#' @examples
#' # 1) generate an vector consisting of random values from beta distribution
#' x <- rbeta(1000, shape1=0.5, shape2=1)
#'
#' # 2) fit a p-value distribution under beta-uniform mixture model
#' fit <- dBUMfit(x)
#'
#' # 3) calculate the scores according to the fitted BUM and fdr=0.01
#' # using "pdf" method
#' scores <- dBUMscore(fit, method="pdf", fdr=0.01)
#' # using "cdf" method
#' scores <- dBUMscore(fit, method="cdf", fdr=0.01)

dBUMscore <- function(fit, method=c("pdf","cdf"), fdr=NULL, scatter.bum=T)
{

    if (class(fit) != "BUM"){
        stop("The funciton must apply to 'BUM' object.\n")
    }

    method <- match.arg(method)

    a <- fit$a
    lambda <- fit$lambda
    pval <- fit$pvalues
    
    ## A function to calculate fdr for each p-value
    f_fdr <- function(x){
        fp <- (lambda+(1-lambda)*a)*x
        tp <- (1-lambda)*(x^a-a*x)
        fp/(fp+tp)
    }
    
    ## the given fdr should be no more than the maximum fdr
    fdr_max <- f_fdr(1)
    if(is.null(fdr)){
        fdr <- fdr_max
    }else{
        if(fdr > fdr_max){
            fdr <- fdr_max
        }else if(fdr==0){
            fdr <- fdr_max
        }
    }
    
    ## calculate tau for a given FDR cutoff
    tmp_nomi <- lambda+(1-lambda)*a-fdr*lambda
    tmp_deno <- fdr*(1-lambda)
    tau <- (tmp_nomi/tmp_deno)^(1/(a-1))
    
    if (method == "cdf"){
        ## Score for the noise component: defined as an integral of the noise density over the internal [0,x]
        f_score_noise <- function(x){
            (lambda+(1-lambda)*a)*x
        }
        ## Score for the signal component: defined as an integral of the signal density over the internal [0,x]
        f_score_signal <- function(x){
            (1-lambda)*(x^a-a*x)
        }

        ## Score is defined as the log-likehood ratio between sigal score and noise score
        f_score <- function(x){
            log2( f_score_signal(x) / f_score_noise(x) )
        }
        
        ## calculate tau_bound when the noise = the signal
        tmp_nomi <- lambda+2*(1-lambda)*a
        tmp_deno <- (1-lambda)
        tau_bound <- (tmp_nomi/tmp_deno)^(1/(a-1))
        if(tau >= tau_bound){
            tau <- tau_bound
            fdr <- f_fdr(tau)
        }
        
    }else if (method == "pdf"){
    
        ## Score for the noise component: defined as the noise density
        f_score_noise <- function(x){
            lambda+(1-lambda)*a
        }
        ## Score for the signal component: defined as the signal density
        f_score_signal <- function(x){
            (1-lambda)*a*(x^(a-1)-1)
        }

        ## Score is defined as the log-likehood ratio between sigal score and noise score
        f_score <- function(x){
            log2(f_score_signal(x) / f_score_noise(x))
        }
        
        ## calculate tau_bound when the noise = the signal
        tmp_nomi <- lambda+2*(1-lambda)*a
        tmp_deno <- a*(1-lambda)
        tau_bound <- (tmp_nomi/tmp_deno)^(1/(a-1))
        if(tau >= tau_bound){
            tau <- tau_bound
            fdr <- f_fdr(tau)
        }
    }
    
    ## calculate score after taking into account the given FDR
    ## p-values below this are considered significant and thus scored positively
    ## Instead, those p-values above the given FDR are considered insigificant and thus scored negatively
    ## also make sure that f_score(tau) is always positive
    scores <- f_score(pval) - f_score(tau)
    
    if(scatter.bum){
        
        grDevices::dev.new()
        
        plot(pval, scores, pch=".", main="Scores vs P-values", xlab="P-values", ylab="Scores")
        
        graphics::abline(v=tau, lty=2, col="darkgray")
        graphics::abline(h=0, lty=2, col="darkgray")
        graphics::points(tau, 0, cex=1.5)
        if(tau == tau_bound | tau >= 1){
            Lines <- list(bquote(tau==.(round(tau,3))), bquote(fdr>=.(round(fdr,3))))
        }else{
            Lines <- list(bquote(tau==.(round(tau,3))), bquote(fdr==.(round(fdr,3))))
        }
        if(tau < 0.5){
            graphics::text(tau, 0+(graphics::strheight("X")*1.5*seq(length(Lines))), do.call(expression, Lines), adj=c(-0.2,0), cex=0.8)
        }else{
            graphics::text(tau, 0+(graphics::strheight("X")*1.5*seq(length(Lines))), do.call(expression, Lines), adj=c(1,0), cex=0.8)
        }
    }
    
    return(scores)
}
