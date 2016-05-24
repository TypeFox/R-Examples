#' Power calculations for chi-squared test
#'
#' @param es              effect size. A numeric value or output of ES.chisq.gof, ES.chisq.assoc
#' @param df              degree of freedom
#' @param n               total number of observations
#' @param power           power of study
#' @param sig.level       significance level
#' @seealso               \code{\link{ES.chisq.gof}}
#' @seealso               \code{\link{ES.chisq.assoc}}
#' @seealso               \code{\link{power.plot.chisq}}
#' @export
#' @examples
#' ## calculate power
#' power.chisq(es=0.16,df=1,n=530,sig.level=0.05)
#' 
#' ## calculate sig.level
#' power.chisq(es=0.16,df=1,n=530,power=0.9576)
#'
#' ## calculate sample size
#' power.chisq(es=0.16,df=1,power=0.9576,sig.level=0.05)
#' 
#' ## calculate effect size
#' power.chisq(df=1,n=530,power=0.9576,sig.level=0.05)
power.chisq <- function(es=NULL,df=NULL,n=NULL,power=NULL,sig.level=NULL){
  pwr=NULL
  if(is.list(es) & length(es)==5){
    myes=es$es
    df=es$df
    n=es$n
    if (sum(sapply(list(power, sig.level), is.null)) != 1){
      stop("exactly one of power and sig.level must be NULL\n")
    }else if(is.null(power)){
      temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
      power=pchisq(temp, df = df, ncp = n * myes^2, lower.tail = FALSE)
    }else{
      func1 <- function(sig.level){
        temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
        pchisq(temp, df = df, ncp = n * myes^2, lower.tail = FALSE) - power
      }
      sig.level=uniroot(func1,lower=1e-15,upper=1-1e-15)$root
    }
  }else{
    if(!is.numeric(df) || df<1){
      stop("df must be at least 1\n")
    }
    if(sum(sapply(list(es, n, power, sig.level), is.null)) != 1){
      stop("exactly one of es, n, power, and sig.level must be NULL\n")
    }
    if(!is.null(es) && (!is.numeric(es) || es<0)){
      stop("es must be positive\n")
    }
    if(!is.null(n) && (!is.numeric(n) || n<1)){
      stop("n must be at least 1\n")
    }
    if(!is.null(power) && (!is.numeric(power) || power>1 || power<0)){
      stop("power must be in [0,1]\n")
    }
    if(!is.null(sig.level) && (!is.numeric(sig.level) || sig.level>1 || sig.level<0)){
      stop("sig.level must be in [0,1]\n")
    }
    
    if(is.null(power)){
      temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
      power=pchisq(temp, df = df, ncp = n * es^2, lower.tail = FALSE)
    }else if(is.null(n)){
      func2 <- function(n){
        temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
        pchisq(temp, df = df, ncp = n * es^2, lower.tail = FALSE) - power
      }
      n=uniroot(func2,lower=1+1e-15,upper=1e+15)$root
    }else if(is.null(es)){
      func3 <- function(es){
        temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
        pchisq(temp, df = df, ncp = n * es^2, lower.tail = FALSE) - power
      }
      es=uniroot(func3,lower=1e-15,upper=1-1e-15)$root
    }else if(is.null(sig.level)){
      func4 <- function(sig.level){
        temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
        pchisq(temp, df = df, ncp = n * es^2, lower.tail = FALSE) - power
      }
      sig.level=uniroot(func4,lower=1e-15,upper=1-1e-15)$root
    }
  }
  n=ceiling(n)
METHOD <- "Chi squared power calculation"
NOTE <- "'n' is the number of observations"
structure(list(Power = power, Effect_Size=es, df=df, n=n,sig.level=sig.level, note=NOTE,method = METHOD), class = "power.htest")
}
