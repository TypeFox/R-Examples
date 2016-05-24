#' Compute the so called "Student bridge" process.
#' @param X MCMC sampling sequence of length N
#' @return Student bridge sequence: \deqn{S=\left\{s_n\right\}_{1\leq n\leq N-1},  s_n=\sqrt{N-2} {\frac{n*(\hat \mu_{1,n}-\hat \mu_{n+1,N})}{\sqrt{\left ({1 \over n} + {1 \over {N-n}} \right) * ((n-1) \hat{\sigma_{1,n}}^2+(N-n-1)\hat{\sigma_{n+1,N}}^2} }}}
#' @examples
#' x = AR1(rho=0)
#' sb = studentbridge(x)
#' plot(sb,type='l',col='blue')
studentbridge <- function(X) {
    N = length(X)
    s = array(NaN, N-3)
    for (n in 2:(N-2)) {
        s[n-1] = sqrt(N-2) * (mean(X[1:n]) - mean(X[(n+1):N])) / sqrt((1/n + 1/(N-n)) * ( (n-1)*var(X[1:n]) + (N-n-1)*var(X[(n+1):N]) ))
    }
    return(s)
}


#' Compute the so called (abusively) "Brownian bridge" process.
#' @param X MCMC sampling sequence of length N
#' @return cumulative normalized sum sequence: \deqn{B=\left\{b_n\right\}_{0\leq n\leq N},  b_n=\frac{n*(\hat \mu_{1,n}-\hat \mu_{1,N})}{\hat{\sigma} \sqrt(N)}}
#' @examples
#' x = AR1(rho=0)
#' bb = brownianbridge(x)
#' plot(bb,type='l',col='red')
brownianbridge <- function(X) {
    N = length(X)
    b = array(NaN, N+1)
    b[1] = 0
    for (n in 1:N) {
        b[n+1] = n * (mean(X[1:n]) - mean(X[1:N])) / (sqrt(N) * sd(X))
    }
    return(b)
}

#' Compute the so called "Log-likelihood bridge" process.
#' @param X MCMC sampling sequence of length N
#' @return log-likelihood sequence \deqn{LL=\left\{ll_n\right\}_{2\leq n\leq N-2},  ll_n=N*ln(\hat \sigma_{1,N}^2)-n*ln(\hat \sigma_{1,n}^2)-(N-n)ln(\hat \sigma_{n+1,N}^2)}
#' @examples
#' x = AR1(rho=0)
#' llb = loglikbridge(x)
#' plot(llb,type='l',col='red')
loglikbridge <- function(X) {
    N = length(X)
    ll = array(NaN, N-3)
    for (n in 2:(N-2)) {
        ll[n-1] = N*log(var(X[1:N])) - n*log(var(X[1:n])) - (N-n)*log(var(X[(n+1):N]))
    }
    return(ll)
}

