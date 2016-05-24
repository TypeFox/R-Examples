## the idea and most commands were provided by Marco J. Maier, Institute for
## Statistics and Mathematics, Vienna University of Economics and Business

library(maxLik)
options(warn = -1, digits = 4 )
set.seed(5)
some_data <- rgamma(1e4, shape = 5, scale = 2)

# log-likelihood function(s)
logLL <- function(x, X)   # per observation for maxLik
   dgamma(x = X, shape = exp(x[1]), scale = exp(x[2]), log = TRUE)
logLL_sum <- function(x, X)   # negative sum for nlm()
   -sum(dgamma(x = X, shape = exp(x[1]), scale = exp(x[2]), log = TRUE))

sum(logLL(log(c(5,2)),some_data))
logLL_sum(log(c(5,2)),some_data)
all.equal( sum(logLL(log(c(5,2)),some_data)), -logLL_sum(log(c(5,2)),some_data))

# gradient of log-likelihood function
d_logLL <- function(x, X){   # analytic 1. derivatives
   cbind(shape=exp(x[1])*(-x[2]-psigamma(exp(x[1]),0)+log(X)),
         scale= X / exp(x[2]) - exp(x[1]))
}

d_logLLNum <- function(x, X){
   numericGradient( logLL, x, X = X )
}

colSums(d_logLL(log(c(5,2)),some_data))
colSums(d_logLLNum(log(c(5,2)),some_data))

all.equal( d_logLL(log(c(5,2)),some_data), d_logLLNum(log(c(5,2)),some_data),
   check.attributes=FALSE)

# Hessian of log-likelihood function
dd_logLL <- function(x, X){   # analytic 2. derivatives
   grad <- d_logLL( x, X )
   hessian <- matrix(0, 2, 2)
   hessian[1,1] <- sum( grad[,1] - exp(x[1])^2 * psigamma(exp(x[1]), 1) )
   hessian[2,2] <- - sum( X / exp(x[2]) )
   hessian[cbind(c(2,1), c(1,2))] <- -exp(x[1]) * length(X)
   return(hessian)
}

dd_logLLNum <- function(x, X){
   numericHessian( function(x,X) sum(logLL(x,X)), t0=x, X = X )
}
dd_logLLNumGrad <- function(x, X){
   numericHessian( function(x,X) sum(logLL(x,X)), 
      grad = function(x,X) colSums(d_logLL(x,X)), x, X = X )
}

dd_logLL(log(c(5,2)),some_data)
dd_logLLNum(log(c(5,2)),some_data)
all.equal(dd_logLL(log(c(5,2)),some_data), dd_logLLNum(log(c(5,2)),some_data))
dd_logLLNumGrad(log(c(5,2)),some_data)
all.equal(dd_logLL(log(c(5,2)),some_data), dd_logLLNumGrad(log(c(5,2)),some_data),
   check.attributes=FALSE)

# estimation with nlm()
t_nlm <- system.time( r_nlm  <- nlm(logLL_sum, c(0,0), X=some_data, hessian=TRUE) )

# estimation with nlm() and gradients
logLL_grad <- function(x, X) {
   result <- logLL_sum( x, X )
   attr( result, "gradient" ) <- - colSums( d_logLL( x, X ) )
   return( result )
}
t_nlmg <- system.time( r_nlmg  <- nlm(logLL_grad, c(0,0), X=some_data, hessian=TRUE) )

# estimation with nlm() and gradients and Hessian
logLL_hess <- function(x, X) {
   result <- logLL_sum( x, X )
   attr( result, "gradient" ) <- - colSums( d_logLL( x, X ) )
   attr( result, "hessian" ) <- - dd_logLL( x, X )
   return( result )
}
t_nlmgh <- system.time( r_nlmgh  <- nlm(logLL_hess, c(0,0), X=some_data, hessian=TRUE) )

# estimation with optim() / BFGS
t_bfgs <- system.time( r_bfgs <- optim(c(0,0), logLL_sum, X=some_data, 
   method="BFGS", hessian=TRUE) )

# estimation with maxLik() / BFGS
t_bfgsM <- system.time( r_bfgsM <- maxLik( logLL, start = c(0,0), 
   method="BFGS", X=some_data ) )

# estimation with maxLik() / BFGS with gradients
t_bfgsMg <- system.time( r_bfgsMg <- maxLik( logLL, d_logLL, start = c(0,0), 
   method="BFGS", X=some_data ) )

# estimation with maxLik() / BHHH
t_bhhh <- system.time( r_bhhh <- maxLik( logLL, start = c(0,0), 
   method="BHHH", X=some_data ) )

# estimation with maxLik() / BHHH with gradients
t_bhhhg <- system.time( r_bhhhg <- maxLik( logLL, d_logLL, start = c(0,0), 
   method="BHHH", X=some_data ) )

# estimation with maxLik() / NR
t_NRn <- system.time( r_NRn <- maxLik( logLL, start = c(0,0), 
   method="NR", X=some_data ) )

# estimation with maxLik() / NR with gradients
t_NRg <- system.time( r_NRg <- maxLik( logLL, d_logLL, start = c(0,0), 
   method="NR", X=some_data ) )

# estimation with maxLik() / NR with gradients and Hessian
t_NRgh <- system.time( r_NRgh <- maxLik( logLL, d_logLL, dd_logLL, start = c(0,0), 
   method="NR", X=some_data ) )

# log likelihood values
rbind(NLM=-r_nlm$minimum, 
      NLM_grad=-r_nlmg$minimum,
      NLM_gradHess=-r_nlmgh$minimum,
      BFGS=-r_bfgs$value,
      maxLikBfgs = logLik( r_bfgsM ),
      maxLikBfgs_grad = logLik( r_bfgsMg ),
      BHHH = logLik( r_bhhh ),
      BHHH_grad = logLik( r_bhhhg ),
      NR_numeric= logLik( r_NRn ),
      NR_grad= logLik( r_NRg ),
      NR_gradHess= logLik( r_NRgh ) )
      

# estimated coefficients
pp <- exp(rbind(NLM=r_nlm$estimate, 
                NLM_grad=r_nlmg$estimate,
                NLM_gradHess=r_nlmgh$estimate,
                BFGS=r_bfgs$par,
                maxLikBfgs = coef( r_bfgsM ),
                maxLikBfgs_grad = coef( r_bfgsMg ),
                BHHH = coef( r_bhhh ),
                BHHH_grad = coef( r_bhhhg ),
                NR_numeric= coef( r_NRn ),
                NR_grad= coef( r_NRg ),
                NR_gradHess= coef( r_NRgh ) ))
colnames(pp) <- c("shape_alpha", "scale_theta")
pp


# some Hessians
-100*round(r_nlm$hessian/100,0)
round(solve(r_nlm$hessian),5)

-100*round(r_nlmg$hessian/100,0)
round(solve(r_nlmg$hessian),5)

-100*round(r_nlmgh$hessian/100,0)
round(solve(r_nlmgh$hessian),5)

-100*round(r_bfgs$hessian/100,0)
round(solve(r_bfgs$hessian),5)

100*round(r_NRn$hessian/100,0)
round(solve(-r_NRn$hessian),5)

100*round(r_NRg$hessian/100,0)
round(solve(-r_NRg$hessian),5)


# standard errors
se <- exp(rbind(NLM=sqrt(diag( solve(r_nlm$hessian) )), 
                NLM_grad=sqrt(diag( solve(r_nlmg$hessian) )),
                NLM_gradHess=sqrt(diag( solve(r_nlmgh$hessian) )),
                BFGS=sqrt(diag( solve(r_bfgs$hessian) )),
                maxLikBfgs = stdEr( r_bfgsM ),
                maxLikBfgs_grad = stdEr( r_bfgsMg ),
                BHHH = stdEr( r_bhhh ),
                BHHH_grad = stdEr( r_bhhhg ),
                NR_numeric= stdEr( r_NRn ),
                NR_grad= stdEr( r_NRg ),
                NR_gradHess= stdEr( r_NRgh ) ))
colnames(se) <- c("shape_alpha", "scale_theta")
se

# execution times
tt <- rbind(t_nlm, t_nlmg, t_nlmgh, t_bfgs, t_bfgsM, t_bfgsMg, 
            t_bhhh, t_bhhhg, t_NRn, t_NRg, t_NRgh )
# tt
