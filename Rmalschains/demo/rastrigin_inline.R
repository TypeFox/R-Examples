# This demo is strongly inspired by an 
# example from the RcppDE package

# It shows how to use the interface to pass the objective function 
# directly as an external pointer, which may speed things up drastically.
# Furthermore it shows how to use an environment to pass 
# additional parameters (e.g., data) to the objective function. This even works
# if the objective function is implemented in C++

library(inline)

inc <- 'double rastrigin(SEXP xs, SEXP env) {
    Rcpp::NumericVector x(xs);
    Rcpp::Environment e(env);
    
    double sum = e["target.value"];

    int n = x.size();
    
    for (int i=0; i<n; i++) {
    sum += x[i]*x[i] - 10*cos(2*M_PI*x[i]) + 10;
    
    }
    return(sum);
    }
    '

# define the function that gets the external pointer to the target function
src.xptr <- '
    typedef double (*funcPtr)(SEXP, SEXP);
    return(XPtr<funcPtr>(new funcPtr(&rastrigin)));
    '
create_xptr <- cxxfunction(signature(), body=src.xptr, inc=inc, plugin="Rcpp")


rastrigin <- function(x) {
  
  dimension <- length(x)
  
  res <- target.value
  for (i in 1:dimension) {
    res <- res + (x[i]*x[i] - 10.0*cos(2.0*pi*x[i]) + 10.0)
  }
  
  res 
}

env <- environment(fun=rastrigin)
env[["target.value"]] <- 20 

n <- 10

library(Rmalschains)

set.seed(5)
time.inline <- system.time(fit.inline <- malschains(fn=create_xptr(), env=env, lower=rep(-25, n), upper=rep(25, n), maxEvals=50000,
    control=malschains.control(popsize=50, istep=300, optimum=20, ls="cmaes")))

set.seed(5)
time.normal <- system.time(fit.normal <- malschains(rastrigin, env=env, lower=rep(-25, n), upper=rep(25, n), maxEvals=50000, 
    control=malschains.control(popsize=50, istep=300, optimum=20, ls="cmaes")))

fit.normal
fit.inline

time.normal
time.inline
