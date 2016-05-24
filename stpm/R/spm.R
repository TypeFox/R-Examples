#'A central function that estimates Stochastic Process Model parameters a from given dataset.
#'@references Yashin, A. et al (2007), Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.
#'@references Akushevich I., Kulminski A. and Manton K. (2005). Life tables with covariates: Dynamic model 
#'for Nonlinear Analysis of Longitudinal Data. Mathematical Popu-lation Studies, 12(2), pp.: 51-80.
#'<DOI: 10.1080/08898480590932296>.
#'@references Yashin, A. et al (2007), Health decline, aging and mortality: how are they related? 
#'Biogerontology, 8(3), 291-302.<DOI:10.1007/s10522-006-9073-3>.
#'@param x A dataset: is the output from prepare_data(...) function and consists of two separate data tables:
#'(1) a data table for continuous-time model and (2) a data table for discrete-time model.
#'@param model A model type. Choices are: "discrete", "continuous" or "time-dependent".
#'@param formulas A list of parameter formulas used in the "time-dependent" model.
#'@param start A starting values of coefficients in the "time-dependent" model.
#'@param tol A tolerance threshold for matrix inversion (NULL by default).
#'@param stopifbound A flag (default=FALSE) if it is set then the optimization stops 
#'when any of the parametrs achives lower or upper boundary.
#'@param algorithm An optimization algorithm used in \code{nloptr} package.
#'Default: \code{NLOPTR_NL_NELDERMEAD}.
#'@param lb Lower boundary, default \code{NULL}.
#'@param ub Upper boundary, default \code{NULL}.
#'@param maxeval Maximum number of evaluations of optimization algorithm. 
#'Default 100.
#'@param pinv.tol A tolerance threshold for matrix pseudo-inverse. Default: 0.01.
#'@param theta.range A user-defined range of the parameter \code{theta} used in 
#'discrete-time optimization and estimating of starting point for continuous-time optimization.
#'@param verbose A verbosing output indicator (FALSE by default).
#'@return For "discrete" and "continuous" model types: 
#'(1) a list of model parameter estimates for the discrete model type described in 
#'"Life tables with covariates: Dynamic Model for Nonlinear Analysis of Longitudinal Data", 
#'Akushevich et al, 2005.<DOI:10.1080/08898480590932296>,  and  
#'(2) a list of model parameter estimates for the continuous model type described in 
#'"Stochastic model for analysis of longitudinal data on aging and mortality", 
#'Yashin et al, 2007, Math Biosci.<DOI:10.1016/j.mbs.2006.11.006>.
#'
#'For the "time-dependent" model (model parameters depend on time): a set of model parameter estimates.
#'@examples \dontrun{ 
#'library(stpm)
#'#Prepare data for optimization
#'data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"), 
#'					   y=system.file("data","vitstat.csv",package="stpm"))
#'#Parameters estimation (default model: discrete-time):
#'p.discr.model <- spm(data)
#'p.discr.model
#'# Continuous-time model:
#'p.cont.model <- spm(data, model="continuous")
#'p.cont.model
#'#Model with time-dependent coefficients:
#'data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"), 
#'					   y=system.file("data","vitstat.csv",package="stpm"), 
#'					   covariates="BMI")
#'p.td.model <- spm(data, model="time-dependent",
#'                  f=list(at="aa*t+bb", f1t="f1", Qt="Q", ft="f", bt="b", mu0t="mu0"),
#'                  start=list(a=-0.001, bb=0.05, f1=80, Q=2e-8, f=80, b=5, mu0=1e-3))
#'p.td.model
#'}
spm <- function(x, model="discrete", formulas = NULL, start=NULL, tol=NULL, 
                stopifbound=FALSE, algorithm="NLOPT_LN_NELDERMEAD", 
                lb=NULL, ub=NULL, maxeval=100,
                pinv.tol = 0.01,
                theta.range=seq(0.01, 0.2, by=0.001),
                verbose=FALSE) {
  
  # List of available models:
  models <- c("discrete", "continuous", "time-dependent")
  
  if(!(model %in% models)) {
    stop(cat(model, " - unknown model type!"))
  }
  
  # Number of variables (dimensions):
  k <- (dim(x[[1]])[2] - 4)/2
  
  
  if(model == "discrete") {
    # Estimation of starting point with discrete optimization:
    pars <- spm_discrete(dat=x[[2]],verbose = verbose, tol = tol, theta_range=theta.range)
    res <- list(Ak2005=list(u=pars$Ak2005$u, 
                            R=pars$Ak2005$R, 
                            b=pars$Ak2005$b, 
                            Q=pars$Ak2005$Q, 
                            Sigma=pars$Ak2005$Sigma,
                            mu0=pars$Ak2005$mu0,
                            theta=pars$Ak2005$theta), 
                Ya2007=list(a=pars$Ya2007$a, 
                            f1=pars$Ya2007$f1,
                            Q=pars$Ya2007$Q,
                            f=pars$Ya2007$f, 
                            b=pars$Ya2007$b, 
                            mu0=pars$Ya2007$mu0, 
                            theta=pars$Ya2007$theta))
    
  }
  
  
  if(model == "continuous") {
    pars <- spm_discrete(dat=x[[2]],verbose = verbose, tol = tol)
    data <- data.frame(x[[1]])
  
    if(verbose) {
      cat("Starting parameters:\n")
      print(pars)
    }
    
    if(det(pars$Ya2007$Q) < 0) {
      cat("Error: determinant of Q < 0\n")
      cat("Q:\n")
      print(pars$Ya2007$Q)
      cat("Det(Q):\n")
      print(det(pars$Ya2007$Q))
      
      res <- NA
    
    } else {
      res.t <- spm_continuous(as.matrix(data), 
                    a=pars$Ya2007$a, 
                    f1=pars$Ya2007$f1, 
                    Q=pars$Ya2007$Q, 
                    f=pars$Ya2007$f, 
                    b=pars$Ya2007$b, 
                    mu0=pars$Ya2007$mu0, 
                    theta=pars$Ya2007$theta, 
                    stopifbound = stopifbound,
                    algorithm = algorithm, 
                    lb = lb, ub = ub,
                    maxeval = maxeval, 
                    pinv.tol = pinv.tol,
                    verbose = verbose)
  
      #res.t <- get("results",envir=.GlobalEnv)
      
      Q.c <- res.t$Q
      R.c <- res.t$a + diag(k)
      Sigma.c <- as.matrix(res.t$b)
      u.c <- (-1)*(t(res.t$f1) %*% res.t$a)
      b.c <- -2*t(res.t$f) %*% res.t$Q
      mu0.c <- res.t$mu0 + t(res.t$f) %*% res.t$Q %*% res.t$f
      theta.c <- res.t$theta
      
      res <- list(Ak2005=list(u=u.c, 
                              R=R.c, 
                              b=b.c, 
                              Q=Q.c, 
                              Sigma=Sigma.c,
                              mu0=mu0.c,
                              theta=theta.c), 
                  Ya2007=list(a=res.t$a, 
                              f1=res.t$f1,
                              Q=res.t$Q,
                              f=res.t$f, 
                              b=res.t$b, 
                              mu0=res.t$mu0, 
                              theta=res.t$theta))
      
    }
  }
  
  if(model == "time-dependent") {
    
    if(k > 1) {
      stop("Number of variables > 1. Model with time-dependent parameters can be used only with one variable!")
    }
    
    #if(length(formulas) != 6) {
    #  stop("It must be 6 equations for corresponding coefficients.")
    #}
    
    # Raw parameters estimates
    #pars <- spm_discrete(dat=x[[2]],verbose = verbose, tol = tol, theta_range=theta.range)
    # Parameter optimization for time-dependent model
    if(is.null(start)) {
      stop("Specify starting values.")
    }
    
    res <- spm_time_dep(x[[1]], 
                        f=formulas,
                        start=start,
                        algorithm=algorithm,
                        lb=lb, ub=ub,
                        verbose=verbose, maxeval=maxeval)
    
    #res <- get("results",envir=.GlobalEnv)
  }
  class(res) <- "spm"
  invisible(res)
}