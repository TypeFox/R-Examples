#'A data projection with previously estimated or user-defined parameters. 
#'Projections are constructed for a cohort with fixed or 
#'normally distributed initial covariates. 
#'@references Yashin, A. et al (2007), Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.
#'@references Akushevich I., Kulminski A. and Manton K. (2005). Life tables with covariates: Dynamic model 
#'for Nonlinear Analysis of Longitudinal Data. Mathematical Popu-lation Studies, 12(2), pp.: 51-80.
#'<DOI: 10.1080/08898480590932296>.
#'@references Yashin, A. et al (2007), Health decline, aging and mortality: how are they related? 
#'Biogerontology, 8(3), 291-302.<DOI:10.1007/s10522-006-9073-3>.
#'@param x A list of parameters from output of the \code{spm(...)} function.
#'@param N A number of individual to simulate, N=100 by default.
#'@param ystart A vector of starting values of covariates (variables), ystart=80 by default.
#'@param model A model type. Choices are: "discrete", "continuous" or "time-dependent".
#'@param f A list of formulas for the time-dependent model (NULL by default).
#'@param tstart Start time (age), default=30.
#'@param tend End time (age), default=105.
#'@param dt A time interval between observations, dt=1 by default.
#'@param sd0 A standard deviation value for simulation of the next value of variable.
#'sd0=4 by default.
#'@return An object of 'spm.projection' class with two elements. 
#'(1) A simulated data set.
#'(2) A summary statistics which includes (i) age-specific means of state variables and
#'(ii) Survival probabilities.
#'@examples \dontrun{ 
#'library(stpm)
#'# Setting up the model
#'model.par <- list()
#'model.par$a <- matrix(c(-0.05, 1e-3, 2e-3, -0.05), nrow=2, ncol=2, byrow=TRUE)
#'model.par$f1 <- matrix(c(90, 35), nrow=1, ncol=2)
#'model.par$Q <- matrix(c(1e-8, 1e-9, 1e-9, 1e-8), nrow=2, ncol=2, byrow=TRUE)
#'model.par$f <- matrix(c(80, 27), nrow=1, ncol=2)
#'model.par$b <- matrix(c(6, 2), nrow=2, ncol=2)
#'model.par$mu0 <- 1e-6
#'model.par$theta <- 0.09
#'# Projection
#'# Discrete-time model
#'data.proj.discrete <- spm_projection(model.par, N=5000, ystart=c(80, 27))
#'plot(data.proj.discrete$stat$srv.prob)
#'# Continuous-time model
#'data.proj.continuous <- spm_projection(model.par, N=5000, ystart=c(80, 27), model="continuous")
#'plot(data.proj.continuous$stat$srv.prob)
#'}
spm_projection <- function(x, 
                           N=100, 
                           ystart=80, 
                           model="discrete", 
                           f = NULL,
                           tstart=30, tend=105, 
                           dt=1, 
                           sd0=4) {
  
  avail.models <- c("discrete", "continuous", "time-dependent")
  if(!(model %in% avail.models)) {
    stop(paste("Provided model", model, "not found in the list of available models."))
  }
  
  res <- list()
  
  if(model == "time-dependent") {
    # Data simulation for time-dependent model
    
    formulas.work <- list(at="-0.05", f1t="80", Qt="2e-8", ft="80", bt="5", mu0t="2e-5")
    
    if (!is.null(f)) {
      for(item in f) {
        formulas.work[[item]] <- f[[item]]
      }
    }
    
    #Simulate (project) data:
    res.time_dep <- simdata_time_dep(N=N,f=f,
                                    step=dt, 
                                    tstart=tstart, 
                                    tend=tend, 
                                    ystart=ystart, 
                                    sd0=sd0)
    #Compute summary statistics:
    # Statistics
    stat <- list()
    # Age-specific means:
    bins<-10
    cutpoints<-quantile(res.time_dep[,3],(0:bins)/bins)
    binned <-cut(res.time_dep[,3],cutpoints,include.lowest=TRUE)
    for(i in seq(5,length(colnames(res.time_dep)),by=2)) {
      mean.cov <- tapply(res.time_dep[,5], binned, mean)
      stat[["mean.by.age"]][[colnames(res.time_dep)[i]]] <- mean.cov
    }
    
    #Survival probabilities:
    srv.prob <- survfit( Surv(t1, xi) ~ 1, conf.type="none", data = data.frame(res.time_dep))
    stat[["srv.prob"]] <- srv.prob
    
    res <- list(data=res.time_dep, stat=stat)
    
  } else if(model == "discrete") {
    
    if(length(x$a) != length(ystart)^2) {
      stop("Number of dimensions does not match with the number of values provided in ystart.")
    }
    
    res.discr <- simdata_discr(N=N, 
                               a=x$a, 
                               f1=x$f1, 
                               Q=x$Q, 
                               f=x$f, 
                               b=x$b, 
                               mu0=x$mu0, 
                               theta=x$theta, 
                               ystart=ystart, 
                               tstart=tstart, tend=tend, 
                               dt=dt)
    
    # Statistics
    stat <- list()
    # Age-specific means:
    bins<-10
    cutpoints<-quantile(res.discr[,3],(0:bins)/bins)
    binned <-cut(res.discr[,3],cutpoints,include.lowest=TRUE)
    for(i in seq(5,length(colnames(res.discr)),by=2)) {
      mean.cov <- tapply(res.discr[,5], binned, mean)
      stat[["mean.by.age"]][[colnames(res.discr)[i]]] <- mean.cov
    }
    
    #Survival probabilities:
    srv.prob <- survfit( Surv(t1, xi) ~ 1, conf.type="none", data = data.frame(res.discr))
    stat[["srv.prob"]] <- srv.prob
    
    res <- list(data=res.discr, stat=stat)
    
  } else if(model == "continuous") {
    
    if(length(x$a) != length(ystart)^2) {
      stop("Number of dimensions does not match with the number of values provided in ystart.")
    }
    
    # Data simulation for discrete and continuous models
    res.cont <- simdata_cont2(N=N, 
                             a=x$a, 
                             f1=x$f1, 
                             Q=x$Q, 
                             f=x$f, 
                             b=x$b, 
                             mu0=x$mu0, 
                             theta=x$theta,
                             dt=dt, 
                             ystart=ystart,
                             tstart=tstart, tend=tend, 
                             sd0=sd0)
    
    
    # Statistics
    stat <- list()
    # Age-specific means:
    bins<-10
    cutpoints<-quantile(res.cont[,3],(0:bins)/bins)
    binned <-cut(res.cont[,3],cutpoints,include.lowest=TRUE)
    for(i in seq(5,length(colnames(res.cont)),by=2)) {
      mean.cov <- tapply(res.cont[,5], binned, mean)
      stat[["mean.by.age"]][[colnames(res.cont)[i]]] <- mean.cov
    }
    
    #Survival probabilities:
    srv.prob <- survfit( Surv(t1, xi) ~ 1, conf.type="none", data = data.frame(res.cont))
    stat[["srv.prob"]] <- srv.prob
    
    res <- list(data=res.cont, stat=stat)
  }
  class(res) <- "spm.projection"
  invisible(res)
}