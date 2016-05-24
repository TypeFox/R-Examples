#' Simulation function for continuous trait with time-dependant coefficients.
#' @param N Number of individuals.
#' @param f a list of formulas that define age (time) - dependency. Default: list(at="a", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)")
#' @param step An interval between two observations, a random uniformally-distributed value is then added to this step.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param ystart A starting value of covariates.
#' @param sd0 A standard deviation for modelling the next covariate value.
#' @return A table with simulated data.
#'@references Yashin, A. et al (2007), Health decline, aging and mortality: how are they related? 
#'Biogerontology, 8(3), 291-302.<DOI:10.1007/s10522-006-9073-3>.
#' @examples
#' library(stpm)
#' dat <- simdata_time_dep(N=100)
#' head(dat)
#'
simdata_time_dep <- function(N=10,f=list(at="-0.05", f1t="80", Qt="2e-5", ft="80", bt="2.5", mu0t="1e-1"),
                         step=1, tstart=30, tend=105, ystart=80, sd0=2) {
  formulas <- f  
  at <- NULL
  f1t <- NULL
  Qt <- NULL
  ft <- NULL
  bt <- NULL
  mu0t <- NULL
  
  comp_func_params <- function(astring, f1string, qstring, fstring, bstring, mu0string) {
    at <<- eval(bquote(function(t) .(parse(text = astring)[[1]])))
    f1t <<- eval(bquote(function(t) .(parse(text = f1string)[[1]]))) 
    Qt <<- eval(bquote(function(t) .(parse(text = qstring)[[1]])))
    ft <<- eval(bquote(function(t) .(parse(text = fstring)[[1]])))
    bt <<- eval(bquote(function(t) .(parse(text = bstring)[[1]])))
    mu0t <<- eval(bquote(function(t) .(parse(text = mu0string)[[1]])))
  }
  
  sigma_sq <- function(t1, t2) {
    # t2 = t_{j}, t1 = t_{j-1}
    ans <- bt(t1)*(t2-t1)
    #ans <- 2.5*(t2-t1)
    ans
  }
  
  m <- function(y, t1, t2) {
    # y = y_{j-1}, t1 = t_{j-1}, t2 = t_{j}
    ans <- y + at(t1)*(y - f1t(t1))*(t2 - t1)
    #ans <- y + (-0.05)*(y - 80)*(t2 - t1)
    ans
  }
  
  mu <- function(y, t) {
    ans <- mu0t(t) + (y - ft(t))^2*Qt(t)
    #ans <- 1e-1 + (y - 80)^2*2e-5
  }
  
  
  comp_func_params(formulas$at, formulas$f1t, formulas$Qt, formulas$ft, formulas$bt, formulas$mu0t)
  
  data <- matrix(nrow=1,ncol=6, NA)
  record <- 1
  id <- 1
  for(i in 1:N) {
    t2 <- runif(1,tstart, tend) # Starting time
    # Starting point
    new_person <- FALSE
    y2 <- rnorm(1,mean=ystart, sd=sd0)  
    while(new_person == FALSE){
      t1 <- t2
      t2 <- t1 + runif(1,0,1) + step
      y1 <- y2
        
      S <- exp(-1*mu(y1,t1)*(t2-t1))
      #S <- exp(-1*mu(y1, t1))
      xi <- 0
      if (S > runif(1,0,1)) {
        xi <- 0
        y2 <- rnorm(1,mean=m(y1, t1, t2), sd=sqrt(sigma_sq(t1,t2)))
        #y2 <- rnorm(1,mean=m(y1, t1, t2), sd=sigma_sq(t1,t2))
        new_person <- FALSE
        cov <- c(y1, y2)
        data <- rbind(data, c(id, xi, t1, t2, cov))
        
      } else {
        xi <- 1
        y2 <- NA
        new_person <- TRUE
        cov <- c(y1, y2)
        data <- rbind(data, c(id, xi, t1, t2, cov))
        id <- id + 1
        
      }
        
      if(t2 > tend & new_person == FALSE) {
        new_person <- TRUE
        id <- id + 1
      }
    }
      
  }
    
  # One last step:
  data <- data[2:dim(data)[1],]
  colnames(data) <- c("id","xi","t1","t2", "y", "y.next")
  invisible(data)
}