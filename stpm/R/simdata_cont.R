#' Multi-dimensional simulation function for continuous trait.
#'@references Yashin, A.I. et al (2007). Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.<DOI:10.1016/j.mbs.2006.11.006>.
#' @param N Number of individuals.
#' @param a A k by k matrix, which characterize the rate of the adaptive response.
#' @param f1 A particular state, which if a deviation from the normal (or optimal). This is a vector with length of k.
#' @param Q A matrix k by k, which is a non-negative-definite symmetric matrix.
#' @param f A vector-function (with length k) of the normal (or optimal) state.
#' @param b A diffusion coefficient, k by k matrix.
#' @param mu0 mortality at start period of time.
#' @param theta A displacement coefficient of the Gompertz function.
#' @param ystart A vector with length equal to number of dimensions used, defines starting values of covariates.
#' @param dt A discrete step size between two observations. A random uniform value is then added to this step size.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param sd0 a standard deviation for modelling the next covariate value.
#' @return A table with simulated data.
#' @examples
#' library(stpm)
#' dat <- simdata_cont(N=50)
#' head(dat)
#'
simdata_cont <- function(N=100, 
                         a=-0.05, 
                         f1=80, 
                         Q=2e-08, 
                         f=80, 
                         b=5, 
                         mu0=2e-05, 
                         theta=0.08,
                         dt=1, 
                         tstart=30, 
                         tend=105, 
                         ystart=80, 
                         sd0=2) {
  
  
  k <- length(ystart)
  
  if ( (dim(as.data.frame(a))[1] != k) & (dim(as.data.frame(a))[2] != k) &
       (dim(as.data.frame(Q))[1] != k) & (dim(as.data.frame(Q))[2] != k) & 
       (dim(as.data.frame(f1))[1] != k) & (dim(as.data.frame(f))[1] != k) &
       (dim(as.data.frame(b))[1] != k) & 
       (dim(as.data.frame(ystart))[1] != k) ) {
    stop("Dimenstions of provided parameters are not equal.")
  }  
  
  aH<-matrix(a,nrow=k,ncol=k)
  f1H<-matrix(f1,nrow=k,ncol=1,byrow=FALSE)
  QH<-matrix(Q,nrow=k,ncol=k)
  fH<-matrix(f,nrow=k,ncol=1,byrow=FALSE)
  bH<-matrix(b,nrow=k,ncol=1,byrow=FALSE)
  ystart<-matrix(ystart,nrow=k,ncol=1,byrow=FALSE)
  
  Q <- function(t) {
    Q <- QH*exp(theta*t)
    Q
  }
    
  mu <- function(t, par) {
    hfH <- fH - par[[1]]
    hf1H <- f1H - par[[1]]
      
    mu0Ht <- mu0*exp(theta*t)
    QH_gamma1 <- QH %*% par[[2]]
    mu <- mu0Ht + t(hfH) %*% QH %*% hfH + sum(diag(QH_gamma1))
    mu
  }
  
  
  func1 <- function(t, y) {
    hfH <- fH - y[[1]]
    hf1H <- f1H - y[[1]]
    dm <- -1.0*aH%*%hf1H + 2.0*y[[2]]%*%Q(t)%*%hfH
    dgamma <- aH%*%y[[2]] + y[[2]]%*%t(aH) + bH%*%t(bH) - 2.0*y[[2]]%*%Q(t)%*%y[[2]] 
    
    list(dm, dgamma)
  }
  
  data <- matrix(nrow=1,ncol=(4+2*k),NA)
  record <- 1
  id <- 1
  for(i in 1:N) {
    t2 <- runif(1,tstart, tend) # Starting time
    # Starting point
    new_person <- FALSE
    y2 <- matrix(unlist(lapply(1:k, function(n) {rnorm(1,mean=ystart[n,1], sd=sd0[n])} )), 
                 nrow=k, ncol=1, byrow=T)
    #out <- list(y2,matrix(0,nrow=k,ncol=k))
    yfin <- list()
    ytmp <- list()
    
    while(new_person == FALSE){
      t1 <- t2
      t2 <- t1 + dt + runif(1,0,1)
      y1 <- y2
        
      nsteps <- 2
      h=(t2-t1)/nsteps
        
      # Integration:
      s <- h/3*(-1)*mu(t1, list(y1,matrix(0,nrow=k,ncol=k)))
        
      t <- t1
      out <- list(y1,matrix(0,nrow=k,ncol=k))
      for(j in 1:nsteps) {
        # Runge-Kutta method:
        k1ar <- func1(t,out)
        yfin[[1]] <- out[[1]] + h/6.00*k1ar[[1]]
        yfin[[2]] <- out[[2]] + h/6.00*k1ar[[2]]
        ytmp[[1]] <- out[[1]] + h/2.00*k1ar[[1]]
        ytmp[[2]] <- out[[2]] + h/2.00*k1ar[[2]]
          
        k2ar <- func1(t,ytmp)
        yfin[[1]] <- yfin[[1]] + h/3.00*k2ar[[1]]
        yfin[[2]] <- yfin[[2]] + h/3.00*k2ar[[2]]
        ytmp[[1]] <- out[[1]] + h/2.00*k2ar[[1]]
        ytmp[[2]] <- out[[2]] + h/2.00*k2ar[[2]]
          
        k3ar <- func1(t,ytmp)
        yfin[[1]] <- yfin[[1]] + h/3.00*k3ar[[1]]
        yfin[[2]] <- yfin[[2]] + h/3.00*k3ar[[2]]
        ytmp[[1]] <- out[[1]] + h*k3ar[[1]]
        ytmp[[2]] <- out[[2]] + h*k3ar[[2]]
          
        k4ar <- func1(t,ytmp)
        out[[1]] <- yfin[[1]] + h/6.00*k4ar[[1]]
        out[[2]] <- yfin[[2]] + h/6.00*k4ar[[2]]
          
        t <- t + h
          
        # Integration:
        if (j == nsteps) {
          ifactor <- 1
        } else {
          if ((j %% 2) == 0) {
            ifactor <- 2
          } else {
            ifactor <- 4
          }
        }
        s <- s + ifactor*h/3.00*(-1)*mu(t,out)
          
      }
        
      m <- out[[1]]
      gamma <- out[[2]]
      
      #S <- exp(s*(t2-t1))
      S <- exp(s)
        
      xi <- 0
      if (S > runif(1,0,1)) {
        xi <- 0
        y2 <- matrix(unlist(lapply(1:k, function(n) {rnorm(1,mean=m[n,1], sd=sqrt(gamma[n,n]))} )), 
                     nrow=k, ncol=1, byrow=FALSE)
        #print(paste(m[1,1], gamma[1,1])) 
        #y2 <- matrix(unlist(lapply(1:k, function(n) {rnorm(1,mean=ystart[n,1], sd=sd0[n])} )), 
        #             nrow=k, ncol=1, byrow=FALSE)
        
        new_person <- FALSE
        cov <- unlist(lapply(seq(1,k), function(n) {c(y1[n,1], y2[n,1])}))
        data <- rbind(data, c(id, xi, t1, t2, cov))
        record <- record + 1
      } else {
        xi <- 1
        y2 <- matrix(nrow=k, ncol=1, NA)
        new_person <- TRUE
        cov <- unlist(lapply(1:k, function(n) {c(y1[n,1], y2[n,1])}))
        data <- rbind(data, c(id, xi, t1, t2, cov))
        record <- record + 1
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
  colnames(data) <- c("id","xi","t1","t2", unlist(lapply(1:k, function(n) {c(paste("y", n, sep=""), paste("y", n, ".next",sep="") )} )) )
  rownames(data) <- 1:dim(data)[1]
  invisible(data)
}

#' Multi-dimensional simulation function for continuous trait.
#' Similar to simdata_cont(...) but much faster.
#'@references Yashin, A.I. et al (2007). Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.<DOI:10.1016/j.mbs.2006.11.006>.
#' @param N Number of individuals.
#' @param a A k by k matrix, which characterize the rate of the adaptive response.
#' @param f1 A particular state, which if a deviation from the normal (or optimal). This is a vector with length of k.
#' @param Q A matrix k by k, which is a non-negative-definite symmetric matrix.
#' @param f A vector-function (with length k) of the normal (or optimal) state.
#' @param b A diffusion coefficient, k by k matrix.
#' @param mu0 mortality at start period of time.
#' @param theta A displacement coefficient of the Gompertz function.
#' @param ystart A vector with length equal to number of dimensions used, defines starting values of covariates.
#' @param dt A discrete step size between two observations. A random uniform value is then added to this step size.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param sd0 a standard deviation for modelling the next covariate value.
#' @return A table with simulated data.
#' @examples
#' library(stpm)
#' dat <- simdata_cont2(N=50)
#' head(dat)
#'
simdata_cont2 <- function(N=10, 
                          a=-0.05, 
                          f1=80, 
                          Q=2e-8, 
                          f=80, b=5, 
                          mu0=1e-5, 
                          theta=0.08, 
                          ystart=80, 
                          tstart=30, 
                          tend=105, 
                          dt=1, 
                          sd0=1) {
  
  k <- length(ystart)
  
  if ( (dim(as.data.frame(a))[1] != k) & (dim(as.data.frame(a))[2] != k) &
       (dim(as.data.frame(Q))[1] != k) & (dim(as.data.frame(Q))[2] != k) & 
       (dim(as.data.frame(f1))[1] != k) & (dim(as.data.frame(f))[1] != k) &
       (dim(as.data.frame(b))[1] != k) & 
       (dim(as.data.frame(ystart))[1] != k) ) {
    stop("Dimenstions of provided parameters are not equal.")
  }  
  
  if ( (dim(as.data.frame(a))[1] != k) & (dim(as.data.frame(a))[2] != k) &
       (dim(as.data.frame(Q))[1] != k) & (dim(as.data.frame(Q))[2] != k) & 
       (dim(as.data.frame(f1))[1] != k) & (dim(as.data.frame(f))[1] != k) &
       (dim(as.data.frame(b))[1] != k) & 
       (dim(as.data.frame(ystart))[1] != k) ) {
    stop("Dimenstions of provided parameters are not equal.")
  }  
  
  aH<-matrix(a,nrow=k,ncol=k)
  f1H<-matrix(f1,nrow=k,ncol=1,byrow=FALSE)
  QH<-matrix(Q,nrow=k,ncol=k)
  fH<-matrix(f,nrow=k,ncol=1,byrow=FALSE)
  bH<-matrix(b,nrow=k,ncol=1,byrow=FALSE)
  ystart<-matrix(ystart,nrow=k,ncol=1,byrow=FALSE)
  
  simulated = .Call("simCont", N, aH, f1H, QH, fH, bH, mu0, theta, tstart, ystart, tend, k, dt, sd0);
  
  data_names <- c()
  for(n in 1:k) {
    data_names <- c(data_names, paste("y",n, sep=''), paste("y",n, ".next", sep=''))
  }
  colnames(simulated) <- c("id", "xi", "t1", "t2", data_names)
  
  invisible(simulated)
}

