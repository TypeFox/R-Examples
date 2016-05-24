########################################################################
# Bayesdouble.R - code to perform Bayesian analysis for comparing two
# bacterial growth curves
# (c) Lydia Rickett
########################################################################

.logllDouble <- function(params,par.no,modelfunc,model,dataset1,dataset2,inc.nd1,inc.nd2,threshold1,threshold2,t.nd1,t.nd2,inf.sigma1,
                  inf.sigma2,sigma1,sigma2,hyp,transformParams,mumax.prior1,mumax.prior2) {
  
  t1 <- dataset1$t1; y1 <- dataset1$y1
  t2 <- dataset2$t2; y2 <- dataset2$y2
  ymax1 <- max(y1); ymin1 <- min(y1); tmax1 <- max(t1); ydiff1 <- (ymax1-ymin1)
  ymax2 <- max(y2); ymin2 <- min(y2); tmax2 <- max(t2); ydiff2 <- (ymax2-ymin2)
  tParams <- transformParams(params)
  
  params1 <- tParams[1:par.no]
  if (hyp == "H1"){
    params2 <- params1
  } else if (hyp == "H2"){
    if (par.no == 4){
      params2 <- c(tParams[5],tParams[6],tParams[3],tParams[7])
    } else if (model == "logistic"){
      params2 <- c(tParams[4],tParams[5],tParams[3])
    } else if (model == "Bar3par"){
      params2 <- c(tParams[4],tParams[2],tParams[5])
    } else { # par.no == 2
      params2 <- c(tParams[3],tParams[2]) 
    }
  } else if (hyp == "H3"){
    params2 <- tParams[(par.no+1):(2*par.no)]
  }
  
  if (inf.sigma1 && inf.sigma2){
    sigma1 <- rep(tParams[(length(tParams)-1)],length=length(y1))
    sigma2 <- rep(tParams[length(tParams)],length=length(y2))
  } else if (inf.sigma1) {
    sigma1 <- rep(tParams[length(tParams)],length=length(y1))
  } else if (inf.sigma2) {
    sigma2 <- rep(tParams[length(tParams)],length=length(y2))
  }
  
  # If using a Cauchy prior for mumax, penalise values outside of a sensible range in the
  # likelihood function (to prevent guesses that are negative or too large)
  if ( ((mumax.prior1 == "Cauchy") && (model == "linear" || model == "Bar3par")
      && (tParams[2] < 0 || tParams[2] > 10 * ydiff1/tmax1)) ||
        ((mumax.prior1 == "Cauchy") && (model == "logistic" || model == "Bar4par") 
         && (tParams[3] < 0 || tParams[3] > 10 * ydiff1/tmax1)) ||
        ((hyp == "H3") && (mumax.prior2 == "Cauchy") && (model == "linear")
         && (tParams[4] < 0 || tParams[4] > 10 * ydiff2/tmax2)) ||
        ((hyp == "H3") && (mumax.prior2 == "Cauchy") && (model == "Bar3par")
         && (tParams[5] < 0 || tParams[5] > 10 * ydiff2/tmax2)) ||
        ((hyp == "H3") && (mumax.prior2 == "Cauchy") && (model == "logistic") 
         && (tParams[6] < 0 || tParams[6] > 10 * ydiff2/tmax2)) ||
        ((hyp == "H3") && (mumax.prior2 == "Cauchy") && (model == "Bar4par") 
         && (tParams[7] < 0 || tParams[7] > 10 * ydiff2/tmax2)) ) {
          llhood <- -10^4
  } else {
    fun1 <- modelfunc(t1,params1); fun2 <- modelfunc(t2,params2)
    res1 <- sum(((y1-fun1)/sigma1)^2); res2 <- sum(((y2-fun2)/sigma2)^2)
    llhood <- log(1/prod(sigma1*(2*pi)^(1/2))) + log(1/prod(sigma2*(2*pi)^(1/2))) - res1/2 - res2/2
    
    if (inc.nd1 && length(t.nd1 > 0)) {
      fun.nd1 <- modelfunc(t.nd1,params1)
      # Approximate censored values using a uniform distribution
      P.nd1 <- rep(1/threshold1,length=length(fun.nd1))
      P.d1 <- numeric(length=length(fun.nd1))
      for (k in 1:length(t.nd1)){
        if ((fun.nd1[k] >= 0) && (fun.nd1[k] <= threshold1)) {
          P.d1[k] <- 0 # In censored region, uniform distribution
        } else {         
          # In uncensored region, add Gaussian tail with very small variance to avoid negative infinity in log likelihood
          P.d1[k] <- ((fun.nd1[k]-threshold1)/10^(-10))^2 
        }
      }
      llhood <- llhood + log(prod(P.nd1)) - sum(P.d1)
    }
    if (inc.nd2 && length(t.nd2 > 0)) {
      fun.nd2 <- modelfunc(t.nd2,params2)
      # Approximate censored values using a uniform distribution
      P.nd2 <- rep(1/threshold2,length=length(fun.nd2))
      P.d2 <- numeric(length=length(fun.nd2))
      for (k in 1:length(t.nd2)){
        if ((fun.nd2[k] >= 0) && (fun.nd2[k] <= threshold2)) {
          P.d2[k] <- 0 # In censored region, uniform distribution
        } else {         
          # In uncensored region, add Gaussian tail with very small variance to avoid negative infinity in log likelihood
          P.d2[k] <- ((fun.nd2[k]-threshold2)/10^(-10))^2 
        }
      }
      llhood <- llhood + log(prod(P.nd2)) - sum(P.d2)
    }
  }
  
  return (llhood)
}

.generateTransformDouble <- function(dataset1,dataset2,par.no,model,inc.nd1,inc.nd2,inf.sigma1,inf.sigma2,hyp,
                              mumax.prior1,mu.mean1,mu.sd1,mumax.prior2,mu.mean2,mu.sd2) {
  
  t1 <- dataset1$t1; t2 <- dataset2$t2; y1 <- dataset1$y1; y2 <- dataset2$y2
  ymax1 <- max(y1); ymin1 <- min(y1); tmax1 <- max(t1); ydiff1 <- (ymax1-ymin1)
  ymax2 <- max(y2); ymin2 <- min(y2); tmax2 <- max(t2); ydiff2 <- (ymax2-ymin2)

  transformParams <- function(uParams) {
    
    tParams = numeric(length = length(uParams))
    
    # when t0!=0 or including censored values, use lower bound y0=0 instead, or lower bound might not be low enough
    if (t1[1] != 0 | inc.nd1){
      fact1 <- 0
    } else {
      fact1 <- 1
    }
    if (t2[1] != 0 | inc.nd2){
      fact2 <- 0
    } else {
      fact2 <- 1
    }
    
    if (hyp == "H1"){
      
      if (par.no == 4) {
        tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
        tParams[2] = UniformPrior(uParams[2], ymax1 - ydiff1/2, ymax1 + ydiff1/2)
        if (mumax.prior1 == "Gaussian") {
          tParams[3] = GaussianPrior(uParams[3], mu.mean1, mu.sd1)
        } else if (mumax.prior1 == "Cauchy") {
          tParams[3] = CauchyPrior(uParams[3], mu.mean1, mu.sd1)
        } else { # mumax.prior1 = Uniform
          tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff1/tmax1)
        }
        tParams[4] = UniformPrior(uParams[4], 0, 9 * ydiff1)
      } else if (par.no == 3) {
        if (model == "logistic"){
          tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
          tParams[2] = UniformPrior(uParams[2], ymax1 - ydiff1/2, ymax1 + ydiff1/2)
          if (mumax.prior1 == "Gaussian") {
            tParams[3] = GaussianPrior(uParams[3], mu.mean1, mu.sd1)
          } else if (mumax.prior1 == "Cauchy") {
            tParams[3] = CauchyPrior(uParams[3], mu.mean1, mu.sd1)
          } else { # mumax.prior1 = Uniform
            tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff1/tmax1)
          }
        } else {
          tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
          if (mumax.prior1 == "Gaussian") {
            tParams[2] = GaussianPrior(uParams[2], mu.mean1, mu.sd1)
          } else if (mumax.prior1 == "Cauchy") {
            tParams[2] = CauchyPrior(uParams[2], mu.mean1, mu.sd1)
          } else { # mumax.prior1 = Uniform
            tParams[2] = UniformPrior(uParams[2], 0, 10 * ydiff1/tmax1)
          }
          tParams[3] = UniformPrior(uParams[3], 0, 9 * ydiff1)
        }
      } else { # par.no = 2
        tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
        if (mumax.prior1 == "Gaussian") {
          tParams[2] = GaussianPrior(uParams[2], mu.mean1, mu.sd1)
        } else if (mumax.prior1 == "Cauchy") {
          tParams[2] = CauchyPrior(uParams[2], mu.mean1, mu.sd1)
        } else { # mumax.prior1 = Uniform
          tParams[2] = UniformPrior(uParams[2], 0, 10 * ydiff1/tmax1)
        }
      }
      
    } else if (hyp == "H2"){
      
      if (par.no == 4) { # params = (y01,ymax1,mu_max,h_01,y02,ymax2,h_02)
        tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
        tParams[2] = UniformPrior(uParams[2], ymax1 - ydiff1/2, ymax1 + ydiff1/2)
        if (mumax.prior1 == "Gaussian") {
          tParams[3] = GaussianPrior(uParams[3], mu.mean1, mu.sd1)
        } else if (mumax.prior1 == "Cauchy") {
          tParams[3] = CauchyPrior(uParams[3], mu.mean1, mu.sd1)
        } else { # mumax.prior1 = Uniform
          tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff1/tmax1)
        }
        tParams[4] = UniformPrior(uParams[4], 0, 9 * ydiff1)
        tParams[5] = UniformPrior(uParams[5], fact2*(ymin2 - ydiff2/2), ymin2 + ydiff2/2)
        tParams[6] = UniformPrior(uParams[6], ymax2 - ydiff2/2, ymax2 + ydiff2/2)
        tParams[7] = UniformPrior(uParams[7], 0, 9 * ydiff2)
      } else if (par.no == 3) { # params = (y01,ymax1,mu_max,y02,ymax2)
        if (model == "logistic"){
          tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
          tParams[2] = UniformPrior(uParams[2], ymax1 - ydiff1/2, ymax1 + ydiff1/2)
          if (mumax.prior1 == "Gaussian") {
            tParams[3] = GaussianPrior(uParams[3], mu.mean1, mu.sd1)
          } else if (mumax.prior1 == "Cauchy") {
            tParams[3] = CauchyPrior(uParams[3], mu.mean1, mu.sd1)
          } else { # mumax.prior1 = Uniform
            tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff1/tmax1)
          }
          tParams[4] = UniformPrior(uParams[4], fact2*(ymin2 - ydiff2/2), ymin2 + ydiff2/2)
          tParams[5] = UniformPrior(uParams[5], ymax2 - ydiff2/2, ymax2 + ydiff2/2)
        } else { # params = (y01,mu_max,h_01,y02,h_01)
          tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
          if (mumax.prior1 == "Gaussian") {
            tParams[2] = GaussianPrior(uParams[2], mu.mean1, mu.sd1)
          } else if (mumax.prior1 == "Cauchy") {
            tParams[2] = CauchyPrior(uParams[2], mu.mean1, mu.sd1)
          } else { # mumax.prior1 = Uniform
            tParams[2] = UniformPrior(uParams[2], 0, 10 * ydiff1/tmax1)
          }
          tParams[3] = UniformPrior(uParams[3], 0, 9 * ydiff1)
          tParams[4] = UniformPrior(uParams[4], fact2*(ymin2 - ydiff2/2), ymin2 + ydiff2/2)
          tParams[5] = UniformPrior(uParams[5], 0, 9 * ydiff2)
        }
      } else { # par.no = 2, params = (y01,mumax,y02)
        tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
        if (mumax.prior1 == "Gaussian") {
          tParams[2] = GaussianPrior(uParams[2], mu.mean1, mu.sd1)
        } else if (mumax.prior1 == "Cauchy") {
          tParams[2] = CauchyPrior(uParams[2], mu.mean1, mu.sd1)
        } else { # mumax.prior1 = Uniform
          tParams[2] = UniformPrior(uParams[2], 0, 10 * ydiff1/tmax1)
        }
        tParams[3] = UniformPrior(uParams[3], fact2*(ymin2 - ydiff2/2), ymin2 + ydiff2/2)
      }
      
    } else{ # hyp = H3
      
      if (par.no == 4) { # params = (y01,ymax1,mu_max1,h_01,y02,ymax2,mu_max2,h_02)
        tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
        tParams[2] = UniformPrior(uParams[2], ymax1 - ydiff1/2, ymax1 + ydiff1/2)
        if (mumax.prior1 == "Gaussian") {
          tParams[3] = GaussianPrior(uParams[3], mu.mean1, mu.sd1)
        } else if (mumax.prior1 == "Cauchy") {
          tParams[3] = CauchyPrior(uParams[3], mu.mean1, mu.sd1)
        } else { # mumax.prior1 = Uniform
          tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff1/tmax1)
        }
        tParams[4] = UniformPrior(uParams[4], 0, 9 * ydiff1)
        tParams[5] = UniformPrior(uParams[5], fact2*(ymin2 - ydiff2/2), ymin2 + ydiff2/2)
        tParams[6] = UniformPrior(uParams[6], ymax2 - ydiff2/2, ymax2 + ydiff2/2)
        if (mumax.prior2 == "Gaussian") {
          tParams[7] = GaussianPrior(uParams[7], mu.mean2, mu.sd2)
        } else if (mumax.prior2 == "Cauchy") {
          tParams[7] = CauchyPrior(uParams[7], mu.mean2, mu.sd2)
        } else { # mumax.prior2 = Uniform
          tParams[7] = UniformPrior(uParams[7], 0, 10 * ydiff2/tmax2)
        }
        tParams[8] = UniformPrior(uParams[8], 0, 9 * ydiff2)
      } else if (par.no == 3) {
        if (model == "logistic"){ # params = (y01,ymax1,mu_max1,y02,ymax2,mu_max2)
          tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
          tParams[2] = UniformPrior(uParams[2], ymax1 - ydiff1/2, ymax1 + ydiff1/2)
          if (mumax.prior1 == "Gaussian") {
            tParams[3] = GaussianPrior(uParams[3], mu.mean1, mu.sd1)
          } else if (mumax.prior1 == "Cauchy") {
            tParams[3] = CauchyPrior(uParams[3], mu.mean1, mu.sd1)
          } else { # mumax.prior1 = Uniform
            tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff1/tmax1)
          }
          tParams[4] = UniformPrior(uParams[4], fact2*(ymin2 - ydiff2/2), ymin2 + ydiff2/2)
          tParams[5] = UniformPrior(uParams[5], ymax2 - ydiff2/2, ymax2 + ydiff2/2)
          if (mumax.prior2 == "Gaussian") {
            tParams[6] = GaussianPrior(uParams[6], mu.mean2, mu.sd2)
          } else if (mumax.prior2 == "Cauchy") {
            tParams[6] = CauchyPrior(uParams[6], mu.mean2, mu.sd2)
          } else { # mumax.prior2 = Uniform
            tParams[6] = UniformPrior(uParams[6], 0, 10 * ydiff2/tmax2)
          }
        } else { # params = (y01,mu_max1,h_01,y02,mu_max1,h_01)
          tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
          if (mumax.prior1 == "Gaussian") {
            tParams[2] = GaussianPrior(uParams[2], mu.mean1, mu.sd1)
          } else if (mumax.prior1 == "Cauchy") {
            tParams[2] = CauchyPrior(uParams[2], mu.mean1, mu.sd1)
          } else { # mumax.prior1 = Uniform
            tParams[2] = UniformPrior(uParams[2], 0, 10 * ydiff1/tmax1)
          }
          tParams[3] = UniformPrior(uParams[3], 0, 9 * ydiff1)
          tParams[4] = UniformPrior(uParams[4], fact2*(ymin2 - ydiff2/2), ymin2 + ydiff2/2)
          if (mumax.prior2 == "Gaussian") {
            tParams[5] = GaussianPrior(uParams[5], mu.mean2, mu.sd2)
          } else if (mumax.prior2 == "Cauchy") {
            tParams[5] = CauchyPrior(uParams[5], mu.mean2, mu.sd2)
          } else { # mumax.prior2 = Uniform
            tParams[5] = UniformPrior(uParams[5], 0, 10 * ydiff2/tmax2)
          }
          tParams[6] = UniformPrior(uParams[6], 0, 9 * ydiff2)
        }
      } else { # par.no = 2, params = (y01,mumax1,y02,mu_max2)
        tParams[1] = UniformPrior(uParams[1], fact1*(ymin1 - ydiff1/2), ymin1 + ydiff1/2)
        if (mumax.prior1 == "Gaussian") {
          tParams[2] = GaussianPrior(uParams[2], mu.mean1, mu.sd1)
        } else if (mumax.prior1 == "Cauchy") {
          tParams[2] = CauchyPrior(uParams[2], mu.mean1, mu.sd1)
        } else { # mumax.prior1 = Uniform
          tParams[2] = UniformPrior(uParams[2], 0, 10 * ydiff1/tmax1)
        }
        tParams[3] = UniformPrior(uParams[3], fact2*(ymin2 - ydiff2/2), ymin2 + ydiff2/2)
        if (mumax.prior2 == "Gaussian") {
          tParams[4] = GaussianPrior(uParams[4], mu.mean2, mu.sd2)
        } else if (mumax.prior2 == "Cauchy") {
          tParams[4] = CauchyPrior(uParams[4], mu.mean2, mu.sd2)
        } else { # mumax.prior2 = Uniform
          tParams[4] = UniformPrior(uParams[4], 0, 10 * ydiff2/tmax2)
        }
      }
      
    }
    
    if (inf.sigma1 && inf.sigma2) { # Extra parameter for inferring sigma - bounds between 0.01 and 10
      tParams[(length(uParams)-1)] = JeffreysPrior(uParams[(length(uParams)-1)], -2, 1)
      tParams[length(uParams)] = JeffreysPrior(uParams[length(uParams)], -2, 1)
    } else {
      if (inf.sigma1 | inf.sigma2){
        tParams[length(uParams)] = JeffreysPrior(uParams[length(uParams)], -2, 1)
      }
    }
    
    return(tParams)
  }
  
  return(transformParams)
}

###
# Functions for growth models:
###

.linear <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the linear model (2-parameters)
  #
  # Arguments: vector t of times and parameters par1 = y0 and par3 = mumax
  
  par1 <- params[1]; par3 <- params[2];
  
  y <- par1 + par3*t
  
  return (y=y)
}

.logistic <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the logistic model (3-parameters, no lag phase)
  #
  # Arguments: vector t of times and parameters par1 = y0, par2 = ymax and par3 = mumax

  par1 <- params[1]; par2 <- params[2]; par3 <- params[3]
  
  Y2 <- rep(1,length(t))
  for (i in 1:length(t)) { 
    if (par3*t[i] > 700) { # Preventing overflow
      Y2[i] <- 0
    }
  }
  y <- par1 + Y2*(par3*t - log(1+(exp(par3*t)-1)/exp(par2-par1)))
  
  return (y=y)
}

.Bar3par <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the 3-parameter Baranyi model (no stationary phase)
  #
  # Arguments: vector t of times and parameters par1 = y0, par3 = mumax and par4 = h0
  
  par1 <- params[1]; par3 <- params[2]; par4 <- params[3]
  
  a1 <- rep(1,length(t)); a2 <- rep(1,length(t)); a3 <- rep(1,length(t))
  for (i in 1:length(t)) {
    if (abs(par3*t[i]-par4) > 700) { # Preventing over/underflow
      if (-(par3*t[i]-par4) > 0) {
        a1[i] <- 0
      } else {
        a3[i] <- 0
      }
    }
    if (-par3*t[i] < -700) {
      a2[i] <- 0
    }
  }
  A <- a1*(t-par4/par3 + log(1-a2*exp(-par3*t)+a3*exp(-(par3*t-par4)))/par3)
  
  Y1 <- rep(1,length(t))
  for (i in 1:length(t)) { 
    if (par3*A[i] < 0) { # Preventing round-off error; mumax*A(t) must be positive
      Y1[i] <- 0
    }
  }
  muA <- Y1*par3*A
  y <- par1 + muA
  
  return (y=y)
}

.Bar4par <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the 4-parameter Baranyi model
  #
  # Arguments: vector t of times and parameters par1 = y0, par2 = ymax, par3 = mumax and par4 = h0
  
  par1 <- params[1]; par2 <- params[2]; par3 <- params[3]; par4 <- params[4]
  
  a1 <- rep(1,length(t)); a2 <- rep(1,length(t)); a3 <- rep(1,length(t))
  for (i in 1:length(t)) {
    if (abs(par3*t[i]-par4) > 700) { # Preventing over/underflow
      if (-(par3*t[i]-par4) > 0) {
        a1[i] <- 0
      } else {
        a3[i] <- 0
      }
    }
    if (-par3*t[i] < -700) {
      a2[i] <- 0
    }
  }
  A <- a1*(t-par4/par3 + log(1-a2*exp(-par3*t)+a3*exp(-(par3*t-par4)))/par3)
  
  Y1 <- rep(1,length(t)); Y2 <- rep(1,length(t))
  for (i in 1:length(t)) { 
    if (par3*A[i] < 0) { # Preventing round-off error; mumax*A(t) must be positive
      Y1[i] <- 0
    }
    if (par3*A[i] > 700) { # Preventing overflow
      Y2[i] <- 0
    }
  }
  muA <- Y1*par3*A
  y <- par1 + Y2*(muA - log(1+(exp(muA)-1)/exp(par2-par1)))
  
  return (y=y)
}

.sortdata <- function(data,inc.nd) {
  # Sort out data for use in nested sampling
  
  t <- numeric(0)
  y <- numeric(0)
  t.nd <- numeric(0) # Times of undetected values
  
  for (i in 1:nrow(data)) {
    if(is.na(data[i,2]) == FALSE){ # Takes care of "NA"/"ND" values
      t <- c(t,data[i,1]); y <- c(y,data[i,2])} else {
        if (inc.nd) {
          t.nd <- c(t.nd,data[i,1])
        } else {
          warning("there are undetected points which are not being used as part of the analysis")
        }
      }
  }
  
  ty <- data.frame(t,y)
  
  return(list(ty=ty,t.nd=t.nd))
  ### Returns: the dataset minus undetected values, (t,y), and list of times of undetected values, t.nd
}

Bayescompare <- structure(function
  ### Perform Bayesian analysis for comparing two bacterial growth curves using the Baranyi model.
  (
  data1,
  data2,
  ### Datafiles of the two curves to be fitted. This should consist of two columns, the first for time and second for logc.  
  ### The bacterial concentration should be given in log_10 cfu and there should be at least 2 data points (the first 
  ### of which may be undetected). Undetected y values should be represented by "NA".
  hyp,
  ### Hypothesis to test. This should be one of "H1" (data curves replicates), "H2" (data curves have same growth rate) and 
  ### "H3" (all data curve parameters are different)
  ### 
  model,
  ### The growth model to be used. This should be one of "linear", "logistic", "Bar3par" and "Bar4par".
  inf.sigma1=TRUE,
  inf.sigma2=TRUE,
  ### (TRUE/FALSE) Choose whether or not to infer the noise levels, sigma1 (for curve 1) and sigma2 (for curve 2), as part of 
  ### the analysis. If FALSE, sigma should be specified (or the default value of sigma, 0.3, will be used).
  inc.nd1=FALSE,
  inc.nd2=FALSE,
  ### Choose whether or not to include undetected points for curves 1 and 2 respectively as part of the analysis. If TRUE, 
  ### threshold should be specified.
  sigma1=0.3,
  sigma2=0.3,
  ### The choice of noise levels, sigma1 and sigma2, in log_10 cfu if not inferred as part of the analysis. Default is 0.3.
  threshold1=NULL,
  threshold2=NULL,
  ### Thresholds in log_10 cfu below which values are considered as undetected.
  mumax.prior1="Uniform",
  mumax.prior2="Uniform",
  ### The type of priors to use for mu_max1 and mu_max2. These should be one of "Uniform", "Gaussian" or "Cauchy" (or the 
  ### default "Uniform" will be used). If "Gaussian" or "Cauchy" are specified for either, mu.mean and mu.sd should be given. 
  mu.mean1=NULL,
  mu.mean2=NULL,
  ### The means to be used when using a Gaussian or Cauchy prior.
  mu.sd1=NULL,
  mu.sd2=NULL,
  ### The standard deviations to be used when using a Gaussian or Cauchy prior.
  tol=0.1,
  ### The termination tolerance for nested sampling
  prior.size=250
  ### The number of prior samples to use for nested sampling
) {

  # Setup (converting to natural log scale for use in the model)
  
  if(model!="linear" && model!="logistic" && model!="Bar3par" && model!="Bar4par") {
    stop("'model' must be one of 'linear', 'logistic', 'Bar3par' or 'Bar4par'")
  }
  
  if(mumax.prior1!="Uniform" && mumax.prior1!="Gaussian" && mumax.prior1!="Cauchy") {
    stop("'mumax.prior1' must be one of 'Uniform', 'Gaussian' or 'Cauchy'")
  }
  if(mumax.prior2!="Uniform" && mumax.prior2!="Gaussian" && mumax.prior2!="Cauchy") {
    stop("'mumax.prior2' must be one of 'Uniform', 'Gaussian' or 'Cauchy'")
  }
  if((hyp=="H1" || hyp=="H2") && mumax.prior2 != "Uniform") {
    warning("using the current hypothesis, mumax.prior1 will be used for both parameter sets and mumax.prior2 
            will not be used")
  }
  if((mumax.prior1=="Gaussian" || mumax.prior1=="Cauchy") && (is.null(mu.mean1) || is.null(mu.sd1))) {
    stop("both 'mu.mean1' and 'mu.sd1' must be given when using a Gaussian or Cauchy prior for mu_max1")
  }
  if((mumax.prior2=="Gaussian" || mumax.prior2=="Cauchy") && (is.null(mu.mean2) || is.null(mu.sd2))) {
    stop("both 'mu.mean2' and 'mu.sd2' must be given when using a Gaussian or Cauchy prior for mu_max2")
  }
  mu.mean1 <- mu.mean1*log(10); mu.sd1 <- mu.sd1*log(10)
  mu.mean2 <- mu.mean2*log(10); mu.sd2 <- mu.sd2*log(10)
  
  if (model == "linear") {
    par.no <- 2
    modelfunc <- .linear
  } else if (model == "logistic") {
    par.no <- 3
    modelfunc <- .logistic
  } else if (model == "Bar3par") {
    par.no <- 3
    modelfunc <- .Bar3par
  } else if (model == "Bar4par") {
    par.no <- 4
    modelfunc <- .Bar4par
  }
  
  if (nrow(data1) == 1) {
    stop("'data1' must include at least 2 data points")
  }
  if (nrow(data2) == 1) {
    stop("'data2' must include at least 2 data points")
  }
  
  if ((nrow(data1) == 2) && (inf.sigma1 == TRUE)) {
    warning("inferring noise level sigma is not recommended with only 2 data points, sigma1 has been prescribed as 0.3")
  }
  if ((nrow(data2) == 2) && (inf.sigma2 == TRUE)) {
    warning("inferring noise level sigma is not recommended with only 2 data points, sigma2 has been prescribed as 0.3")
  }

  sort1 <- .sortdata(data1,inc.nd1); sort2 <- .sortdata(data2,inc.nd2)
  data1 <- sort1$ty; data2 <- sort2$ty
  t1 <- data1[,1]; y1 <- log(10)*data1[,2] # Convert y from log_10 to natural log for use in model
  t2 <- data2[,1]; y2 <- log(10)*data2[,2]
  if (length(t1) ==1) {
    stop("'data1' must include at least 2 detected data points")
  }
  if (length(t2) ==1) {
    stop("'data2' must include at least 2 detected data points")
  }
  dataset1 <- data.frame(t1,y1); dataset2 <- data.frame(t2,y2)
  t.nd1 <- sort1$t.nd; t.nd2 <- sort2$t.nd
  
  if(hyp!="H1" && hyp!="H2" && hyp!="H3") {
  stop("'hyp' must be one of 'H1', 'H2' or 'H3'")
  }
  
  if (hyp == "H1"){
    hypinfo <- "data curves are replicates"
    hyp.par.no <- par.no # Calculate total number of parameters
  } else if (hyp == "H2"){
    hypinfo <- "data curves have same growth rate"
    hyp.par.no <- 2*par.no-1
  } else { # hyp = H3
    hypinfo <- "all data curve parameters are different"
    hyp.par.no <- 2*par.no
  }
  
  if (!inf.sigma1) {
    total.pars <- hyp.par.no
    sigma1.save <- sigma1
    sigma1 <- rep(sigma1*log(10),length(dataset1[,1]))
  } else {
    total.pars <- hyp.par.no + 1
  }
  if (!inf.sigma2) {
    total.pars <- total.pars
    sigma2.save <- sigma2
    sigma2 <- rep(sigma2*log(10),length(dataset2[,1]))
  } else {
    total.pars <- total.pars + 1
  }
  
  if ((inc.nd1) && (is.null(threshold1))) {
    stop("'threshold1' must be specified when including undetected values in curve 1")
  }
  if ((inc.nd2) && (is.null(threshold2))) {
    stop("'threshold2' must be specified when including undetected values in curve 2")
  }
  threshold1 <- threshold1*log(10); threshold2 <- threshold2*log(10)

  # Perform nested sampling

  #tol <- 0.1 # Set the termination tolerance

  # Define transformed priors and log likelihood function
  transformParams <- .generateTransformDouble(dataset1,dataset2,par.no,model,inc.nd1,inc.nd2,inf.sigma1,inf.sigma2,hyp,
                                       mumax.prior1,mu.mean1,mu.sd1,mumax.prior2,mu.mean2,mu.sd2)
  logllfun <- function(params) { 
    return(.logllDouble(params,par.no,modelfunc,model,dataset1,dataset2,inc.nd1,inc.nd2,threshold1,threshold2,t.nd1,t.nd2,inf.sigma1,
                 inf.sigma2,sigma1,sigma2,hyp,transformParams,mumax.prior1,mumax.prior2))
  }
  
  # Call the nested sampling routine, which returns the posterior, log evidence & error, logweights and parameter means 
  # and variances
  ret <- nestedSampling(logllfun, total.pars, prior.size, transformParams, exploreFn=ballExplore, tolerance = tol)
  posterior <- ret$posterior
  logevidence <- ret$logevidence
  means <- ret$parameterMeans
  vars <- ret$parameterVariances

  # Gather data for posterior plots
  
  # Get equally weighted samples from posterior distribution using staircase sampling
  chosenSamples = getEqualSamples(ret$posterior, n = Inf)
  
  # Fit model with posterior samples
  if (length(t1) == 1) {
    tfit1 <- seq(from=0,by=0.01*t1,to=t1)
  } else {
    tfit1 <- seq(from=t1[1],by=0.01*t1[length(t1)],to=t1[length(t1)])
  }
  if (length(t2) == 1) {
    tfit2 <- seq(from=0,by=0.01*t2,to=t2)
  } else {
    tfit2 <- seq(from=t2[1],by=0.01*t2[length(t2)],to=t2[length(t2)])
  }
  posteriorModel1 = apply(chosenSamples[, -1] , 1, modelfunc, t=tfit1)
   if (hyp == "H1") {
     posteriorModel2 = apply(chosenSamples[, -1] , 1, modelfunc, t=tfit2)
   } else if (hyp == "H2") {
     if (model == "linear") {
       posteriorModel2 = apply(chosenSamples[, c(4,3)], 1, modelfunc, t=tfit2)
     } else if (model == "logistic") {
       posteriorModel2 = apply(chosenSamples[, c(5,6,4)], 1, modelfunc, t=tfit2)
     } else if (model == "Bar3par") {
       posteriorModel2 = apply(chosenSamples[, c(5,3,6)], 1, modelfunc, t=tfit2)
     } else { # model == "Bar4par")
       posteriorModel2 = apply(chosenSamples[, c(6,7,4,8)], 1, modelfunc,t=tfit2)
     }
   } else { # hyp == H3
     posteriorModel2 = apply(chosenSamples[, (par.no+2):(2*par.no+1)], 1, modelfunc, t=tfit2)
   }

  # Fit model with mean samples
  meanfit1 <- modelfunc(tfit1,means)
  if (hyp == "H1") {
    meanfit2 <- modelfunc(tfit2,means)
  } else if (hyp == "H2") {
    if (model == "linear") {
      meanfit2 <- modelfunc(tfit2,means[c(3,2)])
    } else if (model == "logistic") {
      meanfit2 <- modelfunc(tfit2,means[c(4,5,3)])
    } else if (model == "Bar3par") {
      meanfit2 <- modelfunc(tfit2,means[c(4,2,5)])
    } else { # model == "Bar4par")
      meanfit2 <- modelfunc(tfit2,means[c(5,6,3,7)])
    }
  } else { # hyp == H3
    meanfit2 <- modelfunc(tfit2,means[(par.no+1):(2*par.no)])
  }

  # Convert y back to log_10 for output
  
  posterior <- cbind(posterior[,1:2],posterior[,3:ncol(posterior)]/log(10))
  means <- means/log(10)
  vars <- vars/log(10)^2
  chosenSamples <- cbind(chosenSamples[,1],chosenSamples[,2:ncol(chosenSamples)]/log(10))
  posteriorModel1 <- posteriorModel1/log(10); posteriorModel2 <- posteriorModel2/log(10)
  meanfit1 <- meanfit1/log(10); meanfit2 <- meanfit2/log(10)

  # Print out results
  
  cat(rep("#",30),collapse='','\n')
  cat('Model = ', model,'\n') 
  cat('mu_max1 prior type =', mumax.prior1, ', mu_max2 prior type =', mumax.prior2, '\n')
  cat('Hypothesis = ', hyp,'(',hypinfo,')','\n')
  cat(rep("#",30),collapse='',"\n\n")
  cat('log evidence = ', logevidence, '\n')
  cat('Means and standard deviations:', '\n')
  if (hyp == "H1"){
    if (model == "linear") {
      cat('Log cell counts at time 0, y_01 = y_02 =',means[1],'+/-',(vars[1])^(1/2),'\n')
      cat('Growth rates, mu_max1 = mu_max2 =',means[2],'+/-',(vars[2])^(1/2),'\n')
    } else if (model == "logistic") {
      cat('Log cell counts at time 0, y_01 = y_02 =',means[1],'+/-',(vars[1])^(1/2),'\n')
      cat('Log final cell counts, y_max1 = y_max2 =',means[2],'+/-',(vars[2])^(1/2),'\n')
      cat('Growth rates, mu_max1 = mu_max2 =',means[3],'+/-',(vars[3])^(1/2),'\n')
    } else if (model == "Bar3par") {
      cat('Log cell counts at time 0, y_01 = y_02 =',means[1],'+/-',(vars[1])^(1/2),'\n')
      cat('Growth rates, mu_max1 = mu_max2 =',means[2],'+/-',(vars[2])^(1/2),'\n')
      cat('h_01 = h_02 =',means[3],'+/-',(vars[3])^(1/2),'\n')
    } else { # model = Bar4par
      cat('Log cell counts at time 0, y_01 = y_02',means[1],'+/-',(vars[1])^(1/2),'\n')
      cat('Log final cell counts, y_max1 = y_max2 =',means[2],'+/-',(vars[2])^(1/2),'\n')
      cat('Growth rates, mu_max1 = mu_max2 =',means[3],'+/-',(vars[3])^(1/2),'\n')
      cat('h_01 = h_02 =',means[4],'+/-',(vars[4])^(1/2),'\n')
    }
  }
  if (hyp == "H2"){
    if (model == "linear") {
      cat('Log cell counts at time 0, y_01 =',means[1],'+/-',(vars[1])^(1/2),', y_02 =',means[3],'+/-',(vars[3])^(1/2),'\n')
      cat('Growth rates, mu_max1 = mu_max2 =',means[2],'+/-',(vars[2])^(1/2),'\n')
    } else if (model == "logistic") {
      cat('Log cell counts at time 0, y_01 =',means[1],'+/-',(vars[1])^(1/2),', y_02 =',means[4],'+/-',(vars[4])^(1/2),'\n')
      cat('Log final cell counts, y_max1 =',means[2],'+/-',(vars[2])^(1/2),', y_max2 =',means[5],'+/-',(vars[5])^(1/2),'\n')
      cat('Growth rates, mu_max1 = mu_max2 =',means[3],'+/-',(vars[3])^(1/2),'\n')
    } else if (model == "Bar3par") {
      cat('Log cell counts at time 0, y_01 =',means[1],'+/-',(vars[1])^(1/2),', y_02 =',means[4],'+/-',(vars[4])^(1/2),'\n')
      cat('Growth rates, mu_max1 = mu_max2 =',means[2],'+/-',(vars[2])^(1/2),'\n')
      cat('h_01 =',means[3],'+/-',(vars[3])^(1/2),', h_02 =',means[5],'+/-',(vars[5])^(1/2),'\n')
    } else { # model = Bar4par
      cat('Log cell counts at time 0, y_01 =',means[1],'+/-',(vars[1])^(1/2),', y_02 =',means[5],'+/-',(vars[5])^(1/2),'\n')
      cat('Log final cell counts, y_max1 =',means[2],'+/-',(vars[2])^(1/2),', y_max2 =',means[6],'+/-',(vars[6])^(1/2),'\n')
      cat('Growth rates, mu_max1 = mu_max2 =',means[3],'+/-',(vars[3])^(1/2),'\n')
      cat('h_01 =',means[4],'+/-',(vars[4])^(1/2),', h_02 =',means[7],'+/-',(vars[7])^(1/2),'\n')
    }
  }
  if (hyp == "H3"){
    if (model == "linear") {
      cat('Log cell counts at time 0, y_01 =',means[1],'+/-',(vars[1])^(1/2),', y_02 =',means[3],'+/-',(vars[3])^(1/2),'\n')
      cat('Growth rates, mu_max1 =',means[2],'+/-',(vars[2])^(1/2),', mu_max2 =',means[4],'+/-',(vars[4])^(1/2),'\n')
    } else if (model == "logistic") {
      cat('Log cell counts at time 0, y_01 =',means[1],'+/-',(vars[1])^(1/2),', y_02 =',means[4],'+/-',(vars[4])^(1/2),'\n')
      cat('Log final cell counts, y_max1 =',means[2],'+/-',(vars[2])^(1/2),', y_max2 =',means[5],'+/-',(vars[5])^(1/2),'\n')
      cat('Growth rates, mu_max1 =',means[3],'+/-',(vars[3])^(1/2),', mu_max2 =',means[6],'+/-',(vars[6])^(1/2),'\n')
    } else if (model == "Bar3par") {
      cat('Log cell counts at time 0, y_01 =',means[1],'+/-',(vars[1])^(1/2),', y_02 =',means[4],'+/-',(vars[4])^(1/2),'\n')
      cat('Growth rates, mu_max1 =',means[2],'+/-',(vars[2])^(1/2),', mu_max2 =',means[5],'+/-',(vars[5])^(1/2),'\n')
      cat('h_01 =',means[3],'+/-',(vars[3])^(1/2),', h_02 =',means[6],'+/-',(vars[6])^(1/2),'\n')
    } else { # model = Bar4par
      cat('Log cell counts at time 0, y_01 =',means[1],'+/-',(vars[1])^(1/2),', y_02 =',means[5],'+/-',(vars[5])^(1/2),'\n')
      cat('Log final cell counts, y_max1 =',means[2],'+/-',(vars[2])^(1/2),', y_max2 =',means[6],'+/-',(vars[6])^(1/2),'\n')
      cat('Growth rates, mu_max1 =',means[3],'+/-',(vars[3])^(1/2),', mu_max2 =',means[7],'+/-',(vars[7])^(1/2),'\n')
      cat('h_01 =',means[4],'+/-',(vars[4])^(1/2),', h_02 =',means[8],'+/-',(vars[8])^(1/2),'\n')
    }
  }
  if (inf.sigma1 && inf.sigma2) {
    cat('Noise level, sigma1 =',means[hyp.par.no+1],'+/-',(vars[par.no+1])^(1/2),'\n')
    cat('Noise level, sigma2 =',means[hyp.par.no+2],'+/-',(vars[par.no+2])^(1/2),'\n\n')
  } else if (inf.sigma1) {
    cat('Noise level, sigma1 =',means[hyp.par.no+1],'+/-',(vars[par.no+1])^(1/2),'\n')
    cat('Noise level, sigma2 = prescribed at',sigma2.save,'\n\n')
  } else if (inf.sigma2) {
    cat('Noise level, sigma1 = prescribed at',sigma1.save,'\n\n')
    cat('Noise level, sigma2 =',means[hyp.par.no+1],'+/-',(vars[par.no+2])^(1/2),'\n\n')
  } else {
    cat('Noise level, sigma1 = prescribed at',sigma1.save,'\n')
    cat('Noise level, sigma2 = prescribed at',sigma2.save,'\n\n')
  }
  
  return(list(posterior=posterior,logevidence=logevidence,means=means,vars=vars,equalposterior=chosenSamples,
              fit.t1=tfit1,fit.y1=posteriorModel1,fit.t2=tfit2,fit.y2=posteriorModel2,fit.y1mean1=meanfit1,fit.y2mean=meanfit2))
  ### Returns: 
  ###
  ### posterior: The samples from the posterior, together with their log weights and 
  ###            log likelihoods as a m x n matrix, where m is the
  ###            number of posterior samples and n is the number of
  ###            parameters + 2. The log weights are the first column and the log likelihood values are the
  ###            second column of this matrix. The sum of the log-weights = logZ.
  ###
  ### logevidence: The logarithm of the evidence, a scalar.
  ###
  ### means: A vector of the mean of each parameter, length = no. of parameters.
  ###
  ### vars: A vector of the variance of each parameter, length = no. of parameters.
  ###
  ### equalposterior: Equally weighted posterior samples together with their 
  ###                 log likelihoods as a m x n matrix, where m is the
  ###                 number of posterior samples and n is the number of
  ###                 parameters + 1. The log likelihood values are the
  ###                 first column of this matrix.
  ###
  ### fit.t1,fit.t2: Vectors of time points at which the model is fitted for data1 and data2 respectively.
  ###
  ### fit.y1,fit.y2: Matrices of fitted model points, y1 (for data1) and y2 (for data 2), using posterior parameter samples in the  
  ###                model. Each column represents a different posterior sample.
  ###
  ### fit.y1mean and fit.y2mean: Vectors of fitted model points, y1 and y2, using the mean of the posterior parameter samples in the  
  ###                            model.
  
}, ex=function() {
  LmH_411.file <- system.file("extdata", "LmH_411.csv", package = "babar")
  LmH_411.data <- read.csv(LmH_411.file, header=TRUE, sep =",",
                           na.strings=c("ND","NA"))
  M126_50.file <- system.file("extdata", "M126_50.csv", package = "babar")
  M126_50.data <- read.csv(M126_50.file, header=TRUE, sep =",",
                           na.strings=c("ND","NA"))
  
  # Get a quick approximation of the evidence/model parameters.
  results.linear.short <- Bayescompare(LmH_411.data, M126_50.data, hyp="H1",
                                       model="linear",tol=100, prior.size=25)
  
  # Compute a better estimate of the evidence/model parameters (slow so not
  # run as part of the automated examples).
  results.linear.full <- Bayescompare(LmH_411.data, M126_50.data, hyp="H1", model="linear")


})

