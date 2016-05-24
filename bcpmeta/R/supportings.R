###################################################################
######                  Supporting functions                #######
###################################################################


###########################################################
## 1. Baic functions: marginal likelihood calculations   ##
###########################################################

#### create meta or eta from non-zero-one sequences
loc2dirac = function(loc, n){
  dirac = rep(0, n);
  dirac[loc] = 1; dirac[1] = 0;
  return(dirac);
}

#### create a (m + 1) dimensional squared matrix A, m != 0
matA = function(phi, m){
  tmp = matrix(0, nrow = m + 1, ncol = m + 1);
  indexoff = cbind(c(1:m, 2:(m + 1)), c(2:(m + 1), 1:m));
  tmp[indexoff] = -phi;
  diag(tmp) = c(1, rep(1 + phi^2, m - 1), 1);
  return(tmp);
}

# Default parameter choices
# nu0 = 5; phi.upper = 0.99; phi.lower = -0.99;
# trend = TRUE or FALSE
lMarLik.normAR = function(Xtilde, eta, nu0, phi.lower, phi.upper, trend){

  n = length(Xtilde);
  m = sum(eta);
  tau = which(eta == 1);
  
  #### integrand p(Xtilde | phi, eta) * p(phi) ####
  ## these are fixed: eta, Xtilde, nu0
  joint.Xtilde.phi = function(phi){
  	
    ## calculate DtAD
    if(m == 0)
      DtAD = matrix((1 - phi^2) + (1 - phi)^2 * (n - 1), nrow = 1, ncol = 1);
    if(m > 0){
      DtAD = matA(phi, m);
      diag(DtAD) = diag(DtAD) + (c(tau, n + 1) - c(1, tau) - 1) * (1- phi)^2;
    }    
    
    ## calculate XtdtAXtd and XtdtAT
    predXtd = c(sqrt(1 - phi^2) * Xtilde[1], Xtilde[-1] - Xtilde[-n] * phi);
    XtdtAXtd = sum(predXtd^2);
    if(trend == TRUE){
      predT = c(sqrt(1 - phi^2) * 1, (2:n) - (1:(n-1)) * phi);
      TtAT = sum(predT^2);
      XtdtAT = sum(predXtd * predT);
    }
  
    ## calculate XtdtAD and TtAD
    if(m == 0){
      XtdtAD = sum(predXtd * c(sqrt(1 - phi^2), rep(1 - phi, n - 1)));
    }
    if(m > 0){
  	  XtdtAD = rep(0, m + 1);
  	  XtdtAD[1] = sum(predXtd[1 : tau[1]] * c(sqrt(1 - phi^2), rep(1 - phi, tau[1] - 2), -phi));
  	  tmp.index = cbind(tau[-m], tau[-1]);
  	  sum1 = function(tauj){
        return(
          sum(predXtd[tauj[1] : tauj[2]] * c(1, rep(1 - phi, tauj[2] - tauj[1] - 1), -phi))
        );  
      }
  	  if(length(tmp.index) > 0){
  	    XtdtAD[2:m] = apply(tmp.index, 1, sum1);
  	  }
  	  XtdtAD[m + 1] = sum(predXtd[tau[m]:n] * c(1, rep(1 - phi, n - tau[m])));
    }
    
    if(trend == TRUE){
      if(m == 0){
        TtAD = sum(predT * c(sqrt(1 - phi^2), rep(1 - phi, n - 1)));
      }
      if(m > 0){
    	TtAD = rep(0, m + 1);
   	    TtAD[1] = sum(predT[1 : tau[1]] * c(sqrt(1 - phi^2), rep(1 - phi, tau[1] - 2), -phi));
  	    tmp.index = cbind(tau[-m], tau[-1]);
   	    sum2 = function(tauj){
          return(
            sum(predT[tauj[1] : tauj[2]] * c(1, rep(1 - phi, tauj[2] - tauj[1] - 1), -phi))
          );  
        }
  	    if(length(tmp.index) > 0){
  	      TtAD[2:m] = apply(tmp.index, 1, sum2);
  	    }
  	    TtAD[m + 1] = sum(predT[tau[m]:n] * c(1, rep(1 - phi, n - tau[m])));
      }
    }
    
    invDtADnu0 = solve(DtAD + diag(1 / nu0, m + 1));
    XtdtBXtd = XtdtAXtd - t(XtdtAD) %*% invDtADnu0 %*% XtdtAD;
    if(trend == TRUE){
      TtBT = TtAT - t(TtAD) %*% invDtADnu0 %*% TtAD;
      XtdtBT = XtdtAT - t(XtdtAD) %*% invDtADnu0 %*% TtAD;
    
      ljoint = 0.5 * (log(1 - phi^2) - log(det(DtAD + diag(1 / nu0, m + 1))) - log(TtBT)) - (n - 1) /2 * log(2 * pi) - (m + 1) / 2 * log(nu0) + lgamma((n - 1) / 2) - (n - 1) / 2 * log( (XtdtBXtd - XtdtBT^2 / TtBT) / 2 );
    }
    
    if(trend == FALSE){
      ljoint = 0.5 * (log(1 - phi^2) - log(det(DtAD + diag(1 / nu0, m + 1)))) - n /2 * log(2 * pi) - (m + 1) / 2 * log(nu0) + lgamma(n / 2) - n / 2 * log( XtdtBXtd / 2 );    	
    }
    
    return(exp(ljoint));  	
  }
  
  lml = integrate(Vectorize(joint.Xtilde.phi), lower = phi.lower, upper = phi.upper);
  
  return(list(lml = log(lml$value), lerr = log(lml$abs.error)));
}

########## Empirical Bayes (EB) estimates of sigmasq and phi ############
# Default parameter choices
# nu0 = 5; phi.upper = 0.99; phi.lower = -0.99;
# trend = TRUE or FALSE
lMarLik.EB.normAR = function(Xtilde, eta, nu0, phi.lower, phi.upper, trend){

  n = length(Xtilde);
  m = sum(eta);
  tau = which(eta == 1);
  
  #### integrand p(S | phi, gamma) * p(phi) ####
  ## these are fixed: gamma, S (which should be tilde(S) = S - mu0), nu0
  joint.Xtilde.phi = function(phi){
  	
    ## calculate DtAD
    if(m == 0)
      DtAD = matrix((1 - phi^2) + (1 - phi)^2 * (n - 1), nrow = 1, ncol = 1);
    if(m > 0){
      DtAD = matA(phi, m);
      diag(DtAD) = diag(DtAD) + (c(tau, n + 1) - c(1, tau) - 1) * (1- phi)^2;
    }    
    
    ## calculate XtdtAXtd and XtdtAT
    predXtd = c(sqrt(1 - phi^2) * Xtilde[1], Xtilde[-1] - Xtilde[-n] * phi);
    XtdtAXtd = sum(predXtd^2);
    if(trend == TRUE){
      predT = c(sqrt(1 - phi^2) * 1, (2:n) - (1:(n-1)) * phi);
      TtAT = sum(predT^2);
      XtdtAT = sum(predXtd * predT);
    }
  
    ## calculate XtdtAD and TtAD
    if(m == 0){
      XtdtAD = sum(predXtd * c(sqrt(1 - phi^2), rep(1 - phi, n - 1)));
    }
    if(m > 0){
  	  XtdtAD = rep(0, m + 1);
  	  XtdtAD[1] = sum(predXtd[1 : tau[1]] * c(sqrt(1 - phi^2), rep(1 - phi, tau[1] - 2), -phi));
  	  tmp.index = cbind(tau[-m], tau[-1]);
  	  sum1 = function(tauj){
        return(
          sum(predXtd[tauj[1] : tauj[2]] * c(1, rep(1 - phi, tauj[2] - tauj[1] - 1), -phi))
        );  
      }
  	  if(length(tmp.index) > 0){
  	    XtdtAD[2:m] = apply(tmp.index, 1, sum1);
  	  }
  	  XtdtAD[m + 1] = sum(predXtd[tau[m]:n] * c(1, rep(1 - phi, n - tau[m])));
    }
    
    if(trend == TRUE){
      if(m == 0){
        TtAD = sum(predT * c(sqrt(1 - phi^2), rep(1 - phi, n - 1)));
      }
      if(m > 0){
    	TtAD = rep(0, m + 1);
   	    TtAD[1] = sum(predT[1 : tau[1]] * c(sqrt(1 - phi^2), rep(1 - phi, tau[1] - 2), -phi));
  	    tmp.index = cbind(tau[-m], tau[-1]);
   	    sum2 = function(tauj){
          return(
            sum(predT[tauj[1] : tauj[2]] * c(1, rep(1 - phi, tauj[2] - tauj[1] - 1), -phi))
          );  
        }
  	    if(length(tmp.index) > 0){
  	      TtAD[2:m] = apply(tmp.index, 1, sum2);
  	    }
  	    TtAD[m + 1] = sum(predT[tau[m]:n] * c(1, rep(1 - phi, n - tau[m])));
      }
    }
    
    invDtADnu0 = solve(DtAD + diag(1 / nu0, m + 1));
    XtdtBXtd = XtdtAXtd - t(XtdtAD) %*% invDtADnu0 %*% XtdtAD;
    if(trend == TRUE){
      TtBT = TtAT - t(TtAD) %*% invDtADnu0 %*% TtAD;
      XtdtBT = XtdtAT - t(XtdtAD) %*% invDtADnu0 %*% TtAD;
    
      ljoint = 0.5 * (log(1 - phi^2) - log(det(DtAD + diag(1 / nu0, m + 1))) - log(TtBT)) - (n - 1) /2 * log(2 * pi / (n - 1) * (XtdtBXtd - XtdtBT^2 / TtBT)) - (m + 1) / 2 * log(nu0) - (n - 1) / 2;
    }
    
    if(trend == FALSE){
      ljoint = 0.5 * (log(1 - phi^2) - log(det(DtAD + diag(1 / nu0, m + 1)))) - n /2 * log(2 * pi / n * XtdtBXtd) - (m + 1) / 2 * log(nu0) - n / 2;    	
    }
    
    return(ljoint);  	
  }
  
  eb = optimize(joint.Xtilde.phi, interval = c(phi.lower, phi.upper), maximum = TRUE);
  return(list(lml = c(eb$objective), phi.eb = c(eb$maximum)));
}

###################################################################
## 2. MCMC single iteration: Metropolis-Hastings update of eta   ##
###################################################################

##########  Update one random component of eta at a time  ##########
## current = list(eta, lml, lpost, lerr, change.eta);
## if EB == TRUE, lerr is replaced by phi = phi.eb

up.eta.MH.normAR = function(Xtilde, meta, current, a1, a2, b1, b2, nu0, phi.lower, phi.upper, trend, EB){
  
  eta = current$eta;

  ## randomly select the component to be flipped 
  n = length(Xtilde);
  j = sample(2:n, 1);
  # j = 51
  eta2 = eta;
  eta2[j] = 1 - eta[j];
  
  m2.out = sum(eta2[(which(meta == 0))] == 1);
  m2.in = sum(eta2[(which(meta == 1))] == 1);

  if(EB == FALSE)
    lmarlik2 = lMarLik.normAR(Xtilde, eta2, nu0, phi.lower, phi.upper, trend); 
  if(EB == TRUE)
    lmarlik2 = lMarLik.EB.normAR(Xtilde, eta2, nu0, phi.lower, phi.upper, trend); 
   
  lpost2 = lmarlik2$lml + lbeta(a1 + m2.out, b1 + n - 1 - sum(meta) - m2.out) + lbeta(a2 + m2.in, b2 + sum(meta) - m2.in) - lbeta(a1, b1) - lbeta(a2, b2);

  loga = lpost2 - current$lpost;
  logu = log(runif(1));
  if(loga > logu){
  	current$eta = eta2;
  	current$lml = lmarlik2$lml;
  	current$lpost = lpost2;
  	current$change.eta = TRUE;
  	if(EB == FALSE)
  	  current$lerr = lmarlik2$lerr;
  	if(EB == TRUE)
  	  current$phi = lmarlik2$phi.eb;  	
  }
  if(loga <= logu)
    current$change.eta = FALSE;
  
  return(current);
}

############  Simple random swapping update of eta  ###########
## current = list(eta, lml, lpost, lerr, change.eta);
## if EB == TRUE, lerr is replaced by phi = phi.eb

up.eta.RS.normAR = function(Xtilde, meta, current, a1, a2, b1, b2, nu0, phi.lower, phi.upper, trend, EB){
  
  eta = current$eta;
  n = length(Xtilde);
  
  if(sum(eta) > 0 && sum(eta) < n - 1){
  	if(sum(eta) == 1)
   	  i = which(eta == 1);
   	if(sum(eta) > 1)
   	  i = sample(which(eta == 1), 1);
   	
   	if(sum(eta) == n - 2)
   	  j = which(eta == 0)[2];
   	if(sum(eta) < n - 2)
      j = sample(which(eta == 0)[-1], 1);

    ## randomly select the component to be flipped 
    eta2 = eta;
    eta2[i] = 1 - eta[i];
    eta2[j] = 1 - eta[j];
  
    m2.out = sum(eta2[(which(meta == 0))] == 1);
    m2.in = sum(eta2[(which(meta == 1))] == 1);
    
    if(EB == FALSE)
      lmarlik2 = lMarLik.normAR(Xtilde, eta2, nu0, phi.lower, phi.upper, trend); 
    if(EB == TRUE)
      lmarlik2 = lMarLik.EB.normAR(Xtilde, eta2, nu0, phi.lower, phi.upper, trend); 
      
    lpost2 = lmarlik2$lml + lbeta(a1 + m2.out, b1 + n - 1 - sum(meta) - m2.out) + lbeta(a2 + m2.in, b2 + sum(meta) - m2.in) - lbeta(a1, b1) - lbeta(a2, b2);

    loga = lpost2 - current$lpost;
    logu = log(runif(1));
    if(loga > logu){
  	  current$eta = eta2;
  	  current$lml = lmarlik2$lml;
  	  current$lpost = lpost2;
  	  current$change.eta = TRUE;
  	  if(EB == FALSE)
  	    current$lerr = lmarlik2$lerr;
  	  if(EB == TRUE)
  	    current$phi = lmarlik2$phi.eb;
    }
    if(loga <= logu)
      current$change.eta = FALSE;
  }
  else{
  	current$change.eta = FALSE;
  }
  
  return(current);
}

###########################################################################
## 3. MCMC: model parameters Mu, Alpha, Sigmasq and Phi (if EB == FALSE) ##
###########################################################################

#############################
############ EB #############
#############################
## find EB estimate for sigamsq, and
## update alpha (if trend == TRUE), and mu (length N) using Gibbs sampler (no need to thin!)
## current = list(phi, sigmasq, alpha, mu);
up.parameters.EB.normAR = function(Xtilde, eta, phi.eb, iter = 1e3, nu0, trend){
  
  n = length(Xtilde);
  m = sum(eta);
  tau = which(eta == 1);
  phi = phi.eb;
  
  if(m == 0)
    DtAD = matrix((1 - phi^2) + (1 - phi)^2 * (n - 1), nrow = 1, ncol = 1);
  if(m > 0){
    DtAD = matA(phi, m);
    diag(DtAD) = diag(DtAD) + (c(tau, n + 1) - c(1, tau) - 1) * (1- phi)^2;
  }    
    
  ## calculate XtdtAXtd and XtdtAT 
  predXtd = c(sqrt(1 - phi^2) * Xtilde[1], Xtilde[-1] - Xtilde[-n] * phi);
  XtdtAXtd = sum(predXtd^2);
  if(trend == TRUE){
    predT = c(sqrt(1 - phi^2) * 1, (2:n) - (1:(n-1)) * phi);
    TtAT = sum(predT^2);
    XtdtAT = sum(predXtd * predT);
  }
  
  ## calculate XtdtAD and TtAD
  if(m == 0){
    XtdtAD = sum(predXtd * c(sqrt(1 - phi^2), rep(1 - phi, n - 1)));
  }
  if(m > 0){
  	XtdtAD = rep(0, m + 1);
  	XtdtAD[1] = sum(predXtd[1 : tau[1]] * c(sqrt(1 - phi^2), rep(1 - phi, tau[1] - 2), -phi));
    tmp.index = cbind(tau[-m], tau[-1]);
  	sum1 = function(tauj){
      return(
        sum(predXtd[tauj[1] : tauj[2]] * c(1, rep(1 - phi, tauj[2] - tauj[1] - 1), -phi))
      );  
    }
  	if(length(tmp.index) > 0){
  	  XtdtAD[2:m] = apply(tmp.index, 1, sum1);
  	}
    XtdtAD[m + 1] = sum(predXtd[tau[m]:n] * c(1, rep(1 - phi, n - tau[m])));
  }
    
  if(trend == TRUE){
    if(m == 0){
      TtAD = sum(predT * c(sqrt(1 - phi^2), rep(1 - phi, n - 1)));
    }
    if(m > 0){
      TtAD = rep(0, m + 1);
   	  TtAD[1] = sum(predT[1 : tau[1]] * c(sqrt(1 - phi^2), rep(1 - phi, tau[1] - 2), -phi));
  	  tmp.index = cbind(tau[-m], tau[-1]);
   	  sum2 = function(tauj){
        return(
          sum(predT[tauj[1] : tauj[2]] * c(1, rep(1 - phi, tauj[2] - tauj[1] - 1), -phi))
        );  
      }
	  if(length(tmp.index) > 0){
  	    TtAD[2:m] = apply(tmp.index, 1, sum2);
      }
	  TtAD[m + 1] = sum(predT[tau[m]:n] * c(1, rep(1 - phi, n - tau[m])));
    }
  }
    
  invDtADnu0 = solve(DtAD + diag(1 / nu0, m + 1));
  XtdtBXtd = XtdtAXtd - t(XtdtAD) %*% invDtADnu0 %*% XtdtAD;  
  
  if(trend == TRUE){
    TtBT = TtAT - t(TtAD) %*% invDtADnu0 %*% TtAD;
    XtdtBT = XtdtAT - t(XtdtAD) %*% invDtADnu0 %*% TtAD;
  	sigmasq = c((XtdtBXtd - XtdtBT^2 / TtBT) / (n - 1));
  	Alpha = rnorm(iter, XtdtBT / TtBT, sqrt(sigmasq / TtBT));
  	Mu = matrix(NA, nrow = iter, ncol = m + 1);
  	for(i in 1:iter)
  	  Mu[i, ] = rmvnorm(1, invDtADnu0 %*% (XtdtAD - Alpha[i] * TtAD), sigmasq * invDtADnu0);
  	#Mu = Mu[, cumsum(eta) + 1];
  	#colnames(Mu) = paste('T', 1:n, sep = '');
  	return(list(sigmasq = sigmasq, Alpha = Alpha, Mu = Mu));
  }
  if(trend == FALSE){
  	sigmasq = c(XtdtBXtd / n);
  	Mu = rmvnorm(iter, invDtADnu0 %*% XtdtAD, sigmasq * invDtADnu0);
  	#Mu = Mu[, cumsum(eta) + 1];
  	#colnames(Mu) = paste('T', 1:n, sep = '');
    return(list(sigmasq = sigmasq, Mu = Mu));
  }
  
}

#############################
############ FB #############
#############################

## Calculate conditional marginal likelihood given phi and eta  
lMarLik.phi.normAR = function(Xtilde, phi, eta, nu0, trend){
  	
  n = length(Xtilde);
  m = sum(eta);
  tau = which(eta == 1);

  ## calculate XtAX
  if(m == 0)
    DtAD = matrix((1 - phi^2) + (1 - phi)^2 * (n - 1), nrow = 1, ncol = 1);
  if(m > 0){
    DtAD = matA(phi, m);
    diag(DtAD) = diag(DtAD) + (c(tau, n + 1) - c(1, tau) - 1) * (1- phi)^2;
  }    
    
  ## calculate XtdtAXtd and XtdtAT 
  predXtd = c(sqrt(1 - phi^2) * Xtilde[1], Xtilde[-1] - Xtilde[-n] * phi);
  XtdtAXtd = sum(predXtd^2);
  if(trend == TRUE){
    predT = c(sqrt(1 - phi^2) * 1, (2:n) - (1:(n-1)) * phi);
    TtAT = sum(predT^2);
    XtdtAT = sum(predXtd * predT);
  }
  
  ## calculate XtdtAD and TtAD
  if(m == 0){
    XtdtAD = sum(predXtd * c(sqrt(1 - phi^2), rep(1 - phi, n - 1)));
  }
  if(m > 0){
  	XtdtAD = rep(0, m + 1);
  	XtdtAD[1] = sum(predXtd[1 : tau[1]] * c(sqrt(1 - phi^2), rep(1 - phi, tau[1] - 2), -phi));
  	tmp.index = cbind(tau[-m], tau[-1]);
  	sum1 = function(tauj){
      return(
        sum(predXtd[tauj[1] : tauj[2]] * c(1, rep(1 - phi, tauj[2] - tauj[1] - 1), -phi))
      );  
    }
  	if(length(tmp.index) > 0){
  	  XtdtAD[2:m] = apply(tmp.index, 1, sum1);
  	}
  	XtdtAD[m + 1] = sum(predXtd[tau[m]:n] * c(1, rep(1 - phi, n - tau[m])));
  }
    
  if(trend == TRUE){
    if(m == 0){
      TtAD = sum(predT * c(sqrt(1 - phi^2), rep(1 - phi, n - 1)));
    }
    if(m > 0){
      TtAD = rep(0, m + 1);
   	  TtAD[1] = sum(predT[1 : tau[1]] * c(sqrt(1 - phi^2), rep(1 - phi, tau[1] - 2), -phi));
  	  tmp.index = cbind(tau[-m], tau[-1]);
   	  sum2 = function(tauj){
        return(
          sum(predT[tauj[1] : tauj[2]] * c(1, rep(1 - phi, tauj[2] - tauj[1] - 1), -phi))
        );  
      }
  	  if(length(tmp.index) > 0){
  	    TtAD[2:m] = apply(tmp.index, 1, sum2);
      }
  	  TtAD[m + 1] = sum(predT[tau[m]:n] * c(1, rep(1 - phi, n - tau[m])));
    }
  }
    
  invDtADnu0 = solve(DtAD + diag(1 / nu0, m + 1));
  XtdtBXtd = XtdtAXtd - t(XtdtAD) %*% invDtADnu0 %*% XtdtAD;
  
  if(trend == TRUE){
    TtBT = TtAT - t(TtAD) %*% invDtADnu0 %*% TtAD;
    XtdtBT = XtdtAT - t(XtdtAD) %*% invDtADnu0 %*% TtAD;
    
    lml.phi = 0.5 * (log(1 - phi^2) - log(det(DtAD + diag(1 / nu0, m + 1))) - log(TtBT)) - (n - 1) /2 * log(2 * pi) - (m + 1) / 2 * log(nu0) + lgamma((n - 1) / 2) - (n - 1) / 2 * log( (XtdtBXtd - XtdtBT^2 / TtBT) / 2 );
    
    return(list(lml.phi = lml.phi, invDtADnu0 = invDtADnu0, XtdtAD = XtdtAD, TtAD = TtAD, XtdtBT = XtdtBT, TtBT = TtBT, XtdtBXtd = XtdtBXtd));  	
    
  }
    
  if(trend == FALSE){
    lml.phi = 0.5 * (log(1 - phi^2) - log(det(DtAD + diag(1 / nu0, m + 1)))) - n /2 * log(2 * pi) - (m + 1) / 2 * log(nu0) + lgamma(n / 2) - n / 2 * log( XtdtBXtd / 2 );    	
    return(list(lml.phi = lml.phi, invDtADnu0 = invDtADnu0, XtdtAD = XtdtAD, XtdtBXtd = XtdtBXtd));  	
  }
    
}
  

##############################################
## update phi using Metropolis-Hastings algorithm
## update sigmasq, alpha (if trend == TRUE), and mu (length N) using Gibbs sampler
## current = list(lml.phi, phi, sigmasq, alpha, mu, change.phi);

up.parameters.normAR = function(Xtilde, meta, eta, current, nu0, phi.lower, phi.upper, trend, sd.xi){
	
  n = length(Xtilde);	
  m = sum(eta);
  tau = which(eta == 1);
  phi = current$phi;
  lml = current$lml.phi;
  
  ##### M-H update of phi
  xi = log((phi - phi.lower) / (phi.upper - phi));
  xi2 = rnorm(1, phi, sd.xi);
  phi2 = (exp(xi2) * phi.upper + phi.lower) / (exp(xi2) + 1);
  lml2 = lMarLik.phi.normAR(Xtilde, phi2, eta, nu0, trend);
  
  loga = lml2$lml.phi - lml$lml.phi + log(phi2 - phi.lower) + log(phi.upper - phi2) - log(phi - phi.lower) - log(phi.upper - phi);
  logu = log(runif(1));
  if(loga > logu){
    current$phi = phi = phi2;
    current$lml.phi = lml = lml2;    
    current$change.phi = current$change.phi + 1;
  }
  
  #### Gibbs update of sigma, alpha, mu (if trend == TRUE)
  if(trend == TRUE){
  	current$sigmasq = 1 / rgamma(1, 0.5 * (n - 1), 0.5 * (lml$XtdtBXtd - lml$XtdtBT^2 / lml$TtBT) );
  	current$alpha = rnorm(1, lml$XtdtBT / lml$TtBT, sqrt(current$sigmasq / lml$TtBT));
  	mu = rmvnorm(1, lml$invDtADnu0 %*% (lml$XtdtAD - current$alpha * lml$TtAD), current$sigmasq * lml$invDtADnu0);
  	current$mu = mu;
  	#current$mu = mu[cumsum(eta) + 1];
  }
  
  #### Gibbs update of sigma, alpha, mu (if trend == FALSE)
  if(trend == FALSE){
  	current$sigmasq = 1 / rgamma(1, 0.5 * n, 0.5 * lml$XtdtBXtd);
  	mu = rmvnorm(1, lml$invDtADnu0 %*% lml$XtdtAD, current$sigmasq * lml$invDtADnu0);
  	current$mu = mu;
  	#current$mu = mu[cumsum(eta) + 1];
  }
  
  return( current );
}


