regression_warp <- function(beta, time, q, y, alpha){
  gam_M = optimum.reparam(beta,time,q,time)
  qM = warp_q_gamma(time, q, gam_M)
  y_M = trapz(time, qM * beta)

  gam_m = optimum.reparam(-1 * beta,time,q,time)
  qm = warp_q_gamma(time, q, gam_m)
  y_m = trapz(time, qm * beta)


  if (y > alpha + y_M){
    gamma_new = gam_M
  }else if (y < alpha + y_m){
    gamma_new - gam_m
  }else{
    gamma_new = zero_crossing(y-alpha, q, beta, time, y_M, y_m, gam_M, gam_m)
  }

  return(gamma_new)
}

logistic_warp <- function(beta, time, q, y){
  if (y==1){
    gamma = optimum.reparam(beta,time,q,time)
  } else if (y== -1){
    gamma = optimum.reparam(-1*beta,time,q,time)
  }

  return(gamma)
}

phi <- function(t){
  idx = t > 0
  out = rep(0, length(t))
  out[idx] = 1/(1+exp(-1*t[idx]))
  exp_t = exp(t[!idx])
  out[!idx] = exp_t / (1+ exp_t)

  return(out)
}

logit_loss <- function(b, X, y){
  z = X %*% b
  yz = y * z
  idx = yz > 0
  out = rep(0, length(yz))
  out[idx] = log(1+ exp(-1*yz[idx]))
  out[!idx] = (-1*yz[!idx] + log(1+exp(yz[!idx])))
  out = sum(out)
  return(out)
}

logit_gradient <- function(b, X, y){
  z = X %*% b
  z = phi(y * z)
  z0 = (z-1) * y
  grad = t(X) %*% z0
  return(grad)
}

logit_hessian <- function(s, b, X, y){
  z = X %*% b
  z = phi(y * z)
  d = z * (1 - z)
  wa = d * (X %*% s)
  Hs = t(X) %*% wa
  return(Hs)
}

mlogit_warp <- function(alpha, beta, time, q, y, max_itr=8000, tol=1e-10,
                        delta=0.008, display=0){
  m1 = length(time)
  m2 = ncol(beta)
  gam1 = seq(0,1,length.out=m1)
  gamout = rep(0,m1)
  alpha = alpha/pvecnorm(alpha,2)
  q = q/pvecnorm(q,2)
  for (ii in 1:m2){
    beta[,ii] = beta[,ii]/pvecnorm(beta[,ii],2)
  }
  beta1 = rep(0,m1*m2)
  for (ii in 1:m2){
    beta1[((ii-1)*m1+1):(ii*m1)] = beta[,ii]
  }
  output = .Call('mlogit_warp_grad_wrap', PACKAGE = 'fdasrvf', m1, m2, alpha,
                 beta1, time, gam1, q, y, max_itr, tol, delta, display, gamout);

  out = output$gamout
  return(out)
}

mlogit_loss <- function(b, X, Y){
  N = nrow(Y)
  m = ncol(Y)
  M = ncol(X)
  B = array(b,c(M, m))
  Yhat = X %*% B
  Yhat = Yhat - apply(Yhat,1,min)
  Yhat = exp(-1*Yhat)
  # l1-normalize
  Yhat = Yhat/apply(Yhat,1,sum)

  Yhat = Yhat * Y
  nll = sum(log(apply(Yhat,1,sum)))
  nll = -1*nll / N
  return(nll)
}

mlogit_gradient <- function(b, X, Y){
  N = nrow(Y)
  m = ncol(Y)
  M = ncol(X)
  B = array(b,c(M, m))
  Yhat = X %*% B
  Yhat = Yhat - apply(Yhat,1,min)
  Yhat = exp(-1*Yhat)
  # l1-normalize
  Yhat = Yhat/apply(Yhat,1,sum)

  Yhat1 = Yhat * Y
  Yhat1 = Yhat1 / apply(Yhat1,1,sum)
  Yhat = Yhat - Yhat1
  grad = t(X) %*% Yhat
  grad = -1*grad/N
  return(grad)
}
