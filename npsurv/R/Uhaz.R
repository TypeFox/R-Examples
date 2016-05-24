############################################
# Estimation of a U-shaped Hazard Function #
############################################

Uhaz = function(data, w=1, deg=1, maxit=100, tol=1e-6, verb=0) {
  x = icendata(data, w)
  h = uh.initial(x, deg)
  attr(h, "ll") = logLikuh(h, x)
  expdH = NULL
  bc = TRUE              # boundary change
  convergence = 1
  for(i in 1:maxit){
    h.old = h
    if(nrow(x$o) > 0) expdH = exp(chazuh(x$o[,1],h) - chazuh(x$o[,2],h))
    maxima = maxgrad(h, x, expdH, bc=bc)
    np1 = maxima$np1
    np2 = maxima$np2
    h = uh(h$alpha, c(h$tau, np1), c(h$nu, double(length(np1))), 
        c(h$eta, np2), c(h$mu, double(length(np2))),
        h$upper, h$deg, collapse=TRUE)
    r = updatemass(h, x, expdH, tol=tol)
    h = r$h
    if(h$deg == 0) {h = simplify(h); attr(h, "ll") = logLikuh(h, x)}
    if(verb>0) { cat("##### Iteration", i, "#####\n")
                 cat("Log-likelihood: ", signif(attr(h,"ll"), 6), "\n")
                 if(verb>1) cat("Gradient values: ", signif(dlogLik(h, x), 6), "\n")
                 if(verb>2) {cat("hazard function:\n"); print(h)} }
    if(r$convergence == 1) bc = FALSE    # backtracking failed.
    else if(attr(h, "ll") <= attr(h.old, "ll") + tol) {convergence = 0; break}
  }
  r = list(convergence=convergence, grad=dlogLik(h, x), numiter=i,
       ll=attr(h, "ll"), h=h)
  class(r) = "Uhaz"
  r
}

# Update masses

updatemass = function(h, x, expdH=NULL, tol=1e-10) {
  tau = h$tau
  k = length(tau)
  j2 = h$eta != h$upper
  eta = h$eta = h$eta[j2]
  h$mu = h$mu[j2]
  m = length(eta)
  p = h$deg
  D1 = D2 = NULL
  t1 = x$t[x$i1]
  n1 = length(t1)
  if(n1 > 0) {
    tau.r = rep(tau, rep.int(n1,k))
    dim(tau.r) = c(n1, k)
    if(p > 0) tau.t = pmax(tau.r - t1, 0)
    if(m > 0) {
      eta.r = rep(eta, rep.int(n1,m))
      dim(eta.r) = c(n1, m)
      if(p > 0) t.eta = pmax(t1 - eta.r, 0)
    }
    D1 = switch(as.character(p),
        "0" = cbind(1, tau.r >= t1, if(m>0) t1 >= eta.r else NULL) / hazuh(t1, h),
        "1" = cbind(1, tau.t, if(m>0) t.eta else NULL) / hazuh(t1, h),
        "2" = cbind(1, tau.t * tau.t,
            if(m>0) t.eta * t.eta else NULL) / hazuh(t1, h),
        cbind(1, tau.t^p, if(m>0) t.eta^p else NULL) / hazuh(t1, h)   )
  }
  n2 = nrow(x$o)
  if(n2 > 0) {
    if(is.null(expdH)) expdH = exp(chazuh(x$o[,1],h) - chazuh(x$o[,2],h))
    delta = sqrt(expdH) / (1 - expdH)
    tau.r1 = rep(tau, rep.int(n2,k))
    dim(tau.r1) = c(n2, k)
    tau.x1 = pmax(tau.r1 - x$o[,1], 0)
    tau.x2 = pmax(tau.r1 - x$o[,2], 0)
    xd0 = x$o[,1] - x$o[,2]
    xd1 = switch(as.character(p),
        "0" = tau.x2 - tau.x1,
        "1" = .5 * (tau.x2 * tau.x2 - tau.x1 * tau.x1),
        "2" = (tau.x2 * tau.x2 * tau.x2 - tau.x1 * tau.x1 * tau.x1) / 3,
        (tau.x2^(p+1) - tau.x1^(p+1)) / (p+1)  )
    if(m > 0) {
      eta.r2 = rep(eta, rep.int(n2,m))
      dim(eta.r2) = c(n2, m)
      x1.eta = pmax(x$o[,1] - eta.r2, 0)
      x2.eta = pmax(x$o[,2] - eta.r2, 0)
      xd2 = switch(as.character(p),
          "0" = x1.eta - x2.eta,
          "1" = .5 * (x1.eta * x1.eta - x2.eta * x2.eta),
          "2" = (x1.eta * x1.eta * x1.eta - x2.eta * x2.eta * x2.eta) / 3,
          (x1.eta^(p+1) - x2.eta^(p+1)) / (p+1) )
    }
    else xd2 = NULL
    D2 = cbind(xd0, xd1, xd2) * delta
    D2[delta == 0] = 0
  }
  D = rbind(D1, D2) * sqrt(c(x$wt[x$i1], x$wo))

  H = crossprod(D)                  # Choleski decomposition
  v = sqrt(diag(H))
  jv = v != 0
  Hv = H[jv,jv,drop=FALSE] / tcrossprod(v[jv])
  diag(Hv) = diag(Hv) + 1e-10
  Rv = chol(Hv)
  gv = dlogLik(h, x, expdH, interior=TRUE)[jv] / v[jv]
  plus = forwardsolve(Rv, gv, upper.tri=TRUE, transpose=TRUE)
  par = c(h$alpha, h$nu, h$mu[j2])[jv] * v[jv]
  w.new = double(length(v))
  w.new[jv] = nnls(Rv, Rv %*% par + plus)$x / v[jv]

  alpha = w.new[1]
  nu = if(k > 0) w.new[2:(k+1)] else numeric()
  mu = if(m > 0) w.new[(k+2):(k+m+1)] else numeric()
  newh = uh(alpha=alpha, tau=h$tau, nu=nu, eta=h$eta, mu=mu,
      upper=x$upper, h$deg, collapse=FALSE)
  if(h$deg == 0) b = backtrack(h, newh, x, expdH, alpha=0)
  else b = backtrack(h, newh, x, expdH) 
  newh = b$h2
  j1 = newh$nu != 0
  j2 = newh$mu != 0
  h2 = uh(newh$alpha, newh$tau[j1], newh$nu[j1], newh$eta[j2], newh$mu[j2],
    upper=h$upper, h$deg, collapse=FALSE)
  h = collapse(h2, x, tol=pmax(tol,1e-10))
  list(h=h, convergence=b$convergence)
}

# Backtracking line search. h and h2 must have the same knots

backtrack = function(h, h2, x, expdH, tol=1e-10, alpha=0.33){
  j = h$eta != h$upper
  h2$eta = h$eta = h$eta[j]
  h$mu = h$mu[j]
  h2$mu = h2$mu[j]
  ll.h = logLikuh(h, x)
  d = c(h2$alpha - h$alpha, h2$nu - h$nu, h2$mu - h$mu)
  g = alpha * sum(dlogLik(h, x, expdH) * d)
  convergence = 0
  r = 1
  repeat {
    hr = uh((1-r) * h$alpha + r * h2$alpha,
      h$tau, (1-r) * h$nu + r * h2$nu, 
      h$eta, (1-r)*h$mu + r * h2$mu,  upper=h$upper, deg=h$deg)
    ll.hr = logLikuh(hr, x)
    if (ll.hr >= ll.h + r * g) {convergence =0; break}
    r = 0.5 * r
    if(r < tol) {r = 0; hr = h; ll.h2 = ll.h; convergence = 1; break}
  }
  attr(h2, "ll") = ll.hr
  list(h2=hr, r=r, convergence=convergence)
}

# Zeroth gradient function

grad0 = function(h, x, expdH=NULL) {
  if(length(x$t) > 0) 
    d0 = sum(x$wt[x$i1] / hazuh(x$t[x$i1], h)) - sum(x$wt * x$t)
  else d0 = 0
  if(nrow(x$o) > 0) {
    if(is.null(expdH)) expdH = exp(chazuh(x$o[,1],h) - chazuh(x$o[,2],h))
    Delta = expdH / (1 - expdH)
    xd = x$o[,1] - x$o[,2]
    xd[Delta == 0] = 0
    d0 = d0 - sum(x$wo * (x$o[,1] + xd * Delta))
  }
  d0
}

# First gradient function

grad1 = function(tau, h, x, expdH=NULL, order=0) {
  g = vector("list", length(order))
  names(g) = paste("d", order, sep="")
  g[1:length(g)] = 0
  if(length(tau) == 0) return(NULL)
  m = length(tau)
  n1 = length(x$t)
  p = h$deg
  if(n1 > 0) {              # for exact observations
    tau.r1 = rep(tau, rep.int(n1,m))
    dim(tau.r1) = c(n1,m)
    ind = tau.r1 >= x$t
    tau.t = pmax(tau.r1 - x$t, 0)
    if(0 %in% order) {
      g$d0 = switch(as.character(p),
          "0" = (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
              ind[x$i1,,drop=FALSE]) - crossprod(x$wt, tau.r1 - tau.t))[1,],
          "1" = (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
              tau.t[x$i1,,drop=FALSE]) - 
              0.5 * crossprod(x$wt, tau.r1 * tau.r1 - tau.t * tau.t))[1,],
          "2" = (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
              tau.t[x$i1,,drop=FALSE]^2) -
              crossprod(x$wt, tau.r1*tau.r1*tau.r1 - tau.t*tau.t*tau.t) / 3)[1,],
          (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
                    tau.t[x$i1,,drop=FALSE]^p) - 
              crossprod(x$wt, tau.r1^(p+1) - tau.t^(p+1))[1,] / (p+1))[1,]  )
    }
    if(1 %in% order) {
      g$d1 = switch(as.character(p),
          "0" = double(m),
          "1" = (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
              ind[x$i1,,drop=FALSE]) - crossprod(x$wt, tau.r1 - tau.t))[1,],
          "2" = (2 * crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
              tau.t[x$i1,,drop=FALSE]) -
              crossprod(x$wt, tau.r1 * tau.r1 - tau.t * tau.t))[1,],
          { tau.t.pm1 = tau.t[x$i1,,drop=FALSE]^(p-1)
            if(p < 1) tau.t.pm1[tau.t.pm1 == Inf] = 0
            (p * crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h), tau.t.pm1) - 
              crossprod(x$wt, tau.r1^p - tau.t^p))[1,] 
          }  )
    }
    if(2 %in% order) {
      g$d2 = switch(as.character(p),
          "0" = double(m),
          "1" = - sum(x$wt) + crossprod(x$wt, ind)[1,],
          "2" = (2 * crossprod(x$wt[x$i1], ind[x$i1,,drop=FALSE]) -
                 crossprod(x$wt, tau.r1 - tau.t))[1,],
          { tau.t1.pm2 = tau.t[x$i1,,drop=FALSE]^(p-2)
            if(p < 2) tau.t1.pm2[tau.t1.pm2 == Inf] = 0
            tau.t.pm1 = tau.t^(p-1)
            if(p < 1) tau.t.pm1[tau.t.pm1 == Inf] = 0
            (p * (p-1) * crossprod(x$wt[x$i1], tau.t1.pm2) -
             p * crossprod(x$wt, tau.r1^(p-1) - tau.t.pm1))[1,]
          }  )
    }
  }
  n2 = nrow(x$o)
  if(n2 > 0) {              # for interval-censored observations
    if(is.null(expdH)) expdH = exp(chazuh(x$o[,1],h) - chazuh(x$o[,2],h))
    Delta = expdH / (1 - expdH)
    tau.r2 = rep(tau, rep.int(n2,m))
    dim(tau.r2) = c(n2,m)
    tau.x1 = pmax(tau.r2 - x$o[,1], 0)
    tau.x2 = pmax(tau.r2 - x$o[,2], 0)
    ind1 = tau.r2 >= x$o[,1]
    ind2 = tau.r2 >= x$o[,2]
    if(0 %in% order) {
      xd0 = switch(as.character(p),
          "0" = {tau.r2.p1 = tau.r2; tau.x2 - (tau.x1.p1 <- tau.x1)},
          "1" = {tau.r2.p1 = tau.r2  * tau.r2
                 tau.x2 * tau.x2 - (tau.x1.p1 <- tau.x1 * tau.x1)},
          "2" = {tau.r2.p1 = tau.r2 * tau.r2 * tau.r2
                 tau.x2 * tau.x2 * tau.x2 -
                     (tau.x1.p1 <- tau.x1 * tau.x1 * tau.x1)},
          { tau.r2.p1 = tau.r2^(p+1)
            tau.x2^(p+1) - (tau.x1.p1 <- tau.x1^(p+1))}  )
      xd0[Delta == 0] = 0
      g$d0 = g$d0 -
        crossprod(x$wo, (tau.r2.p1 - tau.x1.p1 + xd0 * Delta))[1,] / (p+1)
    }
    if(1 %in% order) {
      xd1 = switch(as.character(p),
          "0" = {tau.r2.p = 1; ind2 - (tau.x1.p <- ind1)},
          "1" = {tau.r2.p = tau.r2; tau.x2 - (tau.x1.p <- tau.x1)},
          "2" = {tau.r2.p = tau.r2 * tau.r2;
                 tau.x2 * tau.x2 - (tau.x1.p <- tau.x1 * tau.x1)},
          {tau.r2.p = tau.r2^p; tau.x2^p - (tau.x1.p <- tau.x1^p)}  )
      xd1[Delta == 0] = 0
      g$d1 = g$d1 - crossprod(x$wo, tau.r2.p - tau.x1.p + xd1 * Delta)[1,]
    }
    if(2 %in% order) {
      g$d2 = switch(as.character(p),
          "0" = double(m),
          "1" = g$d2 - sum(x$wo) +
              crossprod(x$wo, ind1 - (ind2 - ind1) * Delta)[1,],
          "2" = { xd2 = tau.x2 - tau.x1
                  xd2[Delta == 0] = 0
                  g$d2 -
                    2 * crossprod(x$wo, tau.r2 - tau.x1 + xd2 * Delta)[1,] },
          { tau.x1.pm1 = tau.x1^(p-1)
            if(p < 1) tau.x1.pm1[tau.x1.pm1 == Inf] = 0
            xdp = tau.x2^(p-1) - tau.x1^(p-1)
            xdp[Delta == 0] = 0
            g$d2 - p * crossprod(x$wo,
                                 tau.r2^(p-1) - tau.x1^(p-1) + xdp * Delta)[1,]
          }  )
    }
  }
  g
}

# Second gradient function

grad2 = function(eta, h, x, expdH=NULL, order=0) {
  g = vector("list", length(order))
  names(g) = paste("d", order, sep="")
  g[1:length(g)] = 0
  if(length(eta) == 0) return(NULL)
  m = length(eta)
  n1 = length(x$t)
  p = h$deg
  if(n1 > 0) {
    eta.r1 = rep(eta, rep.int(n1,m))
    dim(eta.r1) = c(n1,m)
    t.eta = pmax(x$t - eta.r1, 0)
    ind = x$t >= eta.r1
    if(0 %in% order) {
      g$d0  = switch(as.character(p),
          "0" = (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
              ind[x$i1,,drop=FALSE]) - crossprod(x$wt, t.eta))[1,],
          "1" = (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
              t.eta[x$i1,,drop=FALSE]) - 0.5 * crossprod(x$wt, t.eta^2))[1,],
          "2" = (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
              t.eta[x$i1,,drop=FALSE]^2) -
                    crossprod(x$wt, t.eta * t.eta * t.eta) / 3)[1,],
          (crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
                     t.eta[x$i1,,drop=FALSE]^p) -
           crossprod(x$wt, t.eta^(p+1)) / (p+1))[1,]  )
    }
    if(1 %in% order) {
      g$d1 = switch(as.character(p),
          "0" = double(m),
          "1" = (-crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h),
                            ind[x$i1,,drop=FALSE]) +
                 crossprod(x$wt, t.eta))[1,],
          "2" = (-2 * crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h), 
                                t.eta[x$i1,,drop=FALSE]) +
                 crossprod(x$wt, t.eta^2))[1,],
          { t.eta.pm1 = t.eta[x$i1,,drop=FALSE]^(p-1)
            if(p < 1) t.eta.pm1[t.eta.pm1 == Inf] = 0
            (-p * crossprod(x$wt[x$i1] / hazuh(x$t[x$i1], h), t.eta.pm1) +
              crossprod(x$wt, t.eta^p))[1,]
          }  )
    }
    if(2 %in% order) {
      g$d2 = switch(as.character(p),
          "0" = double(m),
          "1" = - crossprod(x$wt, ind)[1,],
          "2" = - 2 * crossprod(x$wt, t.eta)[1,],
          { t.eta.pm1 = t.eta^(p-1)
            if(p < 1) t.eta.pm1[t.eta.pm1 == Inf] = 0
            - p * crossprod(x$wt, t.eta.pm1)[1,]
          }  )
    }
  }
  n2 = nrow(x$o)
  if(n2 > 0) {
    if(is.null(expdH)) expdH = exp(chazuh(x$o[,1],h) - chazuh(x$o[,2],h))
    Delta = expdH / (1 - expdH)
    eta.r2 = rep(eta, rep.int(n2,m))
    dim(eta.r2) = c(n2,m)
    x1.eta = pmax(x$o[,1] - eta.r2, 0)
    x2.eta = pmax(x$o[,2] - eta.r2, 0)
    ind1 = x$o[,1] >= eta.r2
    ind2 = x$o[,2] >= eta.r2
    if(0 %in% order) {
      xd0 = switch(as.character(p),
          "0" = (x1.eta.p1 <- x1.eta) - x2.eta,
          "1" = (x1.eta.p1 <- x1.eta * x1.eta) - x2.eta * x2.eta,
          "2" = (x1.eta.p1 <- x1.eta * x1.eta * x1.eta) -
            x2.eta * x2.eta * x2.eta,
          (x1.eta.p1 <- x1.eta^(p+1)) -  x2.eta^(p+1) )
      xd0[Delta == 0] = 0
      g$d0 = g$d0 - crossprod(x$wo, x1.eta.p1 + xd0 * Delta)[1,] / (p+1)
    }
    if(1 %in% order) {
      xd1 = switch(as.character(p),
          "0" = (x1.eta.p <- ind1) - ind2,
          "1" = (x1.eta.p <- x1.eta) - x2.eta,
          "2" = (x1.eta.p <- x1.eta * x1.eta) - x2.eta * x2.eta,
          (x1.eta.p <- x1.eta^p) - x2.eta^p  )
      xd1[Delta == 0] = 0
      g$d1 = g$d1 + crossprod(x$wo, x1.eta.p + xd1 * Delta)[1,]
    }
    if(2 %in% order) {
      g$d2 = switch(as.character(p),
          "0" = double(m),
          "1" = g$d2 - crossprod(x$wo, ind1 + (ind1 - ind2) * Delta)[1,],
          "2" = { xd2 = x1.eta - x2.eta
                  xd2[Delta == 0] = 0
                  g$d2 - 2 * crossprod(x$wo, x1.eta + xd2 * Delta)[1,] },
          { x1.eta.pm1 = x1.eta^(p-1)
            x2.eta.pm1 = x2.eta^(p-1)
            if(p < 1) {
              x1.eta.pm1[x1.eta.pm1 == Inf] = 0
              x2.eta.pm1[x2.eta.pm1 == Inf] = 0
            }
            xdp = x1.eta.pm1 - x2.eta.pm1
            xdp[Delta == 0] = 0
            g$d2 - p * crossprod(x$wo, x1.eta.pm1 + xdp * Delta)[1,]
          }  )
    }
  }
  g
}

logLikuh = function(h, data) {
  x = icendata(data)
  if(length(x$t) > 0) ll = sum(x$wt[x$i1] * log(hazuh(x$t[x$i1], h))) -
      sum(x$wt * chazuh(x$t, h))
  else ll = 0
  if(nrow(x$o) > 0) {
    ll = ll + sum(x$wo * 
      (log(exp(-chazuh(x$o[,1], h)) - exp(-chazuh(x$o[,2], h)))))
  }
  ll
}

dlogLik = function(h, x, expdH=NULL, interior=FALSE) {
  if(interior) eta = h$eta[h$eta != h$upper]
  else eta = h$eta
  m = length(eta)
  d = c(grad0(h, x, expdH), grad1(h$tau, h, x, expdH)$d0,
    if(m>0) grad2(eta, h, x, expdH)$d0 else NULL)
  names(d) = c("alpha", paste("nu",1:length(h$tau),sep=""),
         if(m>0) paste("mu",1:m,sep="") else NULL)
  d
}

### Finds all local maxima of the gradient functions

# bc   boundary change

maxgrad = function(h, x, expdH=NULL, maxit=100, grid=100, bc=TRUE) {
  if(!is.icendata(x)) x = icendata(x, w=1)
  u = sort(unique(c(x$u, h$tau, h$eta)))
  p = h$deg
  if(p > 0 && p < 0.1)
    u = rep(u, rep.int(21, length(u))) + seq(-h$upper*1e-2, h$upper*1e-2, len=21)
  if(p == 1) maxit = 1
  tau1 = c(0, h$tau[h$tau > 0])
  k = length(tau1) - 1
  eta1 = c(h$eta[h$eta<h$upper], h$upper)
  m = length(eta1) - 1
  if(p < 1) grid = grid * log2(2/(p*p))

  ## grad1
  if(p < 0.1) {                     # optimal points
    if(bc) pt1 = u[u<=eta1[1]]
    else pt1 = u[u<eta1[1]]
    gpt1 = grad1(pt1, h, x, expdH)$d0
  }
  else {
    if(p == 1) u1 = u[u <= eta1[1]]
    else u1 = seq(0, eta1[1], len=grid)
    if(!bc) u1 = u1[u1 < eta1[1]]
    m1 = length(u1)
    pt1 = gpt1 = numeric()
    if(p == 1) jd = rep(TRUE, m1-1)
    else {
      g = grad1(u1, h, x, expdH, order=1)
      jd = g$d1[-m1] > 0 & g$d1[-1] < 0
    }
    if(any(jd)) {
      left = u1[-m1][jd]
      right = u1[-1][jd]
      pt1 = (left + right) * .5
      for(i in 1:maxit) {
        g = grad1(pt1, h, x, expdH, order=1:2)
        left[g$d1>0] = pt1[g$d1>0]
        right[g$d1<0] = pt1[g$d1<0]
        pt1.old = pt1
        pt1 = pt1 - g$d1 / g$d2
        j = is.na(pt1) | pt1 < left | pt1 > right
        pt1[j] = (left[j] + right[j]) * .5
        if( max(abs(pt1 - pt1.old)) <= 1e-14 * h$upper) break
      }
      if(p == 1) pt1 = pt1[!j]
      gpt1 = grad1(pt1, h, x, expdH, order=0)$d0
    }
  }
  i = pt1 > 0 & pt1 <= tau1[k+1]
  pt1i = pt1[i]
  gpt1i = gpt1[i]
  if(k > 0 && length(pt1i) > 0) {
    grp = apply(outer(pt1i, tau1[-(k+1)], ">") & outer(pt1i, tau1[-1], "<="),
        1, which.max)
    r1 = aggregate(gpt1i, by=list(group=grp), which.max)
    j = integer(k)
    j[r1[,1]] = r1[,2]
    j = j + c(0,cumsum(tabulate(grp, nbins=k))[-k])
    np1 = pt1i[j]
    gnp1 = gpt1i[j]
    j0 = gnp1 > 0
    np1 = np1[j0]
    gnp1 = gnp1[j0]
  }
  else np1 = gnp1 = numeric()

  ## grad2
  if(p < 0.1) {
    if(bc) pt2 = u[u>=tau1[k+1] & u<=h$upper]
    else pt2 = u[u>tau1[k+1] & u<h$upper]
    gpt2 = grad2(pt2, h, x, expdH)$d0
  }
  else {
    if(p == 1) u2 = u[u >= tau1[k+1]]
    else u2 = seq(tau1[k+1], h$upper, len=grid)
    if(!bc) u2 = u2[u2 > tau1[k+1]]
    m2 = length(u2)
    pt2 = gpt2 = numeric()
    if(p == 1) jd = rep(TRUE, m2-1)
    else {
      g = grad2(u2, h, x, expdH, order=1)
      jd = g$d1[-m2] > 0 & g$d1[-1] < 0
    }
    if(any(jd)) {
      left = u2[-m2][jd]
      right = u2[-1][jd]
      pt2 = (left + right) * .5
      for(i in 1:maxit) {
        g = grad2(pt2, h, x, expdH, order=1:2)
        left[g$d1>0] = pt2[g$d1>0]
        right[g$d1<0] = pt2[g$d1<0]
        pt2.old = pt2
        pt2 = pt2 - g$d1 / g$d2
        j = is.na(pt2) | pt2 < left | pt2 > right
        pt2[j] = (left[j] + right[j]) * .5
        if( max(abs(pt2 - pt2.old)) <= 1e-14 * h$upper) break
      }
      if(p == 1) pt2 = pt2[!j]
      gpt2 = grad2(pt2, h, x, expdH, order=0)$d0
    }
  }
  i = pt2 >= eta1[1] & pt2 < eta1[m+1]
  pt2i = pt2[i]
  gpt2i = gpt2[i]
  if(m > 0 && length(pt2i) > 0) {
    grp = apply(outer(pt2i, eta1[-(m+1)], ">=") & outer(pt2i, eta1[-1], "<"), 1,
        which.max)
    r2 = aggregate(gpt2i, by=list(group=grp), which.max)
    j = integer(m)
    j[r2[,1]] = r2[,2]
    j = j + c(0, cumsum(tabulate(grp, nbins=m))[-m])
    np2 = pt2i[j]
    gnp2 = gpt2i[j]
    j0 = gnp2 > 0 
    np2 = np2[j0]
    gnp2 = gnp2[j0]
  }
  else np2 = gnp2 = numeric()

  ## grad1 and grad2
  if(max(h$tau) != h$eta[1]) {
    jj1 = pt1 >= tau1[k+1] & pt1 <= eta1[1]
    
    if(p == 0) { uj1 = pt1[jj1]; gj1 = gpt1[jj1] }
    else {
      uj1 = c(tau1[k+1], if(bc) eta1[1] else NULL, pt1[jj1])
      gj1 = c(grad1(c(tau1[k+1], if(bc) eta1[1] else NULL),
          h, x, expdH)$d0, gpt1[jj1])
    }
    jmax = which.max(gj1)
    np31 = uj1[jmax]
    gnp31 = gj1[jmax]
    
    jj2 = pt2 >= tau1[k+1] & pt2 <= eta1[1]
    if(p == 0) { uj2 = pt2[jj2]; gj2 = gpt2[jj2] }
    else {
      uj2 = c(if(bc) tau1[k+1] else NULL, eta1[1], pt2[jj2])
      gj2 = c(grad2(c(if(bc) tau1[k+1] else NULL, eta1[1]),
          h, x, expdH)$d0, gpt2[jj2])
    }
    jmax = which.max(gj2)
    np32 = uj2[jmax]
    gnp32 = gj2[jmax]
    
    if(gnp31 > gnp32) {np1 = c(np1, np31); gnp1 = c(gnp1, gnp31)}
    else {np2 = c(np2, np32); gnp2 = c(gnp2, gnp32)}
  }
  
  list(np1=np1, gnp1=gnp1, np2=np2, gnp2=gnp2)
}

simplify = function(h) {
  i1 = order(h$tau)                     # remove identical knots
  tau = h$tau[i1]
  nu = h$nu[i1]
  i2 = order(h$eta)
  eta = h$eta[i2]
  mu = h$mu[i2]
  if(h$deg != 0 || length(tau) == 0 || length(eta) == 0) return(h)
  if(tau[length(tau)] != eta[1]) return(h)
  if(nu[length(tau)] < mu[1]) {
    if(eta[1] == 0) {
      h$alpha = h$alpha + mu[1]
      eta = eta[-1]
      mu = mu[-1]
    }
    else {
      h$alpha = h$alpha + nu[length(nu)]
      mu[1] = mu[1] - nu[length(nu)]
      tau = tau[-length(tau)]
      nu = nu[-length(nu)]
    }
  }
  else {
    h$alpha = h$alpha + mu[1]
    nu[length(nu)] = nu[length(nu)] - mu[1]
    eta = eta[-1]
    mu = mu[-1]
  }
  uh(h$alpha, tau, nu, eta, mu, h$upper, h$deg, collapse=FALSE)
}

# Collapse similar knots

collapse = function(h, x, tol=0){
  ll = attr(h, "ll")
  i1 = order(h$tau)                     # remove identical knots
  tau = h$tau[i1]
  ind1 = cumsum(!duplicated(c(tau)))
  tau = unique(tau)
  nu = aggregate(h$nu[i1], by=list(group=ind1), sum)[,2]
  nu[tau == 0] = 0
  i2 = order(h$eta)
  eta = h$eta[i2]
  ind2 = cumsum(!duplicated(eta))
  eta = unique(eta)
  mu = aggregate(h$mu[i2], by=list(group=ind2), sum)[,2]
  mu[eta == h$upper] = 0
  h = uh(h$alpha, tau, nu, eta, mu, h$upper, h$deg, collapse=FALSE)
  if(tol > 0) {
    if(is.null(ll)) ll = logLikuh(h, x)
    # if(h$deg < 1) {attr(h, "ll") = ll; return(h)}  # why?
    h2 = h
    break1 = break2 = FALSE
    repeat {
      if(!break1 && length(h2$nu) > 1) {
        j = which.min(diff(h$tau))
        h2$nu[j] = h2$nu[j] + h2$nu[j+1]
        h2$nu = h2$nu[-(j+1)]
        h2$tau[j] = (h2$tau[j] + h2$tau[j+1]) * .5
        h2$tau = h2$tau[-(j+1)]
        ll2 = logLikuh(h2, x)
        if(ll2 + tol >= ll) {h = h2; ll = ll2; break1 = FALSE}
        else {h2 = h; break1 = TRUE}
      }
      else break1 = TRUE
      if(!break2 && length(h2$mu) > 1) {
        j = which.min(diff(h2$eta))
        h2$mu[j] = h2$mu[j] + h2$mu[j+1]
        h2$mu = h2$mu[-(j+1)]
        h2$eta[j] = (h2$eta[j] + h2$eta[j+1]) * .5
        h2$eta = h2$eta[-(j+1)]
        ll2 = logLikuh(h2, x)
        if(ll2 + tol >= ll) {h = h2; ll = ll2; break2 = FALSE}
        else {h2 = h; break2 = TRUE}
      }
      else break2 = TRUE
      if(break1 && break2) break
    }
    attr(h, "ll") = ll
  }
  h
}

# deg - polynomial degree

uh = function(alpha, tau, nu, eta, mu, upper=Inf, deg=1, collapse=TRUE) {
  if(length(tau) == 0) {tau=0; nu=0}
  if(length(eta) == 0) {eta=upper; mu=0}
  i1 = order(tau)
  tau = tau[i1]
  nu = nu[i1]
  i2 = order(eta)
  eta = eta[i2]
  mu = mu[i2]
  h = list(alpha=alpha, tau=tau, nu=nu, eta=eta, mu=mu, upper=upper, deg=deg)
  if(collapse) h = collapse(h)
  class(h) = "uh"
  h
}

print.uh = function(x, ...) {
  cat("$alpha\n")
  print(x$alpha, ...)
  print(cbind(tau=x$tau, nu=x$nu), ...)
  print(cbind(eta=x$eta, mu=x$mu), ...)
  cat("$upper\n")
  print(x$upper, ...)
  cat("$deg\n")
  print(x$deg, ...)
}

# Hazard function

hazuh = function(t, h) {
  p = h$deg
  b = c = 0
  if(length(h$tau) > 0) {
    if(p == 0) d = outer(h$tau, t, ">=")
    else {
      d = pmax(outer(h$tau, t, "-"), 0)
      if(p != 1) d = d^p
    }
    b = (h$nu %*% d)[1,]
  }
  if(length(h$eta) > 0) {
    if(p == 0) d = outer(t, h$eta, ">=")
    else {
      d = pmax(outer(t, h$eta, "-"), 0)
      if(p != 1) d = d^p
    }
    c = (d %*% h$mu)[,1]
  }
  h$alpha + pmax(b, c)
}

chazuh = function(t, h) {
  deg = pmax(h$deg, 1)
  p1 = h$deg + 1
  a = b = c = 0
  if(h$alpha > 0) a = h$alpha * t
  if(length(h$tau) > 0) {
    tau.t = pmax(outer(h$tau, t, "-"), 0)
    b = (h$nu %*% (h$tau^p1  - tau.t^p1))[1,] / p1
  }
  if(length(h$eta) > 0) {
    t.eta = pmax(outer(t, h$eta, "-"), 0)
    c = (t.eta^p1 %*% h$mu)[,1] / p1
  }
  H = a + b + c
  H[t > h$upper] = 1e100
  H
}

# survival function

survuh = function(t, h) exp(-chazuh(t, h))

# density function

denuh = function(t, h) hazuh(t, h) * survuh(t, h)

## plotting functions

plot.Uhaz = function(x, ...) plot(x$h, ...)

plot.uh = function(x, data, fn=c("haz","grad","surv","den","chaz"), ...) {
  fn = match.arg(fn)
  fnR = getFunction(paste("plot",fn,"uh",sep=""))
  switch(fn,
         "haz" =,
         "surv" =,
         "den" =,
         "chaz" = fnR(x, ...),
         "grad" = fnR(x, data, ...)  )
}

plothazuh = function(h, add=FALSE, col="darkblue", lty=1,
    xlim, ylim, lwd=2, pch=19, len=500, vert=FALSE, add.knots=TRUE,
    xlab="Time", ylab="Hazard", ...) {
  p = h$deg
  pt = switch(as.character(p),
      "0" = unique(sort(c(0, h$tau, h$eta, h$upper))),
      "1" = unique(sort(c(0, h$tau, h$eta, h$upper))),
      unique(sort(c(h$tau, h$eta, seq(0, h$upper, len=len)))))
  m = length(pt)
  knots = unique(c(h$tau, h$eta)) 
  haz = hazuh(pt, h)
  max.haz = max(haz)
  if(missing(xlim)) xlim = range(pt)
  if(missing(ylim)) ylim = c(0, max.haz)
  if(!add) plot(pt, haz, type="n", xlim=xlim, ylim=ylim,
                xlab=xlab, ylab=ylab, ...)
  if(vert) {
    lines(rep(max(h$tau),2), ylim, col="grey", lty=2)
    lines(rep(min(h$eta),2), ylim, col="grey", lty=2)
  }
  abline(h=0, col ="black")

  if(p == 0) {
    lines(c(h$tau[length(h$tau)], h$eta[1]), rep(h$alpha,2),
          lwd=lwd, col=col, lty=lty)
    lines(c(rep(rev(h$tau),each=2), 0), 
          c(h$alpha, rep(hazuh(rev(h$tau), h), each=2)),
          lwd=lwd, col=col, lty=lty)
    lines(c(rep(h$eta,each=2), h$upper), 
          c(h$alpha, rep(hazuh(h$eta, h), each=2)),
          lwd=lwd, col=col, lty=lty)
  }
  else lines(pt, haz, lwd=lwd, col=col, lty=lty)
  if(add.knots && length(knots) > 0)
    points(knots, hazuh(knots, h), col=col, pch=pch)
}

plotchazuh = function(h, add=FALSE, lwd=2, len=500, col="darkblue",
    pch=19, add.knots=TRUE, vert=FALSE, xlim, ylim, ...) {
  pt = unique(sort(c(seq(0, h$upper, len=len), h$tau, h$eta)))
  m = length(pt)
  H = chazuh(pt, h)
  max.H = max(H)
  if(missing(xlim)) xlim = range(pt)
  if(missing(ylim)) ylim = c(0, max.H)
  plot(rep(max(h$tau),2), c(0,max.H), type="n",
       xlim=xlim, ylim=ylim, xlab="Time", ylab="Cumulative Hazard", ...)
  if(vert) lines(rep(max(h$tau),2), c(0,max.H), type="l", col ="grey", lty=2)
  if(vert) lines(rep(min(h$eta),2), c(0,max.H), col ="grey", lty=2)
  abline(h=0, col ="black")
  lines(pt, H, type="l", lwd=lwd, col=col)
  if(add.knots) {
    knots = c(h$tau, h$eta)
    knots = knots[knots>0 & knots<h$upper]
    points(knots, chazuh(knots, h), col=col, pch=pch)
  }
}

plotdenuh = function(h, add=FALSE, lty=1, lwd=2,
    col="darkblue", add.knots=TRUE, pch=19,
    ylim, len=500, vert=FALSE, ...) {
  pt = sort(c(seq(0, h$upper, len=len), h$tau, h$eta))
  m = length(pt)
  d = denuh(pt, h)
  max.d = max(d)
  type = if(vert) "l" else "n"
  if(missing(ylim)) ylim = c(0, max.d)
  if(add) lines(rep(max(h$tau),2), ylim, type=type, col ="grey", lty=2)
  else plot(rep(max(h$tau),2), ylim, type=type,
            col ="grey", lty=2, xlim=range(pt), ylim=ylim,
            xlab="Time", ylab="Density", ...)
  if(vert) lines(rep(min(h$eta),2), ylim, col="grey", lty=2)
  abline(h=0, col ="black")
  lines(pt, d, lty=lty, lwd=lwd, col=col)
  if(add.knots) {
    knots = c(h$tau, h$eta)
    knots = knots[knots>0 & knots<h$upper]
    points(knots, denuh(knots, h), col=col, pch=pch)
  }
}

plotsurvuh = function(h, add=FALSE, lty=1, lwd=2, len=500, vert=FALSE,
    col="darkblue", pch=19, add.knots=TRUE, xlim, ylim, ...) {
  pt = sort(c(seq(0, h$upper, len=len), h$tau, h$eta))
  m = length(pt)
  s = survuh(pt, h)
  if(missing(xlim)) xlim = range(pt)
  if(missing(ylim)) ylim = c(0, 1)
  if(add) {
    if(vert) lines(rep(max(h$tau),2), c(0,1), col ="grey", lty=2)
  }
  else {
    if(vert)
      plot(rep(max(h$tau),2), c(0,1), type="l", col ="grey", lty=2,
           xlim=xlim, ylim=ylim, xlab="Time", ylab="Survival Probability",
           ...)
    else plot(0, 0, type="n", xlim=xlim, ylim=ylim,
              xlab="Time", ylab="Survival Probability")
  }
  if(vert) lines(rep(min(h$eta),2), c(0,1), col ="grey", lty=2)
  abline(h=0, col ="black")
  lines(pt, s, lty=lty, lwd=lwd, col=col)
  pt = c(h$tau, h$eta)
  if(add.knots) {
    knots = c(h$tau, h$eta)
    knots = knots[knots>0 & knots<h$upper]
    points(knots, survuh(knots, h), col=col, pch=pch)
  }
}

plotgraduh = function(h, data, w=1, len=500, xlim, ylim, vert=TRUE,
    add=FALSE, xlab="Time", ylab="Gradient",
    col0="red3", col1="blue3", col2="green3", order=0, ...) {
  x = icendata(data, w)
  maxf = h$upper 
  if(missing(xlim)) xlim = c(0, maxf)
  k = length(h$tau)
  tau1 = unique(sort(c(seq(xlim[1], h$eta[1] - if(h$deg!=0) 0 else 1e-4*h$upper,
      len=len), h$tau)))
  eta1 = unique(sort(c(seq(h$tau[k] + if(h$deg!=0) 0 else 1e-4*h$upper, xlim[2],
      len=len), h$eta)))
  g11 = grad1(tau1, h, x, order=order)[[1]]
  g22 = grad2(eta1, h, x, order=order)[[1]]
  g0 = grad0(h, x)
  if(missing (ylim)) ylim= range(g0, g11, g22)
  if(add) lines(c(0,maxf), c(0, 0), col="grey")
  else plot(c(0,maxf), c(0, 0), type="l", col="grey", xlim=xlim, ylim=ylim,
            xlab=xlab, ylab=ylab, ...)
  if(vert) {
    lines(rep(h$tau[k],2), ylim, col="gray", lty=2)
    lines(rep(h$eta[1],2), ylim, col="gray", lty=2)
  }
  lines(eta1, g22, col=col2, lty=1, lwd=1)
  lines(tau1, g11, col=col1, lty=1, lwd=1)
  j1 = h$tau != 0 & h$tau != h$upper
  j2 = h$eta != 0 & h$eta != h$upper
  points(h$tau[j1], grad1(h$tau[j1],h,x)$d0, col=col1, pch=2)
  points(h$eta[j2], grad2(h$eta[j2],h,x)$d0, col=col2, pch=6)
  lines(c(h$tau[k], h$eta[1]), rep(g0, 2), col=col0, lty=1, lwd=1) 
  points((h$tau[k] + h$eta[1])/2, g0, col=col0, pch=20)
  if(!add) legend("bottomright",
                  legend=expression(d[0], d[1], d[2]), col=c(col0, col1, col2),
                  lty=c(1,1,1), pch=c(20,2,6), bty="n", lwd=c(1,1,1))
}

uh.initial = function(x, deg=1) {
  if(!is.icendata(x)) x = icendata(x)
  n1 = length(x$t)
  n2 = nrow(x$o)
  y1 = if(n1 > 0) x$t else numeric()
  y2 = if(n2 > 0) {
    x$o[x$o[,2] == Inf,2] = 1.6 * x$upper
    rowMeans(x$o)
  } 
  else numeric()
  beta = weighted.mean(c(y1,y2), c(x$wt, x$wo))
  uh(alpha=1/beta, tau=NULL, nu=NULL, eta=NULL, mu=NULL, upper=x$upper,
     deg=deg)
}

