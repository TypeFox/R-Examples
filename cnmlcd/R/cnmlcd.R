cnmlcd = function(x, lcd, maxit = 100, tol = 1e-7, quiet = TRUE) {
  X =  x.weight(x)
  x = X$x
  w = X$w
  n = length(x)

  low = x[1]
  up = x[n]

  if(n < 2) {
    lcd = create.lcd(alpha = 0, L = low, U = up)
    return(list(lcd = lcd, ll = logLik.lcd(object = lcd, x, w = w),
                num.iterations = 0, max.gradient = 0, convergence = 0))
  }

  attr(x, "xx") = rev(cumsum(rev(w))) * x -  rev(cumsum(rev(x * w)))

  if(missing(lcd)) {
    k = floor(seq(2, n - 1, length = min(n, 5)))
    p = 0
    theta = x[unique(k)]
    lcd = create.lcd(alpha = 0, theta = theta, p = p, L = low, U = up)
  }
  ll.lcd = logLik.lcd(object = lcd, x, w = w)
  convergence = 1
  ll.lcd1 = -Inf

  for(i in 1:maxit) {
    if(ll.lcd <= ll.lcd1  + tol) {convergence = 0; break}
    lcd1 = lcd
    ll.lcd1 = ll.lcd
    if(!quiet) {
      print(ll.lcd, 16)
      plotgradient(lcd, x, w = w)
    }

    g = maxima.gradient(lcd, x, w = w, tol = 0)
    if(length(g$theta) != 0) {
      nsp = g$theta
      nsp = nsp[! nsp %in% c(lcd$L, lcd$theta, lcd$U)]
      nsl = length(nsp)
      if(nsl >= 1)
        lcd = create.lcd(lcd$alpha, c(lcd$theta, nsp),
                         c(lcd$p, rep(0, nsl)), lcd$L, lcd$U)
    }
    R = cgn.s(x, lcd, ll.lcd = ll.lcd, w = w)
    lcd = R$lcd
    ll.lcd = R$ll
    if(length(lcd$theta) >= 1) {
      j = lcd$p > 0
      t = lcd$theta[j]
      p = lcd$p[j]
      lcd = create.lcd(alpha = lcd$alpha, theta = t, p = p, L = lcd$L,
          U = lcd$U)
    }
  }
  if(!quiet) {
    print(ll.lcd, 16)
    plotgradient(lcd, x, w = w)
  }
  g = lcd.gradient(lcd, x, w = w, theta = lcd$theta)
  list(lcd = lcd, ll = ll.lcd, num.iterations = i,
       max.gradient = max(g), convergence = convergence)
}

x.weight = function (x) {
  x = sort(x)
  n = length(x)
  if(n > 1) {
    i = which(diff(x) > 0)
    i.n = c(i, n)
    w = i.n - c(0, i)
  }
  else w = i.n = 1
  list(w = w, x = x[i.n])
}

## Creates an 'lcd' object
create.lcd = function(alpha, theta = NULL, p = NULL, L, U) {
  if(length(theta) >  0) {
    k = max(length(theta), length(p), na.rm = TRUE)
    theta = rep(theta, len = k)
    p = rep(p, len = k)
    index = order(theta)
    theta = theta[index]
    p = p[index]
  }
  s = alpha - c(0, cumsum(p))
  c = c(0, cumsum(theta * p)) - alpha * L
  tu = c(theta, U)
  tl = c(L, theta)
  Int = integration.lcd(s, c, tu, tl, Order = 2)
  C0 = sum(Int$C)
  Int$C = Int$C / C0
  Int$Cr = Int$Cr / C0
  Int$Cr2 = Int$Cr2 / C0
  Par = list(s = s, c = c, tu = tu, tl = tl, Int = Int)
  lcd = list(alpha = alpha, C = C0, theta = theta, p = p, L = L, U = U)
  attr(lcd, "Par") = Par
  class(lcd) = "lcd"
  lcd
}

## Calculates some integrations
integration.lcd = function(s, c, tu, tl, Order = 0) {
  m = length(s)
  pp = pr = pr2 = double(m)

  j = s == 0
  if(any(j)) {
    E = exp(c[j])
    pp[j] = E *(tu[j] - tl[j])
    if(Order >= 1) pr[j] = E * (tu[j]^2 - tl[j]^2) / 2
    if(Order == 2) pr2[j] = E * (tu[j]^3 - tl[j]^3) / 3
  }

  j2 = !j
  if(any(j2)) {
    x1 = exp(s[j2] * tu[j2] + c[j2]) / s[j2]
    y1 = exp(s[j2] * tl[j2] + c[j2]) / s[j2]
    pp[j2] = x1 - y1
    if(Order >= 1) pr[j2] = tu[j2] * x1 - tl[j2] * y1 - pp[j2] / s[j2]
    if(Order == 2) pr2[j2] = tu[j2]^2 * x1 - tl[j2]^2 * y1 - 2 * pr[j2] / s[j2]
  }
  list(C = pp, Cr = pr, Cr2 = pr2)
}

## the density function return f(x,G)
dlcd = function(x, lcd, log = FALSE) {
  late.part = 0
  nt = length(lcd$theta)
  if(nt > 0) {
    n = length(x)
    xt = pmax(x - rep(lcd$theta, each = n), 0)
    dim(xt) = c(n, nt)
    late.part = xt %*% lcd$p
  }
  index = lcd$alpha * x - late.part - lcd$alpha * lcd$L
  index[x < lcd$L | x > lcd$U] = -Inf
  logd = drop(index - log(lcd$C))
  if(log) logd else exp(logd)
}

## Log-likelihood function
logLik.lcd = function(object, x, w = NULL, ...) {
  if(is.null(w)) {
    X = x.weight(x)
    x = X$x
    w = X$w
  }
  sum(dlcd(x, lcd = object, log = TRUE) * w)
}


## the probability function
plcd = function(q, lcd, lower.tail = TRUE, log.p = FALSE) {
  n = length(q)
  Par = attr(lcd, "Par")
  nkots = length(lcd$theta)
  if(nkots > 0) {
    M = rbind(cbind(lcd$theta, 1), cbind(q, 0))
    I = order(M[, 1])
    OM = M[I, ]
    k = which(OM[, 2] == 0) - 0:(n-1)
    k = k[order(order(q))]
  }
  else k = rep(1, n)
  s = Par$s[k]
  c = Par$c[k]
  tl = c(lcd$L, lcd$theta)[k]
  pb = integration.lcd(s, c, q, tl, Order = 0)$C / lcd$C
  if(nkots > 0) {
    c = cumsum(c(0, Par$Int$C))
    pb = pb + c[k]
   }
  pb = pb
  pb[q <= lcd$L] = 0
  pb[q >= lcd$U] = 1
  if(lower.tail) pb = pb
  else pb = 1 - pb
  if(log.p) log(pb)
  else pb
}

## the calculations for the gradient values.
prob.lcd = function(lcd, x) {
  n = length(x)
  Par = attr(lcd, "Par")
  nkots = length(lcd$theta)
  if(nkots > 0) {
    M = rbind(cbind(lcd$theta, 1), cbind(x, 0))
    I = order(M[, 1])
    OM = M[I, ]
    k = which(OM[, 2] == 0) - 0:(n-1)
    k = k[order(order(x))]
  }
  else k = rep(1, n)
  s = Par$s[k]
  c = Par$c[k]
  tl = c(lcd$L, lcd$theta)[k]
  F = integration.lcd(s, c, x, tl, Order = 1)
  pb = F$C / lcd$C
  pbr = F$Cr / lcd$C
  if(nkots > 0) {
    c = cumsum(c(0, Par$Int$C))
    cr = cumsum(c(0, Par$Int$Cr))
    pb = pb + c[k]
    pbr = pbr + cr[k]
   }
  pbr = sum(Par$Int$Cr) - pbr
  list(pb = 1 - pb, pbr = pbr)
}

## the gradient function
lcd.gradient = function(lcd, x, w = NULL, theta) {
  if(length(theta) < 1) return(NULL)
  if(is.null(w)) {
    X = x.weight(x)
    x = X$x
    w = X$w
  }

  if(!is.numeric(theta)) {
    H = attr(x,"xx")
    pbb = prob.lcd(lcd, x)
    theta = x
   }
  else{
    n = length(x)
    xt = pmax(x - rep(theta, each = n), 0)
    dim(xt) = c(n, length(theta))
    H = -colSums(xt * w)
    pbb = prob.lcd(lcd, theta)
  }
  H + sum(w) * (pbb$pbr - pbb$pb * theta)
}

maxima.gradient = function(lcd, x, w, tol = -Inf) {
  knots = c(lcd$L, lcd$theta, lcd$U)
  nk = length(knots)
  mk = match(knots, x)
  kf = mk[-nk]
  kb = mk[-1]
  k = c(rep(1:(nk - 1), kb - kf), nk - 1)

  Par = attr(lcd, "Par")
  s = Par$s[k]
  c = Par$c[k]
  tl = Par$tl[k]
  F = integration.lcd(s, c, x, tl, Order = 1)
  pb = 1 - F$C / lcd$C
  pbr = sum(Par$Int$Cr) - F$Cr / lcd$C
  if(nk > 2) {
    c = cumsum(c(0, Par$Int$C))
    cr = cumsum(c(0, Par$Int$Cr))
    pb = pb - c[k]
    pbr = pbr - cr[k]
   }

  l.C = attr(x, "xx") + sum(w) * (pbr - pb * x)
  p = g = double(nk - 1)
  for(i in 1:(nk - 1)) {
    ma = which.max(l.C[kf[i]:kb[i]])
    p[i] = x[ma + kf[i] - 1]
    g[i] = max(l.C[kf[i]:kb[i]])
  }
  index = g > tol
  list(theta = p[index], gradient = g[index])
}


## updates lcd$alpha,lcd$pi by using the constrained Newton methods

cgn.s = function(x, lcd, ll.lcd, maxit = 1, tol = 1e-7, w = NULL,
    lambda = 1e-15) {
  if(is.null(w)) {
    X = x.weight(x)
    x = X$x
    w = X$w
  }
  convergence = 1
  Par = attr(lcd, "Par")
  if(missing(ll.lcd)) ll.lcd = logLik.lcd(object = lcd, x, w = w)
  for(i in 1:maxit) {
    lcd1 = lcd
    ll.lcd1 = ll.lcd
    C = rev(Par$Int$C)
    Cr = rev(Par$Int$Cr)
    Cr0 = sum(Cr)
    Cr2 = rev(Par$Int$Cr2)
    Cr20 = sum(Cr2)
    n = sum(w)
    nx = length(x)
    nt = length(C)

    A0 = (x - Cr0) * w
    T = NULL
    T2 = matrix(ncol = nt,nrow = nt)
    T2[1,1] = Cr20 - (Cr0)^2

    if(nt > 1) {
      CC = rev(cumsum(C[-nt]))
      CCr = rev(cumsum(Cr[-nt]))
      CCr2 = rev(cumsum(Cr2[-nt]))
      P = CCr - CC * lcd$theta
      ntt = nt - 1
      xt = -pmax(x - rep(lcd$theta, each =  nx), 0)
      dim(xt) = c(nx, ntt)
      T = (xt + rep(P, each = nx)) * w
      T2[1,-1] = T2[-1,1] = P * Cr0 - CCr2 + CCr * lcd$theta

      mt = rep(lcd$theta, each =  ntt) + lcd$theta
      mt2 = lcd$theta %*% t(lcd$theta)
      dim(mt) = c(ntt, ntt)
      ut = upper.tri(mt)
      mt[ut] = mt2[ut] = 0
      mm = CCr2 - mt * CCr + mt2 * CC
      mm[ut] = 0
      mm1 = t(mm)
      diag(mm1) = 0
      T2[-1,-1] = mm + mm1 - P %*% t(P)
    }

    p = colSums(cbind(A0, T)) / n + drop(c(lcd$alpha, lcd$p) %*% T2)
    eq = eigen(T2)
    v2 = sqrt(eq$values[eq$values >= eq$values[1] * lambda])
    kr = length(v2)
    a = t(eq$vectors[,1:kr,drop = FALSE]) * v2
    b = colSums(eq$vectors[,1:kr,drop = FALSE] * p / rep(v2, each = length(p)))
    R = pnnls(a,b,1)

    lcd$alpha = R$x[1]
    lcd$p = R$x[-1]
    R = line.lcd(lcd1, lcd, x, w = w, ll.lcd1 = ll.lcd1)
    lcd = R$lcd
    ll.lcd = R$ll.lcd
    ## if(ll.lcd < ll.lcd1 + tol) {convergence = 0; break}
  }
  list(lcd = lcd, ll.lcd = ll.lcd, I = i, convergence = convergence)
}


## Line search
line.lcd = function(lcd1, lcd2, x, w = NULL, ll.lcd1, tol = 1e-10) {
  ll.k = function(k) {
    m = create.lcd((1 - k) * lcd1$alpha + k * lcd2$alpha, lcd2$theta,
        (1 - k) * lcd1$p + k * lcd2$p, lcd2$L, lcd2$U)
    ll = logLik.lcd(object = m, x, w = w)
    list(lcd = m, ll = ll)
  }
  convergence = 0
  k = 1
  repeat{
    new = ll.k(k)
    if(new$ll >= ll.lcd1) break
    if(k < tol) {convergence = 1; break}
    k = k / 2
  }
  list(lcd = new$lcd, ll.lcd = new$ll, k = k, convergence = convergence)
}


## Print an 'lcd' object
print.lcd = function(x, ...) {
  a = cbind(alpha = x$alpha, C = x$C)
  print(a, ...)
  if(length(x$theta) > 0) {
    b = cbind(theta = x$theta, p = x$p)
    print(b, ...)
  }
  c = cbind(L = x$L, U = x$U)
  print(c, ...)
}

## plot gradient function
plotgradient = function(lcd, x, w = NULL) {
  if(is.null(w)) {
    X = x.weight(x)
    x = X$x
    w = X$w
  }
  if(is.null(attr(x, "xx")))
    attr(x,"xx") = rev(cumsum(rev(w))) * x -  rev(cumsum(rev(x * w)))

  plot(x, lcd.gradient(lcd, x, w = w, theta = "x"), type = "l",
       ylab = "Gradient", xlab = expression(theta), cex.axis= 1.5, cex.lab = 1.5)
  abline(h = 0, col = 'lightgrey', lty = 1)
  points(c(lcd$L,lcd$theta,lcd$U),
         lcd.gradient(lcd, x, w = w, c(lcd$L, lcd$theta, lcd$U)),
         col = 'red', pch = 16)
}


## plot density function
plot.lcd = function(x, data, log = FALSE, line.col = "blue", knot.col = "red", ...) {
  p = sort( c(seq(x$L, x$U, len = 200), x$theta) )
  sp = c(x$L, x$theta, x$U)
  np = length(p)
  d = dlcd(p, lcd = x)
  sd = dlcd(sp, lcd = x)

  if(!missing(data)) h = hist(data, breaks = 30, plot = FALSE)
  else h = hist(p, breaks = 30, plot = FALSE)

  if(log) {
    d = log(d)
    ylim = c(min(d), max(d))
    if(missing(data)) {
      plot(p, d, xlab = "Data", ylab = "Log-density", ylim = ylim,
           main = "Log density estimation",
           type = "l", lwd = 2, col = line.col, ...)
      base = min(d)
    }
    else{
      logDensity = log(h$density)
      y0 = range(logDensity, finite = TRUE)
      ylim[1] = min(y0[1], ylim[1])
      ylim[2] = max(y0[2], ylim[2])
      base = ylim[1] - 0.25 * diff(ylim)
      loghist(data, h = h, main = "Log density estimation", ylim = ylim,
              xlab = "Data", base = base, ...)
      lines(p, d, col = line.col, lwd = 2)
    }
    segments(p[1], base, p[1], d[1], col = line.col, lwd = 2)
    segments(p[np], base, p[np], d[np], col = line.col, lwd = 2)
    points(sp, log(sd), col = knot.col, pch = 16)
  }
  else{
    ylim = c(0, max(d, h$density))
    if(!missing(data)) {
      hist(data, breaks = 30, freq = F, ylim = ylim, xlab = "Data",
           main = "Log-concave density estimation", ...)
      lines(p, d, col = line.col, lty = 1, lwd = 2)

    }
    else {
      plot(p, d, col = line.col, type = "l", lwd = 2,
           xlab = "Data", ylab = "Density",
           main = "Log-concave density estimation", ...)
    }
    segments(p[1], 0, p[1], d[1], col = line.col, lwd = 2)
    segments(p[np], 0, p[np], d[np], col = line.col, lwd = 2)
    points(sp, sd, col = knot.col, pch = 16)
  }
}

loghist = function(x, h = NULL, breaks = 30, main = NULL, ylim = NULL, lty = 1,
    xlim = NULL, base = NULL, xlab = "x", ylab = "Log-density") {
  if(is.null(h)) h = hist(x, breaks = breaks, plot = FALSE)
  logd = log(h$density)
  if(is.null(ylim)) ylim = range(logd, finite = TRUE)
  if(is.null(base)) base = ylim[1] - 0.5 * diff(ylim)
  breaks = h$breaks
  if (is.null(xlim)) xlim = range(breaks)
  mids = h$mids
  nb = length(breaks)
  heights = rep(NA, nb)
  for (j in 2:(nb - 1))
    if (is.finite(max(logd[j - 1], logd[j])))
      heights[j] = max(logd[j - 1], logd[j])
  if(is.finite(logd[1])) heights[1] = logd[1]
  if(is.finite(logd[nb - 1])) heights[nb] = logd[nb - 1]
  plot(mids, logd, xlim = xlim, ylim = ylim, type = "n",
       xlab = xlab, ylab = ylab, main = main)
  for(i in 1:nb) {
    segments(breaks[i], logd[i], breaks[i + 1], logd[i], lty = lty)
    segments(breaks[i], heights[i], breaks[i], base, lty = lty)
    segments(breaks[nb], heights[nb], breaks[nb], base, lty = lty)
  }
}


