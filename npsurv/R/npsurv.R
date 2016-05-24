# ----------------------------------------------------------------------- #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #

npsurv = function(data, w=1, maxit=100, tol=1e-6, verb=0) {
  # plot = match.arg(plot)
  x2 = icendata(data, w)
  if(nrow(x2$o) == 0 || all(x2$o[,2] == Inf)) { # exact or right-censored only
    r0 = km(x2)
    r = list(f=r0$f, upper=max(x2$t, x2$o[,1]), convergence=TRUE, ll=r0$ll,
        maxgrad=0, numiter=1)
    class(r) = "npsurv"
    return(r)
  }
  x = rbind(cbind(x2$t, x2$t), x2$o)
  nx = nrow(x)
  w = c(x2$wt, x2$wo)
  wr = sqrt(w)
  n = sum(w)
  L = x[,1]
  R = x[,2]
  upper = x2$upper
  dmat = Deltamatrix(x)
  left = dmat$left
  right = dmat$right
  D = dmat$Delta
  m = length(left)
  p = double(m)
  i = rowSums(D) != 1
  j = colSums(D[!i,,drop=FALSE]) > 0
  j[c(1,m)] = TRUE
  repeat {                 # Initial p must ensure P > 0
    jm = which.max(colSums(D[i,,drop=FALSE]))
    j[jm] = TRUE
    i[D[,jm]] = FALSE
    if( sum(i) == 0 ){ p[j] = rep(1/sum(j), sum(j)); break }
  }
  P = (D %*% p)[,1]
  ll = sum( w * log(P) )
  converge = FALSE
  for(i in 1:maxit) {
    p.old = p
    ll.old = ll
    S = D / P
    d = crossprod(w, S)[1,]
    dmax = max(d) - n
    if(verb>0) { cat("##### Iteration", i, "#####\n")
                 cat("Log-likelihood: ", signif(ll, 6), "\n")
                 if(verb>1) cat("Maximum gradient: ", signif(dmax, 6), "\n")
                 if(verb>2) {cat("Probability vector:\n"); print(p)} }
    j[which(j)-1 + aggregate(d, by=list(group=cumsum(j)), which.max)[,2]] = TRUE
    pj = nnls(rbind(wr * (S[,j,drop=FALSE] - 2), rep(1,sum(j))),
        c(rep(0,nx),1))$x
    p[j] = pj / sum(pj)
    alpha = 1                # line search
    pd = p - p.old
    lld = sum(d * pd)
    p.alpha = p
    repeat {
      P.alpha = (D %*% p.alpha)[,1]
      ll.alpha = sum(w * log(P.alpha))
      if(ll.alpha >= ll + alpha * lld * .33)
        { p = p.alpha; P = P.alpha; ll = ll.alpha; break }
      alpha = alpha * .5
      if(alpha < 1e-10) break
      p.alpha = p.old + alpha * pd
    }
    j = p > 0
    if( ll <= ll.old + tol ) {converge=TRUE; break}
  }
  f = idf(left[j], right[j], p[j])
  r = list(f=f, upper=upper, convergence=converge, ll=ll,
      maxgrad=max(crossprod(w/P, D))-n, numiter=i)
  class(r) = "npsurv"
  r
}

# LR    matrix of intervals

Deltamatrix = function(LR) {
  L = LR[,1]
  R = LR[,2]
  ic = L != R          # inverval-censored
  nc = sum(ic)
  if(nc > 0) {
    L1 = L[ic] + max(R[R!=Inf]) * 1e-8          # open left endpoint
    R1 = R[ic]                                  # closed right endpoint
    LRc = cbind(c(L1, R1), c(rep(0,nc), rep(1,nc)), rep(1:nc, 2))
    LRc.o = LRc[order(LRc[,1]),]
    j = which(diff(LRc.o[,2]) == 1)
    left = L[ic][LRc.o[j,3]]
    right = R[ic][LRc.o[j+1,3]]
  }
  else left = right = numeric(0)
  if(nrow(LR) - nc > 0) {
    ut = unique(L[!ic])
    jin = colSums(outer(ut, left, ">") & outer(ut, right, "<=")) > 0
    left = c(ut, left[!jin])   # remove those that contain exact obs.
    right = c(ut, right[!jin])
    o = order(left, right)
    left = left[o]
    right = right[o]
  }
  D = outer(L, left, "<=") & outer(R, right, ">=")
  dimnames(D) = names(left) = names(right) = NULL
  list(left=left, right=right, Delta=D)
}

# interval distribution function, i.e., a distribution function defined on
# a set of intervals.

# left      Left endpoints of the intervals
# right     Right endpoints of the intervals
# p         Probability masses allocated to the intervals

idf = function(left, right, p) {
  if(length(left) != length(right)) stop("length(left) != length(right)")
  names(left) = names(right) = names(p) = NULL
  p = rep(p, length=length(left))
  f = list(left=left, right=right, p=p/sum(p))
  class(f) = "idf"
  f
}

print.idf = function(x, ...) {
  print(cbind(left=x$left, right=x$right, p=x$p), ...)
}

# Kaplan-Meier estimate of the survival function for right-censored data

km = function(data, w=1) {
  x = icendata(data, w)
  if(any(x$o[,2] != Inf))
    stop("Not all observations are exact or right-censored")
  if(nrow(x$o) == 0) {              # no right-censored observations
    f = idf(x$t, x$t, x$wt)
    ll = sum(x$wt * log(f$p))
    return(list(f=f, ll=ll))
  }
  c = colSums(x$wo * outer(x$o[,1], x$t, "<"))
  n = sum(x$wt, x$wo)                            # number of observations
  r = n - c - c(0,cumsum(x$wt))[1:length(x$t)]   # no. at risk
  S = cumprod(1 - x$wt/r)                        # survival prob.
  # tab = cbind(x$t, x$wt, c, r, S)
  p = rev(diff(rev(c(1,S,0))))
  dc = x$wt + c
  if(max(x$t) > max(x$o[,1])) {
    f = idf(x$t, x$t, p[-length(p)])
    ll = sum( x$wt * log(f$p) )
  }
  else {
    f = idf(c(x$t,max(x$o[,1])), c(x$t,Inf), p)
    ll = sum(c(x$wt, n - sum(x$wt)) * log(f$p))
  }
  list(f=f, ll=ll)
}


####  Plot functions

plot.npsurv = function(x, ...) plot(x$f, ...)

plot.idf = function(x, data, fn=c("surv","grad"),
    ...) {
  fn = match.arg(fn)
  fnR = getFunction(paste("plot",fn,"idf",sep=""))
  switch(fn, "surv" = fnR(x, ...), "grad" = fnR(x, data, ...)  )
}

plotgradidf = function(f, data, w=1, col1="red3", col2="blue3", 
    xlab="Survival Time", ylab="Gradient", xlim, ...) {
  x2 = icendata(data, w)
  x = rbind(cbind(x2$t, x2$t), x2$o)
  w = c(x2$wt, x2$wo)
  dmat = Deltamatrix(x)
  D = dmat$Delta
  if(missing(xlim)) {
    upper = max(dmat$left, dmat$right[f$right<Inf])
    xlim = range(0, upper * 1.05)
  }
  m = length(dmat$left)
  p = double(m)
  p[dmat$left %in% f$left & dmat$right %in% f$right] = f$p
  # g = colSums(w * D / (D %*% p)[,1]) - sum(w)
  P = (D %*% p)[,1]
  g = crossprod(w/P, D)[1,] - sum(w)
  plot(dmat$left, g, type="h", col=col2, xlab=xlab, ylab=ylab, xlim=xlim, ...)
  lines(xlim, c(0,0), lty=1)
  j = p > 0
  ms = sum(j)
  points(dmat$left[!j], rep(0,m-ms), pch=1, col=col2, cex=1)
  points(dmat$right[!j], rep(0, m-ms), pch=20, col=col2, cex=.8)
  segments(dmat$left[!j], rep(0, m-ms),
           pmin(dmat$right[!j], xlim[2]), rep(0, m-ms),
           col=col2, lwd=3)
  points(dmat$left[j], rep(0,ms), pch=1, col=col1, cex=1)
  points(dmat$right[j], rep(0, ms), pch=20, col=col1, cex=.8)
  segments(dmat$left[j], rep(0, ms), pmin(dmat$right[j], xlim[2]), rep(0, ms),
           col=col1, lwd=3)
} 

plotsurvidf = function(f, style=c("box","uniform","left","right","midpoint"),
    xlab="Time", ylab="Survival Probability", col="blue3", fill=NA,  
    add=FALSE, lty=1, lty.inf=2, xlim, ...) {
  style = match.arg(style)
  k = length(f$left)
  S = 1 - cumsum(f$p)
  upper = max(f$left, f$right[f$right != Inf])
  if(max(f$right) == Inf) point.inf = upper * 1.2
  else point.inf = upper
  if( missing(xlim) ) xlim = c(0, point.inf)
  m = length(f$p)
  if(is.na(fill)) {
    fill.hsv = drop(rgb2hsv(col2rgb(col))) * c(1, .3, 1)
    fill = hsv(fill.hsv[1], fill.hsv[2], fill.hsv[3], .3)
  }
  switch(style,
         box = {
           d = c(f$left[1], rep(f$right, rep(2,k)), f$right[k]) # right
           s = rep(c(1,S), rep(2,k+1))
           if(f$right[k] == Inf) d[2*k] = upper
           else d[2*k+2] = upper
           if( !add ) plot(d, s, type="n", col=col, xlim=xlim,
                           xlab=xlab, ylab=ylab, lty=lty, ...)
           if(style == "box") {
             Sc = c(1, S)
             for(i in 1:m) {
               if(f$right[i] - f$left[i] > 0) 
                 rect(f$left[i], Sc[i+1], f$right[i], Sc[i], border=col,
                      col=fill)
             }
           }
           lines(d, s, col=col, lty=lty, ...)
           lines(c(upper, point.inf), c(S[k-1],S[k-1]), col=col,
                 lty=lty.inf)
           if(f$right[k] != Inf) {       # left
             d = rep(c(f$left,f$right[k]), rep(2,k+1))
             s = c(1,rep(S, rep(2,k)),0)
           }
           else {
             d = rep(f$left, c(rep(2,k-1), 1))
             s = c(1,rep(S[-k], rep(2,k-1)))
           }
           add = TRUE
         }, 
         left = { d = rep(c(f$left,f$right[k]), rep(2,k+1))
                  s = c(1,rep(S, rep(2,k)),0)
                  d[2*k+2] = upper
                },
         right = { d = c(f$left[1], rep(f$right, rep(2,k)), f$right[k])
                   s = rep(c(1,S), rep(2,k+1))
                   if(f$right[k] == Inf) d[2*k] = upper
                   else d[2*k+2] = upper
                 },
         midpoint = { d1 = (f$left + f$right) / 2
                      d = c(f$left[1], rep(d1, rep(2,k)), f$right[k])
                      if(f$right[k] == Inf) d[2*k] = upper
                      else d[2*k+2] = upper
                      s = rep(c(1,S), rep(2,k+1))
                    },
         uniform = { d = c(rbind(f$left,f$right), rep(f$right[k],2))
                     if(f$right[k] == Inf) d[2*k] = upper
                     else d[2*k+2] = upper
                     s = c(1,rep(S, rep(2,k)),S[k])
                   }     )
  if( add ) lines(d, s, col=col, lty=lty,  ...)
  else plot(d, s, type="l", col=col, xlim=xlim, xlab=xlab, ylab=ylab,
            lty=lty, ...)
  abline(h=0, col="black")
  lines(c(0,f$left[1]), c(1,1), col=col)
  if(f$right[k] < Inf)
    lines(c(upper, point.inf), rep(0,2), col=col, lty=lty)
  else points(upper, S[k-1], col=col, pch=20)
}

