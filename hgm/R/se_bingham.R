### Pfaffian for Bingham distribution
hgm.se.dG.fun.Bingham = function(th, G, fn.params=list(d=rep(1,length(th)+1), logarithm=FALSE)){
  d = fn.params$d
  logarithm = fn.params$logarithm
  thn = length(th)  # length(G) should be equal to (thn + 1) = p
  Gn = length(G)
  G.res = G[1] - sum(G[1+(1:thn)])
  if(logarithm) G.res = 1 - sum(G[1+(1:thn)])
  dG = array(0,c(thn, Gn))
  dG[1:thn,1] = G[1+(1:thn)]
  for(i in 1:thn){
    dG[i, 1+i] = G[1+i]
    for(j in 1:thn){
      if(j != i){
        dG[i,1+j] = (d[j] * G[1+i] - d[i] * G[1+j]) / 2 / (th[i] - th[j])
        dG[i,1+i] = dG[i,1+i] - dG[i,1+j]
      }
    }
    dG[i,1+i] = dG[i,1+i] - (d[thn+1] * G[1+i] - d[i] * G.res) / 2 / th[i]
  }
  if(logarithm){
    arg = 1:thn
    dG[arg, 1+arg] = dG[arg, 1+arg] - outer(G[1+arg], G[1+arg])
  }
  dG
}

### Normalising constant and its derivatives (by power series expansion)
# degenerate case refers to Kume and Wood (2007)
hgm.se.C.Bingham.power = function(th.all, v=rep(0,length(th.all)), d=rep(1,length(th.all)), Ctol=1e-6, thtol=1e-10, Nmax=1000){
  th = th.all
  p = length(th)  # not (p-1)
  thsum = sum(abs(th))
  N0 = max(ceiling(thsum), 1)
  for(N in N0:Nmax){
    logep = N*log(thsum) - lfactorial(N) + log((N+1)/(N+1-thsum))
    if(logep < log(Ctol)) break
  }
  if(N > Nmax) stop("too large value: th.all.")
  f = function(k){
    kn = length(k)
    if(kn == p){
      w = which(k > 0)
      a = prod( th[w]^k[w] )
      b1 = sum( - lfactorial(k) + lgamma(k + v + (d/2)) )
      b2 = - lgamma(sum(k + v + (d/2))) + lgamma(sum(d)/2) - sum(lgamma(d/2))
      return( a * exp(b1+b2) )
    }else{ # recursive part
      knxt = kn + 1
      if(abs(th[knxt]) < thtol){  # speed-up
        return( f(c(k,0)) )
      }else{
        a = 0
        imax = N - sum(k)
        for(i in 0:imax) a = a + f(c(k,i))
        return(a)
      }
    }
  }
  f(c())
}

### Initical value for Pfaffian system (by power series expansion)
hgm.se.G.Bingham.power = function(th, d=rep(1,length(th)+1), ...){
  p = length(th) + 1
  th.all = c(th,0)
  C= hgm.se.C.Bingham.power(th.all, d=d, ...)
  dC = numeric(p-1)
  for(i in 1:(p-1)){
    e = rep(0,p); e[i] = 1
    dC[i] = hgm.se.C.Bingham.power(th.all, v=e, d=d, ...)
  }
  c(C,dC)
}

### Initial value (wrapper)
hgm.se.G.Bingham = function(th, d=rep(1,length(th)+1), method="power", withvol=FALSE, ...){
  dsum = sum(d)
  if(withvol) v = 2 * pi^(dsum/2)/gamma(dsum/2)
  else v = 1
  if(method == "power"){
  	return( v * hgm.se.G.Bingham.power(th, d=d, ...) )
  }
  # TODO: prepare Monte Carlo method "MC"
  stop("'method' not found")
}

### hg method for computing G(th)
hgm.se.hgm.Bingham = function(th, d=rep(1,length(th)+1), logarithm=FALSE, ini.method="power", times=NULL, withvol=FALSE, ...){
  th0 = th / sum(abs(th)*d[1:length(th)]) / 2
  G0 = hgm.se.G.Bingham(th0, d=d, method=ini.method, withvol=withvol, ...)
  # TODO: prepare initial values to reduce computing time
  p = length(th) + 1
  if(logarithm){
    G0[2:p] = G0[2:p] / G0[1]
    G0[1] = log(G0[1])
  }
  dG.fun = hgm.se.dG.fun.Bingham
  hgm.se.hgm(th0, G0, th, dG.fun, times=times, fn.params=list(d=d, logarithm=logarithm))
}
