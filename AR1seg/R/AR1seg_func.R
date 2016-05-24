library(Segmentor3IsBack)

AR1seg_func <-
function (y,Kmax=15,rho=TRUE){
  
  l = length(y)
 
  if (rho)  rho = median((diff(y,lag=2))^2)/median(diff(y)^2) - 1
 
  x = y[2:l] - rho*y[1:(l-1)]
  S = Segmentor(x, model = 2, Kmax = Kmax)
  breaks = S@breaks
  for (i in 1:Kmax){
    for (j in 1:i) breaks[i,j] = breaks[i,j]+1
  }
  rm(i,j)
  parameters = S@parameters
  
  PP = function(t){
    x=t
    l = length(x)
    i = 2
    while (l > 2 && i < l){
      if (x[i] == x[i-1] + 1 && x[i] != x[i+1] - 1){
        x = c(x[1:(i-1)],x[(i+1):l])
        l = l-1
      } else i=i+1
    }
    if (l>1 && x[l-1] == x[l]-1) x = x[1:(l-1)]
    x
  }
  
 
  PPbreaks = matrix(0,nrow = Kmax, ncol=Kmax, dimnames = dimnames(breaks))
  PPbreaks[1,]=breaks[1,]
  for (i in 2:Kmax){
    t = PP(breaks[i,1:(i-1)])
    PPbreaks[i,] = c(t,l,rep(0,Kmax-length(t)-1))
  }
  rm(i,t)
  
  fMa = function(t,mu) {
    M = c()
    t = c(0,t)
    for ( i in 2:length(t) ) {
      M = c(M, rep(mu[i-1], t[i] - t[i-1]))
    }
    M
  }
  
  sswg = function(br, param, series){
    sum((series - fMa( br, param))^2)
  }
  
  sswgseg = function(seg,seri){
    res = c()
    for (i in 1:(Kmax)){
      res = c(res, sswg(seg@breaks[i,1:i], seg@parameters[i,1:i], seri))
    }
    res
  }
  
  minushalflogB = function(t,u){
    t = t[t!=0]
    l = length(t)
    b = log(t[1]/u)
    if (l > 1){
      for (i in 2:l){
        b = b + log(t[i] - t[i-1])
      }
    }
    b = -b/2
  }
  
  
  ZS = function(seg,seri){
    u = length(seg@data)
    Kmax = seg@Kmax
    f = function(t) minushalflogB(t, u)
    wg = sswgseg(seg,seri)
    criterion = -(((u+1):(u-Kmax+2))/2)*log(wg) + lgamma(((u+1):(u-Kmax+2))/2) - (0:(Kmax-1))*log(u) + apply(seg@breaks, 1, f)
    selected = which.max(criterion)
    selected
  }
  
  

  selected = ZS(S,x)
  SelectedBreaks = breaks[selected,1:selected]
  PPSelectedBreaks = PPbreaks[selected,]
  PPSelectedBreaks=PPSelectedBreaks[PPSelectedBreaks!=0]
  PPselected=length(PPSelectedBreaks)

  vec1=c(1,PPSelectedBreaks[1:(PPselected-1)]+1)
  vec2=PPSelectedBreaks[1:(PPselected)]
  m=c()
 for (i in 1:PPselected){
  m[i]=mean(y[vec1[i]:vec2[i]])  
}

  list(data = y, rho = rho, decorrelated = x, breaks=breaks, PPbreaks = PPbreaks, selected=selected, SelectedBreaks = SelectedBreaks,  PPSelectedBreaks = PPSelectedBreaks, PPselected=PPselected,PPmean=m)
}
