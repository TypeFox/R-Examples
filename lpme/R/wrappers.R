.rlaplace=function (use.n, location = 0, scale = 1) 
{
  location <- rep(location, length.out = use.n)
  scale <- rep(scale, length.out = use.n)
  rrrr <- runif(use.n)
  location - sign(rrrr - 0.5) * scale * (log(2) + ifelse(rrrr < 0.5, log(rrrr), log1p(-rrrr)))
}

.dMISE = function(W, sig, error="laplace"){
  n = length(W);
  FK  = function(t) { ifelse( (t<=1 & t>=-1), (1-t^2)^8, 0) }
  if(error=="laplace"){
    hROT = (5*sig^4/n)^(1/9);
    hlist= seq(hROT*0.05, hROT*1.9, 0.01);
    FU  = function(t) {1/(1+(sig*t)^2/2)};
  }else{
    hROT = sqrt(2)*sig/sqrt(log(n));
    hlist= seq(hROT*0.5, hROT*1.9, 0.01);
    FU = function(t) { exp(-(sig*t)^2/2) };
  }
  dd  = density(W);
  x   = dd$x;
  y   = dd$y; 
  nx  = length(x);
  mu2 = 16;
  fW2 = c(0, (y[-c(1,2)]-2*y[-c(1,nx)]+y[-c(nx-1,nx)])/(x[2]-x[1])^2, 0)
  Rf2 = sum(fW2[-1]^2*(x[-1]-x[-nx]));
  mise = function(h){
    integral=integrate( function(x) (FK(x)^2)/min(FU(x/h)^2, 1e30), -1, 1 )$value
    1/(2*pi*n*h)*integral + h^4/4*Rf2*mu2
  }
  hlist[which.min(sapply(hlist, mise))];
}

.dmise = function(W, FK, FU, mu2, hlist){
  n = length(W);
  dd=density(W);
  x = dd$x;
  y = dd$y; nx = length(x);
  fW2 = c(0, (y[-c(1,2)]-2*y[-c(1,nx)]+y[-c(nx-1,nx)])/(x[2]-x[1])^2, 0)
  Rf2 = sum(fW2[-1]^2*(x[-1]-x[-nx]));
  mise = function(h) {
    integral=integrate( function(x) (FK(x)^2)/min(FU(x/h)^2, 1e300), -1, 1 )$value
    1/(2*pi*n*h)*integral + h^4/4*Rf2*mu2
  }
  hlist[which.min(sapply(hlist, mise))];
}
