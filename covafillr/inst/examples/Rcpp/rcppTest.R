library(covafillr)
library(inline)

cftest <- '
  cVector x1 = as<cVector>(x);
  cMatrix coord1 = as<cMatrix>(coord);
  cVector obs1 = as<cVector>(obs);
  int p1 = as<int>(p);
  cVector h1 = as<cVector>(h);
  covafill<double> cf(coord1,obs1,h1,p1);
  return wrap(cf(x1));'

fun <- cxxfunction(signature(x='numeric',
                             coord = 'matrix',
                             obs = 'numeric',
                             p = 'integer',
                             h = 'numeric'),
                   body = cftest,
                   plugin = 'covafillr'
                   )    



coord <- as.matrix(expand.grid(seq(-10,10,0.2),seq(-10,10,0.2)))
ftrue <- function(x)sum(x^3 - prod(x))
obs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))

fn0 <- function(x,h=c(1,1))fun(x,coord,obs,2L,h)
fn0(c(2.1,0))
system.time(fn0(c(2.1,0)))

fn <- Vectorize(function(x,y,d=1)fun(c(x,y),coord,obs,2L,c(1,1))[d],c('x','y'))
x <- y <- seq(-5,5,0.1)
system.time(ztrue <- outer(x,y,Vectorize(function(x,y)ftrue(c(x,y)))))
system.time(zfit <- outer(x,y,fn))

par(mfrow=c(1,2))
image(x,y,ztrue)
image(x,y,zfit)


sttest <- '
  cVector x1 = as<cVector>(x);
  cMatrix coord1 = as<cMatrix>(coord);
  cVector obs1 = as<cVector>(obs);
  int p1 = as<int>(p);
  cVector h1 = as<cVector>(h);
  covafill<double> cf(coord1,obs1,h1,p1);
  covatree<double> ct(100,&cf);
  return wrap(ct(x1));'

tfun <- cxxfunction(signature(x='numeric',
                             coord = 'matrix',
                             obs = 'numeric',
                             p = 'integer',
                             h = 'numeric'),
                   body = sttest,
                   plugin = 'covafillr'
                   )    


tfn <- Vectorize(function(x,y,d=1)tfun(c(x,y),coord,obs,2L,c(1,1))[d],c('x','y'))
system.time(tzfit <- outer(x,y,fn))

par(mfrow=c(1,3))
image(x,y,ztrue)
image(x,y,zfit)
image(x,y,tzfit)
