#############################################################################
# Random KNN Properties                                           
# File:   protert.R                                               
# Author: Shengqiao Li                                            
# Date:   June 24, 2008 (initial)                                 
# Dependency: class, drep, Hmisc                                  
# Change Log:                                                     
#           March 02, 2010 - using gmp to compute eta             
#           2013-08-06 -- using gmp 0.5-5 (chooseZ, pow.bigz and crossprod)         
############################################################################

#Overloading problem %*% and prod in this package
#explicitly call %*%.bigq, %*%.bigz, prod.bigz to solve

rbyb<- function(p, m, eta)
{
  log(1-eta^(1/p)) / log(1-m/p)
}

rbyp<- function(p, m, eta)
{
  (log(-log(eta)) - log(p)) / log(1-m/p)

}
rbyv<- function(p, m, nu)
{
  nu*p/m
}

rbyz<- function(p, m)
{
 
  #if(!require(gmp)) stop("gmp package is required. Otherwise use other method!")
  
  Cpm<- chooseZ(p,m);
  j<- 1:p;

  res<- crossprod((-1)^(j+1) * Cpm, cbind(chooseZ(p, j)/(Cpm-chooseZ(p-j, m))))

  round(as.double(res))
}

rbyz.sim<- function(p, m, nsim=1000)
{

  r<- function(p, m)
  {
    x<- NULL;
    i<- 0;
    repeat{
     if( length(setdiff(1:p, x)>0)){
        i<- i+1;
        x<- union(x, sample(p, m))
      }
      else break;
  }

    return(i)
  }

  Z<- numeric(nsim)
  for(i in 1:nsim) Z[i]<- r(p, m)

  return(round(mean(Z)))

}

rbyz.geo<- function(p, m=floor(sqrt(p)), rmax=p)
{
  #approximate
  #max of iid geometric r.v.
  #P(rf<=r) = P(If=0) = 1-q^r = 1 - (1-m/p)^r
  #P(rz<=r) = [1 - (1-m/p)^r]^p

  #P(rz>r) = 1-[1 - (1-m/p)^r]^p
  #E[rz] = sum(P(rz>r)) over r=0:inf

  r<- 0:rmax;
  rmax+1 - sum((1-(1-m/p)^r)^p) #don't know why square is required

}
rbylambda<- function(p, m, lambda=1)
{

  log(lambda/p)/log(1-m/p)

}
#end nonvisible functions


#eta<- function(p, m, r)
#{

#  j<- 0:p;
#  sum((-1)^(p-j)*choose(p, j)*(choose(j, m)/choose(p,m))^r)
#
#}
 


eta<- function(p, m, r, method=c("binomial", "poisson", "exact"))
{
  #coverage probability
  method <- match.arg(method);
  switch(method,
    exact = {
#        if(!require(gmp)) stop("gmp package is required. Otherwise use other method!")
        j<- 0:p;
        res<- crossprod((-1)^(p-j), cbind(chooseZ(p, j)*pow.bigz(chooseZ(j, m), r))) / pow.bigz(chooseZ(p,m), r)
        as.double(res);
    },
    binomial = {
        #binominal approximation
        (1-(1-m/p)^r)^p     
    },
    poisson = {
       #poisson approximation
       exp(- p*(1-m/p)^r)
    }
  )
}

r<- function(p, m=floor(sqrt(p)), eta=.99, nu=20, rmax=p, nsim=1000, lambda=0.01,
          method=c("binomial", "poisson", "nu", "geo.independent", "geo.sim", "geo.dependent", "lambda"))
{
  # 7 methods to choose number of KNNs
  #eta -- coverage probability
  #nu -- mean mutiplicity of a feature
  #lambda -- mean number of silient features
  #rmax  -- number of series terms for independent geometric approximation
  #nsim -- number of simulations for geometric simulation
  
  method <- match.arg(method);
  r<- switch(method,
        binomial =  rbyb(p, m, eta),
        poisson =   rbyp(p, m, eta),
        nu =        rbyv(p, m, nu),
        geo.independent =  rbyz.geo(p, m, rmax),
        geo.sim = rbyz.sim(p, m, nsim),
        geo.dependent =  rbyz(p, m),
        lambda = rbylambda(p, m, lambda)  
           
  );
  ceiling(r);
}

lambda<- function(p, m, r)
{
  #mean nubmer of silient features
  p*(1-m/p)^r
}

 