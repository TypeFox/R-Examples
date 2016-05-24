# functions for positive exchangeable MVN and derivative with respect
# to a_k, b_k or rho

exchmvn<- function(lb, ub, rh, mu=0, scale=1, eps = 1.e-06)
{ #eps=1.e-6 is recommended value because of bound on Romberg iterations
  #if(!is.loaded(symbol.C("exchmvn"))) dyn.load("./r_exchmvn.so")
  if(rh<0 || rh>=1) stop("0<=rh<1")
  m <- length(ub)
  if(m!=length(lb)) stop("lengths of lb and ub must be the same")
  tem <- scale
  a <- (lb - mu)/tem
  b <- (ub - mu)/tem
  out <- .C("r_exchmvn",
    as.integer(m), as.double(a), as.double(b), as.double(rh),
    as.double(eps), pr = as.double(eps))
  out$pr
}

# lb=lower limit, ub=upper limit, rh>0 is correlation,
# eps is tolerance for Romberg integration,
# k=argument of lb/ub for deriv, ksign=1 for upper, -1 for lower
exchmvn.deriv.margin<- function(lb, ub, rh, k, ksign, eps = 1.e-06)
{ #eps=1.e-6 is recommended value because of bound on Romberg iterations
  #if(!is.loaded(symbol.C("emvnd"))) dyn.load("./r_exchmvn.so")
  if(rh<0 || rh>=1) stop("0<=rh<1")
  m <- length(lb)
  if(m!=length(ub)) stop("lengths of lb and ub must be the same")
  if(k<1) k<-1
  if(k>m) k<-m
  if(abs(ksign)!=1) ksign<-1
  out <- .C("r_emvnd",
    as.integer(m), as.double(lb), as.double(ub), as.double(rh),
    as.integer(k), as.integer(ksign),
    as.double(eps), deriv = as.double(eps))
  out$deriv
}


# deriv of exchmvn with respect to rho
exchmvn.deriv.rho<- function(lb, ub, rh, eps = 1.e-06)
{ #eps=1.e-6 is recommended value because of bound on Romberg iterations
  #if(!is.loaded(symbol.C("emvndrh"))) dyn.load("./r_exchmvn.so")
  if(rh<0 || rh>=1) stop("0<=rh<1")
  m <- length(lb)
  if(m!=length(ub)) stop("lengths of lb and ub must be the same")
  out <- .C("r_emvndrh",
    as.integer(m), as.double(lb), as.double(ub), as.double(rh),
    as.double(eps), deriv = as.double(eps))
  out$deriv
}

# test case 
# m<-5
# a<-rep(-1,m)
# b<-rep(2,m)
# rh<-.6
# print(exchmvn(a,b,rh))

# answer should be 0.5437533059

# names of functions made longer and self-explanatory

#b[m]=1.5
#print(emvnd(a,b,rh,1,-1))
#print(emvnd(a,b,rh,1,1))
#print(emvnd(a,b,rh,m,-1))
#print(emvnd(a,b,rh,m,1))
#print(emvndrh(a,b,rh))

# answers should be 
#integ:     dera1=-0.084819, derb1=0.025738
#integ:     dera1=-0.085593, derb1=0.093674
#integ.   : derrh=0.417630
