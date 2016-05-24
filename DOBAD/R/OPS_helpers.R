## "_helpers" not exported in namespace


##########   the argument 'n' does NOTHING! kept around for backwards compatibility
########## eps.t is epsilon for time; eps.params is epsilon for parameters.
########## they behave differently, although perhaps they should behave the same?
########## as eps.params gets smaller the function gets faster but less safe.
########## 1e-12 is right around the boundary if time is say .1<t<1
########## As eps.t gets smaller the function gets slower but more accurate.
### i'm not sure what diff b/t eps.t and eps.params is. they should be 1 param.
process.prob.one.singlearg <- function(t,lambda,mu,nu,X0, Xt,
                                       ##eps.t=.Machine$double.eps^.65,
                                       eps.t=1e-10,
                                       eps.params=1e-10, n=-111){
  ##NOTE: for errorchecking whether res < 0: process.prob.one.fft always returns a result >=0.
  if (t < eps.t) { ##special cases first that cause over/underflow
    if (X0==Xt) return(1)
    else return(0)
  }
  else if ((X0*(lambda+mu)+nu)*t<eps.t ){
    if (X0==Xt) return(1)
    else return(0)
  }
  else if ((X0*lambda+nu)*t<eps.t ){##Now do several possible time-saving cases
    ##print("A")
    if (X0<Xt) return(0)
    else return(process.prob.one.fft(t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt))

  }
  else if (X0*mu*t < eps.t){
   ##print("B")
    if (X0>Xt) return(0)
    else return(process.prob.one.fft(t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt))
  }
  else if (lambda<eps.params || mu<eps.params){
##    print("C")
    ##the ".fft" is much slower, but can handle param-values around 0!
    return(process.prob.one.fft(t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt))
  }
  else if ( abs(lambda-mu) < eps.params){
    print("process.prob.one.singlearg: lambda equals mu (approx) so using fft; K&McG formulas apply but have not yet been adapted!")
    return(process.prob.one.fft(t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt))
  } 
  ##If didn't  return out above, then do the actual computations, and check whether result is >=0.
  res <- 0
  if (nu >eps.params){
##    print("G")
    beta <- nu/lambda;
    if (lambda<mu){res <- Pij.F(i=X0,j=Xt,t=t,L=lambda,m=mu,beta=beta)}
    else res <- Pij.E(i=X0,j=Xt,t=t,L=mu,m=lambda,beta=beta)
  }
  else {
##    print("H")
    nu <- 0; ## in case it were epsilon.
    if (X0==0) {
      if (Xt==0) res <- 1
      else res <- 0
    }
    else if (Xt==0) res <- process.prob.one.fft(t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt)
    else if (lambda<mu) res <- Pij.C(i=X0-1,j=Xt-1,t=t,L=lambda,m=mu,beta=2)
    else res <- Pij.D(i=X0-1,j=Xt-1,t=t,L=mu,m=lambda,beta=2)
  }
  if (res<0){
    print(paste("process.prob.one.singlearg error: Result is ", res, ", which is <0.  Returning 0.", sep=""))
    return(0)
  }
  else return(res)
}


##This is for L >= mu.
## note that i switch what L and mu are versus the K&G paper, so that
## L is for birth and mu for death always.
Pij.E <- function(i, j, t, L, m, beta){
  if (beta <= 0) {print("Pij.E: Error: beta <= 0 not allowed in E model");}
  if (L >= m) {print("Pij.E: Error: L >= m not allowed in E model");}
  a <- min(i,j); #BE CAREFUL TO USE a,b, NOT i,j TIL END
  b <- max(i,j); #Note the formula is for i<j. then use stationary dist'n to get Pji
  gamma <- L/m;
  elmt <- exp( (L-m)*t); #'emlt' is what the expression looks like ...
  sigma <- gamma * elmt;
  pi.b <- pi.j.E(L=L,m=m,beta=beta,j=b);  # "\pi_j, as in stationary dist'n (not "p_ij")
  Pab <- pi.b* (L/m)^(a+b) * factorial(b) / poch(beta,b) * exp((L-m)*beta*t) *
    getCommonTerms(a=a,b=b,elmt=elmt,gamma=gamma,sigma=sigma,beta=beta);
  if (i>j) { #Pab = Pji
    Pij <- pi.j.E(L=L,m=m,beta=beta,j=j)/pi.j.E(L=L,m=m,beta=beta,j=i) * Pab
  }
  else Pij <- Pab;
  Pij;
}

##recall: proces.prob.one(X0,Xt, nu=0) corresponds to PijD(X0-1,Xt-1,beta=2)
Pij.D<- function(i, j, t, L, m, beta){
  if (beta <= 0) {print("Pij.D: Error: beta <= 0 not allowed in D model?");}
  if (L >= m) {print("Pij.D: Error: L >= m not allowed in D model?");}
  a <- min(i,j); #BE CAREFUL TO USE a,b, NOT i,j TIL END
  b <- max(i,j); #Note the formula is for i<j. then use stationary dist'n to get Pji
  gamma <- L/m;
  elmt <- exp( (L-m)*t); #'emlt' is what the expression looks like ...
  sigma <- gamma * elmt;
  pi.b <- pi.j.D(L=L,m=m,beta=beta,j=b);  # "\pi_j, as in stationary dist'n (not "p_ij")

  Pab <- pi.b* poch(beta,a) / factorial(a) * (L/m)^(a+b) * exp((L-m)*t) *
    getCommonTerms(a=a,b=b,elmt=elmt,gamma=gamma,sigma=sigma,beta=beta);
  if (i>j) { #Pab = Pji
    Pij <- pi.j.D(L=L,m=m,beta=beta,j=j)/pi.j.D(L=L,m=m,beta=beta,j=i) * Pab
  }
  else Pij <- Pab;
  Pij;
}

Pij.C<- function(i, j, t, L, m, beta){
  if (beta <= 0) {print("Pij.C: Error: beta <= 0 not allowed in C model?");}
  if (L >= m) {print("Pij.C: Error: L >= m not allowed in C model?");}
  a <- min(i,j); #BE CAREFUL TO USE a,b, NOT i,j TIL END
  b <- max(i,j); #Note the formula is for i<j. then use stationary dist'n to get Pji
  gamma <- L/m;
  elmt <- exp( (L-m)*t); #'emlt' is what the expression looks like ...
  sigma <- gamma * elmt;
  pi.b <- pi.j.C(L=L,m=m,beta=beta,j=b);  # "\pi_j, as in stationary dist'n (not "p_ij")
  Pab <- pi.b*exp((L-m)*(beta-1)*t) * poch(beta,a) / factorial(a) *
    getCommonTerms(a=a,b=b,elmt=elmt,gamma=gamma,sigma=sigma,beta=beta);
  if (i>j) { #Pab = Pji
    Pij <- pi.j.C(L=L,m=m,beta=beta,j=j)/pi.j.C(L=L,m=m,beta=beta,j=i) * Pab
  }
  else Pij <- Pab;
  Pij;
}


  
# ".F" means this is the formula from the "F" process in Karlin-McGregor '58, "Linear growth, birth
#and death processes".
                                        # NOTE it has the modification that i flipped  an exponent
# -beta became beta (where nu = beta*L)
#NOTE this only works if BETA IS POSITIVE
Pij.F <- function(i, j, t, L, m, beta){
  if (beta <= 0) {print("Pij.F: Error: beta <= 0 not allowed in F model?");}
  if (L >= m) {print("Pij.F: Error: L >= m not allowed in F model?");} 
  a <- min(i,j); #BE CAREFUL TO USE a,b, NOT i,j TIL END
  b <- max(i,j); #Note the formula is for i<j. then use stationary dist'n to get Pji
  gamma <- L/m;
  elmt <- exp( (L-m)*t); #'emlt' is what the expression looks like ...
  sigma <- gamma * elmt;
  pi.b <- pi.j.F(L,m, beta, b);  # "\pi_j, as in stationary dist'n (not "p_ij")

  
#need simplify these , the calculations are excessive, many things cancel
# note that the beta function (maybe missing one extra factor) shows up; R might compute it faster
  # than my dummy code
  Pab <- 
    pi.b * ( factorial(b) / poch(beta, b)) * 
        getCommonTerms(a,b,elmt,gamma,sigma,beta);
                 
  Pij <- Pab; #unless:
  if (i > j){ # Pab = Pji
    Pij <- ( pi.j.F(L,m,beta,j=j) / pi.j.F(L,m,beta,j=i) ) * Pab;
  }
  Pij;
}



##########


## terms that all version of the probabilities share.
getCommonTerms <- function(a,b, elmt,gamma,sigma, beta){
  (1-gamma)^beta * #modification from formula in k&G
    (1-elmt)^a * (1-sigma)^-(a+beta) *
      getSumTerm(a,b,elmt,gamma,sigma,beta);
}
getSumTerm <- function(a,b,elmt,gamma,sigma,beta){
  k <- seq(0, a, 1);
  sumTerm <- choose(a, k) * (-1)^k * ((1-elmt/gamma)/(1-elmt))^k *
    ((1-elmt)/(1-sigma))^(b-k) *  poch(a+beta, b-k) / factorial(b-k)
  sumTerm <- sum(sumTerm);
}
                                        #Compute for the "F" process in Karlin McGregor '58
# stationarydistn(j), i.e. pi_j
pi.j.F <- function(L,m, beta, j){
  if (L >= m) {print("pi.j.F: Error: L >= m not allowed in F model")}
  (L/m)^j  *  poch(beta, j)/factorial(j);
}
pi.j.E <- function(L,m, beta, j){ ## redundant, just like F...
  if (L >= m) {print("pi.j.E: Error: L >= m not allowed in E model")}
  (m/L)^j  *  poch(beta, j)/factorial(j);
}
pi.j.D <- function(L,m, beta, j){ ## redundant, just like F...
  if (L >= m) {print("pi.j.D: Error: L >= m not allowed in D model?")}
  (m/L)^j  *  factorial(j) / poch(beta, j);
}
pi.j.C <- function(L,m, beta, j){ ## redundant, just like F...
  if (L >= m) {print("pi.j.C: Error: L >= m not allowed in C model?")}
  (L/m)^j  *  factorial(j) / poch(beta, j);
}

## compute
##  gamma(A+n)/gamma(A);
## but so that it doesn't overflow with large values unnecessarily
poch <- function(A, n){
  ##return(gamma(A+n)/gamma(A))
  poch1 <- function(oneN){prod(seq(from=A, by=1, length=oneN));}
  apply(matrix(n),1,poch1)
}

