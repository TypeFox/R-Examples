legendre.symbol <- function(a, p){
  if(isprime(p)){
    K <- as.bigz(a)^(as.bigz((p-1)/2))
    modulus(K) <- p
    K[K==p-1] <- -1
    k <- integer(length(K))
    k[K==-1] <- -1
    k[K==1] <- 1
    A <- as.bigz(a)
    modulus(A) <- p
    k[A==0] <- 0
    k
  } else stop("p not prime")
}

# ... TODO:
#kronecker.symbol <- function(a, n){
#  if(n<0) au=1
#  else { if(a<0) au=-1 else au=1}
#  F <- factorize(n)
#  ext.legendre.symbol(a,p){
#    if(p>2)
#      legendre.symbol(a,p)
#    else{
#      if(a%%2==0) 0
#      else {
#        if(a%%8==1 | a%%8==7) 1
#        else -1
#      }
#    }
#  }
#  k <- 
#
#
#    K <- as.bigz(a^((p-1)/2))
#  modulus(K) <- p
#  K[K==p-1] <- -1
#  k <- integer(length(K))
#  k[K==-1] <- -1
#  k[K==1] <- 1
#  k[K==0] <- 0
#  k
#}

jacobsthal <- function(p){
  if(isprime(p)){
    Q=outer(1:p,1:p,function(i,j)legendre.symbol(i-j, p))
    return(Q)
  }else
  stop("p not prime")
}

paley <- function(n)
{
  if(isprime(n-1) && (n-1)%%4==3) {
                                        # Jacobsthal Matrix Q:
    Q=jacobsthal(n-1)
                                        # H =( 1 1^T) 
                                        #    ( 1 Q-I) 
    H=rbind(c(0,rep(1,n-1)),
            cbind(rep(-1,n-1), Q))+diag(1,n)
    return(H)
  }else stop("n-1 is not prime or n%%4 != 0") 
}

paley2 <- function(n){
  
  if( n%%4==0) {
    N=n/2
    if(isprime(N-1) && (N-1)%%4==1) {
      #Jacobsthal Matrix Q:
      Q=jacobsthal(N-1);
      # H =( 1 1^T)
      #    ( 1 Q-E)

      S0=matrix(c(1,-1,-1,-1),2,2,byrow=TRUE)
      S1=matrix(c(1,1,1,-1),2,2,byrow=TRUE)
      S=rbind(c(0,rep(1,N-1)),
              cbind(rep(1,N-1), Q))
      
      H=(S==0)%x%S0 + S%x%S1
      
      return(H);
    }else stop("n/2-1 is not prime or n/2%%4 != 2") 
  }else stop("n%%4 != 0") 
}



walsh <- function(H)
{
  N=dim(H)[1];
  if(N!=dim(H)[2])
    stop("H not square")
  H2=sylvester(1)%x%H

  return(H2)
}

sylvester <- function(k)
{       
  n=2^k;
  if(n==1)
    return(matrix(1,1,1));
  if(n==2)
    return(matrix(c(1, 1, 1, -1),2,2,byrow=TRUE));
  N=2
  H=sylvester(1)
  while(N<n){
    H=walsh(H)
    N=2*N
  }
  return(H)
}

hadamard.matrix <- function(n, verbose=FALSE, getnearest=TRUE)
{
  if(n==1)
    return(sylvester(1))
  if(n==2)
    return(sylvester(2))
  
  f <- factorize(n)
  if(all(f==2)){
    pow <- length(f)
    if(verbose)cat("Hadamard matrix of sylvester type\n")
    return(sylvester(pow))
  } else {
    if(n%%4==0){
      if(isprime(n-1)){
        if(verbose)cat("Hadamard matrix of paley I type\n")
        H=paley(n)
        return(H)
      }else{
        if(isprime(n/2-1) && (n/2-1)%%4==1) {
          if(verbose)cat("Hadamard matrix of paley II type\n")
          H=paley2(n)
          return(H)
        }else{
                                        # try n/2
          if(verbose)cat("try recursive call (walsh step)\n")
          H2=hadamard.matrix(n/2,getnearest=FALSE)
          if(!is.null(H2)){
            H=walsh(H2)
            return(H)
          }else{
            if(verbose) cat(paste("use stored matrix for size=",n,"\n"))
            hadamard.table(n)
          }
        }
      }
    }else{
      if(getnearest){
        warning(paste(n," != 4*k, get the smallest Hadamard matrix with at least dimension ",n))
        hadamard.matrix((n%/%4 +1)*4)
      }else{
        warning(paste(n," != 4*k"))
        return(NULL)
      }
    }
  }
}

is.hadamard <- function(H){
  n <- dim(H)[1]
  all((H%*%t(H)-n*diag(1,n))==0)
}


