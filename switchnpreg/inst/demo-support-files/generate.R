######################################
## Functions to generate the data  ###
######################################

true.functions <- function(X, S1, S2){
    N <- length(X)
    XDif <- outer(X, X, FUN='-')  
    K1 <- (1 / (sqrt(2*pi) * S1)) * exp(-1/2 * ((XDif)^2 / (S1^2)))
    K2 <- (1 / (sqrt(2*pi) * S2)) * exp(-1/2 * ((XDif)^2 / (S2^2))) 
    f <- cbind(mvrnorm(1, rep(0, N), Sigma=K1),mvrnorm(1, rep(0, N), Sigma=K2))
    return(f)
}


z.indep <- function(n.tmp,alpha.tmp){
            number.states <- length(alpha.tmp)
            z <- sample(x=1:number.states,size=n.tmp,
                        replace=TRUE,prob=alpha.tmp)
            return(z)
            }

data.y <- function(f.tmp,n.tmp,s2.tmp,z.tmp){
          Y <- vector(mode='numeric',length=n.tmp)
          for (j in 1:dim(f.tmp)[2]) {
              Y[z.tmp==j] <- f.tmp[,j][z.tmp==j] + rnorm(sum(z.tmp==j),mean=0,sd=sqrt(s2.tmp[j]))
              }
          return(Y)
          }
          
make.z.y <-  function(f,sigma2,alpha,z.make.function){
              N <- dim(f)[1]
              Z <- z.make.function(n.tmp=N,alpha.tmp=alpha)
              Y <- data.y(f.tmp=f,n.tmp=N,s2.tmp=sigma2,z.tmp=Z)
              return <- list(z=Z,y=Y)
              }

z.markov <- function(n.tmp,alpha.tmp){
      ## Create a Markov chain object with no data (NULL)
      MarkovChain <- HiddenMarkov::mchain(NULL, Pi=alpha.tmp$A, delta=alpha.tmp$PI,nonstat=FALSE)
      ## Simulate some data
      z <- simulate(MarkovChain, nsim=n.tmp)$mc
      return(z)
          }
          
data.z.markov.sim <- function(n,NSIM,alpha){
MarkovChain <- HiddenMarkov::mchain(NULL, Pi=alpha$A, delta=alpha$PI,nonstat=FALSE)
Zmatrix <- matrix(NA,nrow=n,ncol=NSIM)
  for (j in 1:NSIM){
    Zmatrix[,j] <-  simulate(MarkovChain, nsim=n)$mc
    }
    
return(Zmatrix)
}
          

          
        
