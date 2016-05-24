"logLik" <-
function(y, z=matrix(1, ncol=2,nrow=length(y)),off=rep(0,length(y)), theta, beta,density.name="negbin",link="log", modify=unity, gamma=NULL, mod.Z=NULL,zero.start=FALSE){
  Lik <- matrix(1,ncol=1,nrow=1)
  chopped <- chop(y,z,off)
  y <- chopped$y
  z <- chopped$z
  off <- chopped$off
  len.y <- length(y)
  if( any( c(is.na(unlist( z[-len.y,])) , is.na(off[-len.y])))){  stop("Missing values in explanatory variables: clean your data or do multiple imputation")
                                                                }

  mu <- linkfn( off+z%*%beta,family=link)
  
	if( zero.start){
          y <- c(0,y)
        }
  if( any( diff( y[!is.na(y)])<0)){ warning("some of your increment values were set to zero")}
  dy <- pmax(diff( y[!is.na(y)]),0)
  k <- 1
  temp.y <- y[1]
  if( length(y)<2) stop("each patient must have at least 2 observations")
  for( i  in 2:(length(y))){
    if( deparse(substitute(modify))=="unity"){
      temp.mu <- mu[i-1]
    }
    else{
#      print( length(temp.y))
#      print( (i-1):(i-2+length(temp.y)))
#      print( dim(mod.Z))
        
      
    temp.mu <-linkfn( invlinkfn( mu[i-1], family=link)+eval( call(modify, x=temp.y, y=gamma ,mod.Z=mod.Z[(i-1),, drop=FALSE])), family=link)
  } 
	if( !is.na( y[i] ) ){
      if( !is.na( y[i-1]) ){
        Lik <- Lik%*%densityfn( dy[k],density.name,  size=theta, mu=temp.mu)
      }
      if( is.na( y[i-1])){
        Lik <- Lik %*% densityfn( dy[k]:0, density.name, size=theta,mu=temp.mu)
        
      }
      k <- k+1
      temp.y <- y[i]
    }
    if( is.na(y[i])){
      if( !is.na(y[i-1])){
        Lik <- Lik%*%densityfn( 0:dy[k], density.name, size=theta, mu=temp.mu)
        temp.y <- y[i-1]+0:dy[k]
      }
      if( is.na( y[i-1])){
        S <- dy[k]+1
        P <- matrix(0, nrow=S, ncol=S)
        for( row in 1:S){
          for( col in row:S){
            P[row,col] <- densityfn( col-row, density.name, size=theta, mu=temp.mu[row])
          }
        }
      
	Lik <- Lik%*%P
      }
    }
  }
 
  log(Lik)
}

