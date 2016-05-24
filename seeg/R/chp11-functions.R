semimarkov.rain <- function(P, Hk, Ha, tsim){

# arguments: markov matrix P 
# Ha rates of holding times, Hk order of holding times
# simulation time

# initial
x <- array(); tau <- array()
i=1;t=0
# start with state at random
y <- runif(1,0,1)
if(y > 0.5) x[i] <- 1 else x[i] <- 0 

while(t< tsim){
 i<- i+1
# apply embedded markov
  y <- runif(1,0,1)
  if(x[i-1]==0) {
   if(y > P[1,1]) {x[i] <- 1; tau[i] <- rgamma(1,shape=Hk[1,1],rate=Ha[1,1])}
   else           {x[i] <- 0; tau[i] <- rgamma(1,shape=Hk[2,1],rate=Ha[2,1])}
  }
  else {
   if(y > P[1,2]) {x[i] <- 1; tau[i] <- rgamma(1,shape=Hk[1,2],rate=Ha[1,2])}
   else           {x[i] <- 0; tau[i] <- rgamma(1,shape=Hk[2,2],rate=Ha[2,2])}
  }
 t <- t+ tau[i]
}

return(list(x=round(x,2),tau=round(tau,2),t=t))
}

semimarkov <- function(P, Hk, Ha, tsim, xinit){

# arguments: markov matrix P 
# Ha rates of holding times, Hk order of holding times
# simulation time

# dim state
nstate <- dim(P)[1]
# initial
x <- array(); tau <- array(); t <- array(); z <- array()
i=1;t[i]=0; x[1] <- xinit
while(t[i]< tsim){
 i<- i+1
 y <- runif(1,0,1)
 z[1] <- 0; for(j in 2:(nstate+1)) z[j] <- z[j-1] + P[j-1,x[i-1]]
 for(j in 2:(nstate+1)){
  if(y>z[j-1]&& y<=z[j]) {x[i] <- j-1;tau[i]<- rgamma(1,shape=Hk[j-1,x[i-1]],rate=Ha[j-1,x[i-1]])}
}
 t[i] <- t[i-1]+ tau[i]
} # end while
return(list(x=round(x,0),tau=round(tau,2),t=t))
}


