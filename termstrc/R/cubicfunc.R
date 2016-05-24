gi <- function(t,T,i,s){
  g <- rep(0,length(t))
  for(j in seq_along(t)){
    if(i==1){
    if(T[i]<=t[j]&t[j]<T[i+1]){
     g[j] <- (T[i])^2/6 + ((T[i])*(t[j]-T[i]))/2 + (t[j]-T[i])^2/2 - (t[j]-T[i])^3/(6*(T[i+1]-T[i]))
    }
    if(t[j]>=T[i+1]){
     g[j] <- (T[i+1])*((2*T[i+1]-T[i])/6 + (t[j]-T[i+1])/2)
    }   
  }
  if(i>1&i<length(T)){
    if(t[j]<T[i-1]){
     g[j] <- 0
    }
    if(T[i-1]<=t[j]&t[j]<T[i]){
     g[j] <- (t[j]-T[i-1])^3/(6*(T[i]-T[i-1]))
    }
    if(T[i]<=t[j]&t[j]<T[i+1]){
     g[j] <- (T[i]-T[i-1])^2/6 + ((T[i]-T[i-1])*(t[j]-T[i]))/2 + (t[j]-T[i])^2/2 - (t[j]-T[i])^3/(6*(T[i+1]-T[i]))
    }
    if(t[j]>=T[i+1]){
     g[j] <- (T[i+1]-T[i-1])*((2*T[i+1]-T[i]-T[i-1])/6 + (t[j]-T[i+1])/2)
    }
  }
   if(i==length(T)){
    if(t[j]<T[i-1]){
     g[j] <- 0
    }
    if(T[i-1]<=t[j]&t[j]<=T[i]){
     g[j] <- (t[j]-T[i-1])^3/(6*(T[i]-T[i-1]))
    }
  } 
  if(i==s){
    g[j] <- t[j]
  }
}
  g
}
