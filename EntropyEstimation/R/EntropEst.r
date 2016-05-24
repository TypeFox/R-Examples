KL.Plugin <- function(x,y){
   if(length(x) != length(y)){
      cat("Arrays must have same sizes.")
   }
  else .C("KlPlugin", as.integer(x), as.integer(y), as.integer(length(x)), as.double(0.0) )[[4]]
}

KL.z<- function(x,y){
  if(length(x) != length(y)){
    cat("Arrays must have the same sizes.")
  }
  else .C("KlSharp", as.integer(x), as.integer(y), as.integer(length(x)), as.double(0.0) )[[4]]
}

SymKL.z <- function(x,y){
  if(length(x) != length(y)){
    cat("Arrays must have the same sizes.")
  }
  else  .5*(.C("KlSharp", as.integer(x), as.integer(y), as.integer(length(x)), as.double(0.0) )[[4]] + .C("KlSharp", as.integer(y), as.integer(x), as.integer(length(x)), as.double(0.0) )[[4]])


  #  .C("SymSharp", as.integer(x), as.integer(y), as.integer(length(x)), as.double(0.0) )[[4]]
  
}

SymKL.Plugin <- function(x,y){
  if(length(x) != length(y)){
    cat("Arrays must have the same sizes.")
  }
else    .5*(.C("KlPlugin", as.integer(x), as.integer(y), as.integer(length(x)), as.double(0.0) )[[4]] + .C("KlPlugin", as.integer(y), as.integer(x), as.integer(length(x)), as.double(0.0) )[[4]])
  #  .C("SymPlugin", as.integer(x), as.integer(y), as.integer(length(x)), as.double(0.0) )[[4]]
  
}

KL.sd<- function(x,y){
  if(length(x) != length(y)){
    cat("Arrays must have the same sizes.")
  }
  else .C("KlSd", as.integer(x), as.integer(y), as.integer(length(x)), as.double(0.0) )[[4]]
  
}

SymKL.sd<- function(x,y){
  if(length(x) != length(y)){
    cat("Arrays must have the same sizes.")
  }
  else .C("SymSd", as.integer(x), as.integer(y), as.integer(length(x)), as.double(0.0) )[[4]]
}

Entropy.z<- function(x){
        .C("EntropySharp", as.integer(x), as.integer(length(x)), as.double(0.0) )[[3]]
    
}

GenSimp.z<- function(x,r){
    if(r<1){
        cat("r must be a positive integer.")
    }
    else{
        if(r>=sum(x)){
            cat("r must be strictly less than sum(x).")
        }
        else .C("GenSimpSharp", as.integer(x), as.integer(length(x)), as.integer(r), as.double(0.0) )[[4]]
    }
}

GenSimp.sd<- function(x,r){
    if(r<1){
        cat("r must be a positive integer.")
    }
    else{
        if(r>=sum(x)){
            cat("r must be strictly less than sum(x).")
        }
        else .C("GenSimpSd", as.integer(x), as.integer(length(x)), as.integer(r), as.double(0.0) )[[4]]
    }
}

RenyiEq.z<- function(x,r){
    if(r <= 0){
        cat("r must be greater than zero")
    }
    else{
        if(r == 1){
            r
        }
        else
        {
            .C("RenyiEqEntropySharp", as.integer(x), as.integer(length(x)), as.double(r), as.double(0.0) )[[4]]
        }
    }
    
}

Renyi.z<- function(x,r){
    if(r <= 0){
        cat("r must be greater than zero")
    }
    else{
        if(r==1){
            Entropy.z(x)
        }
        else
        {
            log(RenyiEq.z(x,r))/(1-r)
        }
    }
    
}

Tsallis.z<- function(x,r){
    if(r <= 0){
        cat("r must be greater than zero")
    }
    else{
        if(r==1){
            Entropy.z(x)
        }
        else
        {
            (RenyiEq.z(x,r)-1)/(1-r)
        }
    }
}

Hill.z<- function(x,r){
    if(r <= 0){
        cat("r must be greater than zero")
    }
    else{
        if(r==1){
            exp(Entropy.z(x))
        }
        else
        {
            (RenyiEq.z(x,r))^(1/(1-r))
        }
    }
}

Entropy.sd<- function(x){
    .C("EntropySd", as.integer(x), as.integer(length(x)), as.double(0.0) )[[3]]
}

RenyiEq.sd<- function(x,r){
    if(r <= 0){
        cat("r must be greater than zero")
    }
    else{
        if(r==1){
            0
        }
        else
        {
            .C("RenyiEqSd", as.integer(x), as.integer(length(x)), r, as.double(0.0) )[[4]]
        }
    }
}

MI.z<- function(x){
    Entropy.z(rowSums(x)) +Entropy.z(colSums(x)) -Entropy.z(as.vector(x))
}

MI.sd<- function(y){
    x=y
    
    if(x[length(x[,1]),length(x[1,])]==0){
        loc = which(x>0,arr.ind = T)[1,]
        x = rbind(x[-loc[1],], x[loc[1],])
        x = cbind(x[,-loc[2]], x[,loc[2]])
    }
    
    r = rowSums(x)
    c = colSums(x)
    g = matrix(rep(0,length(x)),length(r),length(c),byrow=TRUE)
    
    for(i in 1:length(r)){
        for(j in 1:length(c)){
            if(x[i,j]==0) g[i,j]=0
            else{
                if(i< length(r)){
                    if(j< length(c)){
                        g[i,j] = log(r[length(r)]*c[length(c)]*x[i,j]) - log(r[i]*c[j]*x[length(r),length(c)])
                    }
                    else{
                        if(j== length(c)){
                            g[i,j] = log(r[length(r)]*x[i,j]) - log(r[i]*x[length(r),length(c)])
                        }
                    }#else
                }#if i
                else{
                    if(j< length(c)){
                        g[i,j] = log(c[length(c)]*x[i,j]) - log(c[j]*x[length(r),length(c)])     
                    }
                }
            }#else
        }#for j
    }#for i
    x=as.vector(x)
    g = as.vector(g)
    xx=x[x!=0]
    gg = g[x!=0]
    .C("MISd", as.integer(as.vector(xx)), as.integer(length(xx)), as.double(as.vector(gg)), as.double(0.0))[[4]]
}


Renyi.sd<- function(x,r){
    if(r <= 0){
        cat("r must be greater than zero")
    }
    else{
        if(r==1){
            cat("r should not equal 1, use Entropy.sd instead.")
        }
        else{
            RenyiEq.sd(x,r)/(abs(1-r)*RenyiEq.z(x,r))
        }
    }
}

Tsallis.sd<- function(x,r){
    if(r <= 0){
        cat("r must be greater than zero")
    }
    else{
        if(r==1){
            cat("r should not equal 1, use Entropy.sd instead.")
        }
        else{
            RenyiEq.sd(x,r)/(abs(1-r))
        }
    }
}

Hill.sd<- function(x,r){
    if(r <= 0){
        cat("r must be greater than zero")
    }
    else{
        if(r==1){
            Entropy.sd(x)*exp(Entropy.z(x))
        }
        else{
            RenyiEq.sd(x,r)*(RenyiEq.z(x,r))^(r/(1-r))/abs(1-r)
        }
    }
}