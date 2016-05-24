icasamp <-
  function(dname,query=c("rnd","pdf","kur"),nsamp=NULL,data=NULL){
    
    # initial checks
    dname <- dname[1]
    didx <- match(dname,letters[1:18])
    if(is.na(didx)){stop("Input 'dname' must be letter between 'a' and 'r'.")}
    query <- query[1]
    qidx <- match(query,c("rnd","pdf","kur"))
    if(is.na(qidx)){stop("Input 'query' must be 'rnd', 'pdf', or 'kur'.")}
    if(qidx==1L){
      if(is.null(nsamp[1])){stop("Input 'nsamp' must be provided.")}
      nsamp <- as.integer(nsamp[1])
      if(nsamp<=0){stop("Input 'nsamp' must be positive integer.")}
    } else if(qidx==2L){
      if(is.null(data[1])){stop("Input 'data' must be provided.")}
      data <- as.numeric(data)
    } 
    
    # sample data
    if(dname=="a"){
      # Student t with df=3
      if(query=="rnd"){
        return(rt(nsamp,3))
      } else if(query=="pdf"){
        return(dt(data,3))
      } else if(query=="kur"){
        return(Inf)
      } 
    } else if(dname=="b"){
      # double exponential
      if(query=="rnd"){
        return(sign(runif(nsamp)-0.5)*rexp(nsamp,rate=sqrt(2)))
      } else if(query=="pdf"){
        return(exp(-sqrt(2)*abs(data))/sqrt(2))
      } else if(query=="kur"){
        return(3)
      }  
    } else if(dname=="c"){
      # uniform
      if(query=="rnd"){
        return(runif(nsamp)*2*sqrt(3)-sqrt(3))
      } else if(query=="pdf"){
        return(1/2/sqrt(3)*(data<sqrt(3))*(data>-sqrt(3)))
      } else if(query=="kur"){
        return(-1.2)
      }  
    } else if(dname=="d"){
      # Student t with df=5
      if(query=="rnd"){
        return(rt(nsamp,5))
      } else if(query=="pdf"){
        return(dt(data,5))
      } else if(query=="kur"){
        return(6)
      }  
    } else if(dname=="e"){
      # Exponential
      if(query=="rnd"){
        return(-1+rexp(nsamp))
      } else if(query=="pdf"){
        return(dexp(data+1))
      } else if(query=="kur"){
        return(6)
      }  
    } else if(dname=="f"){
      # Mixture 2 Double Exponential
      prop <- rep(0.5,2)
      mus <- c(-1,1)
      covs <- rep(0.5,2)
      if(query=="rnd"){
        idx <- sample(1:2,nsamp,replace=TRUE)
        return(sign(runif(nsamp)-0.5)*rexp(nsamp,sqrt(2))*covs[idx]+mus[idx])
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:2){ myden <- myden + prop[i]/covs[i]*exp(-sqrt(2)*abs(data-mus[i])/covs[i])/sqrt(2) }
        return(myden)
      } else if(query=="kur"){
        mus <- mus*covs
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i])
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]+mus[i]^3)
          x4 <- x4 + prop[i]*(6*covs[i]^2+6*covs[i]*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="g"){
      # Mixture 2 Gaussian (symmetric & multimodal)
      prop <- rep(0.5,2)
      mus <- c(-0.5,0.5)
      covs <- rep(.15,2)
      if(query=="rnd"){
        idx <- sample(1:2,nsamp,replace=TRUE)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:2){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="h"){
      # Mixture 2 Gaussian (symmetric & transitional)
      prop <- rep(0.5,2)
      mus <- c(-0.5,0.5)
      covs <- rep(0.4,2)
      if(query=="rnd"){
        idx <- sample(1:2,nsamp,replace=TRUE)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:2){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="i"){
      # Mixture 2 Gaussian (symmetric & unimodal)
      prop <- rep(0.5,2)
      mus <- c(-0.5,0.5)
      covs <- rep(0.5,2)
      if(query=="rnd"){
        idx <- sample(1:2,nsamp,replace=TRUE)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:2){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="j"){
      # Mixture 2 Gaussian (nonsymmetric & multimodal)
      prop <- c(1,3)
      prop <- prop/sum(prop)
      mus <- c(-0.5,0.5)
      covs <- rep(0.15,2)
      if(query=="rnd"){
        idx <- sample(1:2,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:2){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="k"){
      # Mixture 2 Gaussian (nonsymmetric & transitional)
      prop <- c(1,2)
      prop <- prop/sum(prop)
      mus <- c(-0.7,0.5)
      covs <- rep(0.4,2)
      if(query=="rnd"){
        idx <- sample(1:2,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:2){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="l"){
      # Mixture 2 Gaussian (nonsymmetric & unimodal)
      prop <- c(1,2)
      prop <- prop/sum(prop)
      mus <- c(-0.7,0.5)
      covs <- rep(0.5,2)
      if(query=="rnd"){
        idx <- sample(1:2,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:2){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="m"){
      # Mixture 4 Gaussian (symmetric & multimodal)
      prop <- c(1,2,2,1)
      prop <- prop/sum(prop)
      mus <- c(-1,-.33,.33,1)
      covs <- rep(0.16,4)
      if(query=="rnd"){
        idx <- sample(1:4,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:4){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="n"){
      # Mixture 4 Gaussian (symmetric & transitional)
      prop <- c(1,2,2,1)
      prop <- prop/sum(prop)
      mus <- c(-1,-.2,.2,1)
      covs <- c(.2,.3,.3,.2)
      if(query=="rnd"){
        idx <- sample(1:4,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:4){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="o"){
      # Mixture 4 Gaussian (symmetric & unimodal)
      prop <- c(1,2,2,1)
      prop <- prop/sum(prop)
      mus <- c(-.7,-.2,.2,.7)
      covs <- c(.2,.3,.3,.2)
      if(query=="rnd"){
        idx <- sample(1:4,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:4){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="p"){
      # Mixture 4 Gaussian (nonsymmetric & multimodal)
      prop <- c(1,1,2,1)
      prop <- prop/sum(prop)
      mus <- c(-1,.3,-.3,1.1)
      covs <- c(.2,.2,.2,.2)
      if(query=="rnd"){
        idx <- sample(1:4,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:4){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="q"){
      # Mixture 4 Gaussian (nonsymmetric & transitional)
      prop <- c(1,3,2,.5)
      prop <- prop/sum(prop)
      mus <- c(-1,-.2,.3,1)
      covs <- c(.2,.3,.2,.2)
      if(query=="rnd"){
        idx <- sample(1:4,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:4){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } else if(dname=="r"){
      # Mixture 4 Gaussian (nonsymmetric & unimodal)
      prop <- c(1,2,2,1)
      prop <- prop/sum(prop)
      mus <- c(-.8,-.2,.2,.5)
      covs <- c(.22,.3,.3,.2)
      if(query=="rnd"){
        idx <- sample(1:4,nsamp,replace=TRUE,prob=prop)
        return(rnorm(nsamp,mus[idx],covs[idx]))
      } else if(query=="pdf"){
        myden <- 0
        for(i in 1:4){ myden <- myden + prop[i]*dnorm(data,mus[i],covs[i]) }
        return(myden)
      } else if(query=="kur"){
        mu <- x2 <- x4 <- x3 <- 0
        for(i in 1:length(prop)){
          mu <- mu + prop[i]*mus[i]
          x2 <- x2 + prop[i]*(mus[i]^2+covs[i]^2)
          x3 <- x3 + prop[i]*(3*mus[i]*covs[i]^2+mus[i]^3)
          x4 <- x4 + prop[i]*(3*covs[i]^4+6*covs[i]^2*mus[i]^2+mus[i]^4)
        }
        return((x4-4*mu*x3+6*mu^2*x2-3*mu^4)/(x2-mu^2)^2-3)
      }  
    } # end if(dname=="a")
    
  }