
###############################################################
#          function for adding noise
###############################################################

setGeneric("noisy.sampling",
           function(x,var.adj=0,rng="rnorm",mean.adj=0,...,end.coef=0,n,order.adj=0,znoise)
           standardGeneric("noisy.sampling"))

setMethod("noisy.sampling",signature(x="yuima"),
          function(x,var.adj=0,rng="rnorm",mean.adj=0,...,end.coef=0,n,order.adj=0,znoise)
            noisy.sampling(x@data,var.adj=var.adj,rng=rng,mean.adj=mean.adj,...,
                           end.coef=end.coef,n=n,order.adj=order.adj,znoise=znoise))

setMethod("noisy.sampling",signature(x="yuima.data"),
          function(x,var.adj=0,rng="rnorm",mean.adj=0,...,
                   end.coef=0,n,order.adj=0,znoise){
            
            data <- get.zoo.data(x)
            
            d.size <- length(data)
            ser.times <- lapply(data,"time")
            total.samp <- unique(sort(unlist(ser.times)))
            samp.size <- length(total.samp)
            
            if(missing(n)) n <- sapply(data,"length")
            
            n <- as.vector(matrix(n,d.size,1))
            
            if(missing(znoise)) znoise <- as.list(double(d.size))
            
            result <- vector(d.size,mode="list")
            
            if(any(var.adj!=0)){ # exogenous noise is present
              rn <- n^(order.adj/2)*(matrix(do.call(rng,list(d.size*samp.size,...)),d.size,samp.size)-mean.adj)
              if(is.list(var.adj)){
                total.noise <- matrix(0,d.size,samp.size)
                for(d in 1:d.size){
                  idx <- is.element(time(var.adj[[d]]),total.samp)
                  total.noise[d,] <- sqrt(var.adj[[d]][idx])*rn[d,]
                }
              }else{
                if(d.size>1){
                  tmp <- svd(as.matrix(var.adj))
                  total.noise <- (tmp$u%*%diag(sqrt(tmp$d))%*%t(tmp$v))%*%rn
                }else{
                  total.noise <- sqrt(var.adj)*rn
                }
              }
              if(any(end.coef!=0)){
                if(is.list(end.coef)){
                  n.tmp <- n^((1-order.adj)/2)
                  for(d in 1:d.size){
                    noise <- subset(total.noise[d,],is.element(total.samp,ser.times[[d]]))
                    result[[d]] <- data[[d]]+noise+znoise[[d]]+
                      n.tmp*c(0,end.coef[[d]]*diff(as.numeric(data[[d]])))
                  }
                }else{
                  psi <- n^((1-order.adj)/2)*end.coef
                  for(d in 1:d.size){
                    noise <- subset(total.noise[d,],is.element(total.samp,ser.times[[d]]))
                    result[[d]] <- data[[d]]+noise+znoise[[d]]+
                      psi[d]*c(0,diff(as.numeric(data[[d]])))
                  }
                }
              }else{
                for(d in 1:d.size){
                  noise <- subset(total.noise[d,],is.element(total.samp,ser.times[[d]]))
                  result[[d]] <- data[[d]]+noise+znoise[[d]]
                }
              }
            }else{ # exogenous noise is absent
              if(any(end.coef!=0)){
                if(is.list(end.coef)){
                  n.tmp <- n^((1-order.adj)/2)
                  for(d in 1:d.size){
                    result[[d]] <- data[[d]]+znoise[[d]]+
                      n.tmp*c(0,end.coef[[d]]*diff(as.numeric(data[[d]])))
                  }
                }else{
                  end.coef <- n^((1-order.adj)/2)*end.coef
                  for(d in 1:d.size){
                    result[[d]] <- data[[d]]+znoise[[d]]+
                      end.coef[d]*c(0,diff(as.numeric(data[[d]])))
                  }
                }
              }else{
                for(d in 1:d.size){
                  result[[d]] <- data[[d]]+znoise[[d]]
                }
              }
            }
            
            return(setData(result))
          })
