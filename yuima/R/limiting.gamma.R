setGeneric("limiting.gamma",
           function(obj, theta, verbose=FALSE)
           standardGeneric("limiting.gamma")
           )
setMethod("limiting.gamma", "yuima",
          function(obj, theta, verbose=FALSE){
            limiting.gamma(obj@model, theta, verbose=verbose)
          })
setMethod("limiting.gamma", "yuima.model",
          function(obj, theta, verbose=FALSE){
            
            ## error check 1
            if(missing(obj)){
              stop("main object is missing.")
            }
            if(missing(theta)){
              stop("theta is missing.")
            }
            if(is.list(theta)==FALSE){
              stop("theta required list.\ntheta <- list(theta1, theta2)")
            }
            
            if(verbose){
              yuima.warn("Initializing... ")
            }
            
            r.size <- obj@noise.number
            d.size <- obj@equation.number
            state <- obj@state.variable
            THETA.1 <- obj@parameter@diffusion
            THETA.2 <- obj@parameter@drift
            mi.size <- c(length(THETA.1), length(THETA.2))
            
            ## error check 2
            if(d.size!=1){
              stop("this program is 1-dimention yuima limitation.")
            }
            if(mi.size[1]!=length(theta[[1]])){
              stop("the length of m1 and theta1 is different.")
            }
            if(mi.size[2]!=length(theta[[2]])){
              stop("the length of m2 and theta2 is different.")
            }
            
            
            Differentiation.vector <- function(myfunc, mystate, dim1, dim2){
              tmp <- vector(dim1*dim2, mode="expression")
              for(i in 1:dim1){
                for(j in 1:dim2){
                  tmp[(i-1)*dim2+j] <- parse(text=deparse(D(myfunc[i], mystate[j])))
                }
              }
              return(tmp)
            }
            
            Differentiation.scalar <- function(myfunc, mystate, dim){
              tmp <- vector(dim, mode="expression")
              for(i in 1:dim){
                tmp[i] <- parse(text=deparse(D(myfunc[i], mystate)))
              }
              return(tmp)
            }
            
            
            ## assign theta
            for(i in 1:mi.size[1]){
              assign(THETA.1[i], theta[[1]][i])
            }
            for(i in 1:mi.size[2]){
              assign(THETA.2[i], theta[[2]][i])
            }            
            if(verbose){
              yuima.warn("Done")
            }            
            ## p(x)            
            if(verbose){
              yuima.warn("get C ... ")
            }
            
            ## a part of p(x) function  
            tmp.y <- function(x.arg){
              assign(obj@state.variable, x.arg)
              a.eval <- eval(obj@drift)
              b.eval <- numeric(d.size)
              for(i in 1:d.size){
                b.eval[i] <- eval(obj@diffusion[[1]][i])
              }
              return(2*a.eval/as.double(t(b.eval)%*%b.eval))
            }
            
            ## a part of p(x) function for integrate 
            integrate.tmp.y <- function(x.arg){
              tmp <- numeric(length(x.arg))
              for(i in 1:length(x.arg)){
                tmp[i]<-tmp.y(x.arg[i])
              }
              return(tmp)
            }
            
            ## p(x) function without normalize const C 
            p0.x <- function(x.arg){
              assign(obj@state.variable, x.arg)
              b.eval <- numeric(d.size)
              for(i in 1:d.size){
                b.eval[i] <- eval(obj@diffusion[[1]][i])
              }
              return(exp(as.double(integrate(integrate.tmp.y, 0, x.arg)$value))/as.double(t(b.eval)%*%b.eval))
            }
            
            ## p0.x function for integrate
            integrate.p0.x <- function(x.arg){
              tmp <- numeric(length(x.arg))
              for(i in 1:length(x.arg)){
                tmp[i]<-p0.x(x.arg[i])
              }
              return(tmp)
            }
            
            ## normalize const
            C <- 1/integrate(integrate.p0.x, -Inf, Inf)$value            
            if(verbose){
              yuima.warn("Done")
            }
            
            ## gamma1            
            if(verbose){
              yuima.warn("Get gamma1 ... ")
            }
            counter.gamma1 <- 1
            
            ## gamma1 function
            get.gamma1 <- function(x.arg){
              assign(obj@state.variable, x.arg)
              b.eval <- numeric(d.size)
              for(i in 1:d.size){
                b.eval[i] <- eval(obj@diffusion[[1]][i])
              }
              Bi <- solve(t(b.eval)%*%b.eval)
              dTHETA.1.B <- array(0,c(d.size, d.size, mi.size[1]))
              tmp1.B <- array(0, c(d.size, d.size, mi.size[1]))
              tmp2.B <- numeric(mi.size[1])
              for(k in 1:mi.size[1]){
                dTHETA.1.b <- Differentiation.scalar(obj@diffusion[[1]], THETA.1[k], d.size)
                dTHETA.1.b.eval <- numeric(d.size)
                for(i in 1:d.size){
                  dTHETA.1.b.eval[i] <- eval(dTHETA.1.b[i])
                }
                tmp <- b.eval%*%t(dTHETA.1.b.eval)
                dTHETA.1.B[, , k] <- tmp+t(tmp)
                tmp1.B[, , k] <- dTHETA.1.B[, , k]%*%Bi
                tmp2.B[k] <- tmp1.B[, , k]  #sum(diag(tmp1.B[,,k])) 1-dimention limitation
              }
              tmp <- tmp2.B %*% t(tmp2.B) * p0.x(x.arg)
              return(tmp[((counter.gamma1-1)%%mi.size[1])+1, ((counter.gamma1-1)%/%mi.size[1])+1])
            }
            
            
            ## gamma1 function for integrate
            integrate.get.gamma1 <- function(x.arg){
              tmp <- numeric(length(x.arg))
              for(i in 1:length(x.arg)){
                tmp[i] <- get.gamma1(x.arg[i])
              }
              return(tmp*C/2)
            }
            
            ## calculating gamma1
            gamma1 <- matrix(0, mi.size[1], mi.size[1])
            for(i in 1:(mi.size[1]*mi.size[1])){
              gamma1[((i-1)%%mi.size[1])+1,((i-1)%/%mi.size[1])+1] <- integrate(integrate.get.gamma1, -Inf, Inf)$value
              counter.gamma1 <- counter.gamma1+1
            }            
            if(verbose){
              yuima.warn("Done")
            }
            
            ## gamma2            
            if(verbose){
              yuima.warn("Get gamma2 ... ")
            }
            counter.gamma2 <- 1
            
            ## gamma2 function
            get.gamma2 <- function(x.arg){
              assign(obj@state.variable, x.arg)
              if(mi.size[2]==1){
                dTHETA.2.a <- Differentiation.scalar(obj@drift, THETA.2, d.size)
              }else{
                dTHETA.2.a <- Differentiation.vector(obj@drift, THETA.2, d.size, mi.size[2])
              }
              dTHETA.2.a.eval <- numeric(mi.size[2])
              for(i in 1:mi.size[2]){
                dTHETA.2.a.eval[i] <- eval(dTHETA.2.a[i])
              }
              b.eval <- numeric(d.size)
              for(i in 1:d.size){
                b.eval[i] <- eval(obj@diffusion[[1]][i])
              }
              Bi <- solve(t(b.eval)%*%b.eval)
              tmp <- dTHETA.2.a.eval %*% Bi %*% t(dTHETA.2.a.eval) * p0.x(x.arg)
              return(tmp[((counter.gamma2-1)%%mi.size[2])+1, ((counter.gamma2-1)%/%mi.size[2])+1])
            }
            
            ## gamma2 function for intefrate
            integrate.get.gamma2 <- function(x.arg){
              tmp <- numeric(length(x.arg))
              for(i in 1:length(x.arg)){
                tmp[i]<-get.gamma2(x.arg[i])
              }
              return(tmp*C)
            }
            
            ## calculating gamma2
            gamma2 <- matrix(0, mi.size[2], mi.size[2])
            for(i in 1:(mi.size[2]*mi.size[2])){
              gamma2[((i-1)%%mi.size[2])+1, ((i-1)%/%mi.size[2])+1] <- integrate(integrate.get.gamma2, -Inf, Inf)$value
              counter.gamma2 <- counter.gamma2+1
            }            
            if(verbose){
              yuima.warn("Done")
            }
            
            ## make list for return 
            poi1 <- list(gamma1, gamma2)
            names(poi1) <- c("gamma1", "gamma2")
            poi2 <- append(gamma1, gamma2)            
            ret <- list(poi1, poi2)
            names(ret) <- c("list", "vec")
            
            return(ret)
          })
