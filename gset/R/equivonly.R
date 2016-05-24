equivonly <-
function(l, u, sigma, n1, n2, t.vec, type1, gamma=-4,crange=c(-10,10), 
         plot=TRUE, ll=3, ul=6, n.sim=1e4,seed=NULL) 
{
  HSD <- function(t, error, gamma) error*(1-exp(-gamma*t))/(1-exp(-gamma))
  I.vec <- HSD(t.vec,type1, gamma)  
  K <- length(t.vec) 
  
  nn1 <- ceiling(n1 * t.vec); nn1[K] <- n1
  nn2 <- ceiling(n2 * t.vec); nn2[K] <- n2
  
  if(!is.null(seed)) set.seed(seed)
  
  simdata <- data.frame(matrix(0,n.sim*K,3))
  colnames(simdata)<- c("stage", "t.L", "t.U")
  x1 <- matrix(rnorm(n1*n.sim, 0, sigma), n.sim, n1)
  x2 <- matrix(rnorm(n2*n.sim, l, sigma), n.sim, n2)             
  
  for(k in 1:K){
    x1.stage   <- x1[,1:nn1[k]]   
    mean1.temp <- apply(x1.stage, 1, mean)
    v1.temp   <- apply(x1.stage, 1, var)
    
    x2.stage   <- x2[,1:nn2[k]]   
    mean2.temp <- apply(x2.stage, 1, mean)
    v2.temp   <- apply(x2.stage, 1, var)  
    
    v.pool <- (v1.temp*(nn1[k]-1)+ v2.temp*(nn2[k]-1))/(nn1[k]+ nn2[k]-2)
    se.pool <- sqrt(v.pool*(1/nn1[k]+1/nn2[k]))
    diff <- mean2.temp -  mean1.temp
    
    t.L <-  (diff - l) / se.pool
    t.U <-  (diff - u) / se.pool 
    simdata[(k-1)*n.sim + 1:n.sim,] <- cbind(k, t.L, t.U)		
  }
  
  ct2.L <- rep(-99, K)
  ct2.U <- rep(-99, K)
  
  f <- function(x){
    temp1 <- simdata[(simdata$stage == 1),]$t.L > x
    temp2 <- simdata[(simdata$stage == 1),]$t.U < -x
    result<- sum(temp1 & temp2) / n.sim - I.vec[1]
    return(result)
  }
  ct2.L[1] <- uniroot(f, crange)$root
  ct2.U[1] <- -ct2.L[1]      
  
  for (k in 2:K) {   
    stagealpha <- I.vec[k] - I.vec[k-1]    
    
    temp1 <- (simdata[(simdata$stage == k - 1),]$t.L > ct2.L[k - 1])
    temp2 <- (simdata[(simdata$stage == k - 1),]$t.U < ct2.U[k - 1])
    reject.prev <- temp1 & temp2 
    if(any(is.na(reject.prev)))  reject.prev[is.na(reject.prev)] <- TRUE
    if(any(reject.prev)){
      simdata[(simdata$stage == k),][reject.prev,]$t.L <- NA
      simdata[(simdata$stage == k),][reject.prev,]$t.U <- NA
    }
    
    if(all(is.na(simdata$t.L)) == FALSE){  	
      simdata.k <- simdata[(simdata$stage == k),]
      left.k <- simdata.k[!is.na(simdata.k$t.L), ] 
      if(nrow(left.k) < stagealpha*n.sim) 
        stop(paste("the study has stopped at stage ", as.character(k)))
      else{ 
        f <- function(x) {
          temp1 <- left.k$t.L > x 
          temp2 <- left.k$t.U < -x
          result <- sum( temp1 & temp2) /n.sim - stagealpha
          return(result)
        }
        ct2.L[k] <- uniroot(f, crange, tol=0.001)$root
        ct2.U[k] <- -ct2.L[k]
      }
    }
  }
  
  result <- list(typeI= I.vec, equivL=ct2.L, equivU=ct2.U)
  if(plot){
    if(K<=3) par(mfrow=c(1,K), mar=c(4,2,1,1), ask=T)
    else if(K==4)  par(mfrow=c(2,2), mar=c(4,2,1,1), ask=T)
    else if(K>=5 & K<=6)  par(mfrow=c(2,3), mar=c(4,2,1,1), ask=T)
    figureE(result, K, ll, ul)
  }
  return(result)
}
