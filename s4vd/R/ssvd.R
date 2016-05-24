ssvd <-
function(X, threu = 1, threv = 1, gamu = 0, gamv =0,  u0, v0,
                   merr = 1e-4, niter = 100)
{
    if(missing(u0) || missing(v0))
    {
##
## Use irlba efficient partial svd for large problems
      if(min(dim(X)) < 500) S <- svd(X[],nu=1,nv=1)
      else S <- irlba(X,nu=1,nv=1)
      u0 <- S$u
      v0 <- S$v
    }
    i <- NULL #to pass CRAN checks without NOTE
    n    <- dim(X)[1]
    d    <- dim(X)[2]
    stop <- FALSE
    ud   <- 1
    vd   <- 1
    iter <- 0
    SST  <- sum(X^2)
    while (ud > merr || vd > merr) {
        iter <- iter+1
        cat("iter: ",iter,"\n")
        cat("v: ",length(which(v0!=0)),"\n")
        cat("u: ",length(which(u0!=0)),"\n")
        # Updating v
##        z =  t(X)%*% u0;
## avoid explicitly forming the transpose of X:
        z <- t(t(u0) %*% X)[,drop=FALSE]
        winv <- abs(z)^gamv
        sigsq <- abs(SST - sum(z^2))/(n*d-d)

        tv <- sort(c(0, abs(z*winv)))
        rv <- sum(tv>0)
        Bv <- rep(1,d+1)*Inf
    
#        for (i in 1:rv){
#            lvc  =  tv[d+1-i]
#            temp1 = which(winv!=0)
#            
#            temp2 = s4vd:::thresh(z[temp1], type = threv, delta = lvc/winv[temp1])
#            vc = rep(0,d)
#            vc[temp1] = temp2
##            Bv[i] = sum((X - u0%*%t(vc))^2)/sigsq + i*log(n*d)
#            Bv[i] = sum((z - vc)^2)/sigsq + i*log(n*d)
# Note: The modified Bv[i] differs from Bv[i] by a constant, which
# doesn't change the BIC optimization problem. This observation also
# applies to Bu[i] below...
#        }
        temp1 <- which(winv!=0)
        Bv[1:rv] <- foreach(i=1:rv,.combine=c,.multicombine=TRUE) %dopar% {
          lvc  <-  tv[d+1-i]
          delta <- lvc/winv[temp1]
          temp2 <- sign(z[temp1]) * (abs(z[temp1]) >= delta) * (abs(z) - delta)
          vc <- rep(0,d)
          vc[temp1] <- temp2
          sum((z - vc)^2)/sigsq + i*log(n*d)  ## Note
        }

        Iv <- min(which(Bv==min(Bv)))
        temp <- sort(c(0, abs(z*winv)))
        lv <- temp[d+1-Iv]
    
        temp2 <- thresh(z[temp1],type = threv, delta = lv/winv[temp1])
        v1 <- rep(0,d)
        v1[temp1] <- temp2
        v1 <- v1/sqrt(sum(v1^2)) #v_new
    
        cat("v1", length(which(v1!=0)) ,"\n" )
        # Updating u
        z <- (X %*% cbind(v1))[,drop=FALSE]      ## Note
        winu <- abs(z)^gamu
        sigsq <- abs(SST - sum(z^2))/(n*d-n)

        tu <- sort(c(0,abs(z*winu)))
        ru <- sum(tu>0)
        Bu <- rep(1,n+1)*Inf

#        for (i in 2:ru){
#            luc  =  tu[n+1-i]
#            temp1 = which(winu!=0)
#            temp2 = s4vd:::thresh(z[temp1], type = threu, delta = luc/winu[temp1])
#            uc = rep(0,n)
#            uc[temp1] = temp2
##            Bu[i] = sum((X - uc%*%t(v1))^2)/sigsq + i*log(n*d)
#            Bu[i] = sum((z - uc)^2)/sigsq + i*log(n*d)
#        }
        temp1 <- which(winu!=0)
        Bu[1:ru] <- foreach(i=1:ru,.combine=c) %dopar% {
          luc  <-  tu[n+1-i]
          delta <- luc/winu[temp1]
          temp2 <- sign(z[temp1])*(abs(z[temp1]) >= delta)*(abs(z[temp1])-delta)
          uc <- rep(0,n)
          uc[temp1] <- temp2
          sum((z - uc)^2)/sigsq + i*log(n*d)
        }

        Iu <- min(which(Bu==min(Bu)))
        temp <- sort(c(0, abs(z*winu)))
        lu <- temp[n+1-Iu]
        temp2 <- thresh(z[temp1],delta = lu/winu[temp1])
        u1 <- rep(0,n)
        u1[temp1] <- temp2
        u1 <- u1/sqrt(sum(u1^2))

        ud <- sqrt(sum((u0-u1)^2))
        vd <- sqrt(sum((v0-v1)^2))

        if (iter > niter){
          print("Fail to converge! Increase the niter!")
          stop = TRUE
          break
        }
        
        u0 <- u1
        v0 <- v1
    }
    return(list(u = u1, v = v1, iter = iter,stop=stop))
}
