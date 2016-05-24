sumeventsf <-
function(i,X,R,Ftime,theta,dl,n,p,v,eventtimes){
    xj <- matrix(X[i:n,],nrow=(n-i+1))
    rj <- matrix(R[i:n,],nrow=(n-i+1))
    Fi <- Ftime[eventtimes==i,]
    rr <-exp(rj%*%dl+xj %*% theta %*% Fi) #relative risk
    hi <- 1/sum(rr)  #baseline hazard
    xmean <- hi*t(xj)%*%rr
    rmean <- hi*t(rj)%*%rr  # first dev wrt d
    xdif <- xj -matrix(rep(xmean,n-i+1),ncol=p, byrow=T)
    rdif <- rj -matrix(rep(rmean,n-i+1),ncol=v, byrow=T)
    XX <- t(xdif)%*% (xdif*matrix(rep(rr,p),ncol=p,byrow=F))
    FF <- Fi %*% t(Fi)
    RR <- t(rdif) %*%(rdif*matrix(rep(rr,v),ncol=v,byrow=F)) #second der wrt dxd
    RB <- t(xdif) %*% (rdif * matrix(rep(rr, v), ncol = v, byrow = FALSE)) #second der wrt dx
    j1 <- cbind(kronecker(FF,XX),kronecker(FF[,1],RB))
    j2 <- rbind(j1,matrix(c(kronecker(FF[,1],RB),RR),nrow=v, byrow = FALSE))
    list (hi,                      # hazard hi
          -xmean%*% t(Fi),            # part of term i of score
          -rmean,
          hi*j2)}
