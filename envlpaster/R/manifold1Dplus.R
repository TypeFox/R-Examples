#####################################################
###                                               ###
###  Functions for the 1D algortihm               ###
###                                               ###
#####################################################


###################################################
#                projection                       #
###################################################
projection <- function(a){
  d <- dim(a)[2]
  if(sum(t(a)%*%a)==0){
    return(0)
  }
  pa <- a%*%matpower(t(a)%*%a,-1)%*%t(a)
  return(pa)
}

##################################################
#                  matrix power                  #
##################################################
matpower <- function(a,alpha){
  small <- 0.000001
  p1<-nrow(a)
  eva<-eigen(a)$values
  eve<-as.matrix(eigen(a)$vectors)
  eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
  index<-(1:p1)[eva>small]
  evai<-eva
  evai[index]<-(eva[index])^(alpha)
  foo <- NULL
  if(length(evai) == 1) foo <- diag(evai, nrow = 1)
  else foo <- diag(evai)
  ai<-as.matrix(eve)%*%foo%*%t(as.matrix(eve))
  return(ai)
}

##################################################
#         1D objective function                  #
##################################################
get1Dobj <- function(w,A,B){
  small <- 0.000001
  p <- dim(A)[1]
  foo <- eigen((A+B), symmetric = TRUE)
  if(p == 1) B.int <- foo$vec %*% 1/foo$val %*% t(foo$vec)
  else B.int <- foo$vec %*% diag(1/foo$val) %*% t(foo$vec)
  Fw <- log(t(w)%*%A%*%w + small) + log(t(w)%*%B.int%*%w + small) - 2*log(t(w)%*%w)
  return(Fw)
}

##################################################
#    get initial value for 1D algorithm          #
##################################################
get1Dini <- function(A,B){
  p <- dim(A)[1]
  vecs <- cbind(eigen(A, symmetric = TRUE)$vectors,
    eigen(A+B, symmetric = TRUE)$vectors)
  idx <- order(apply(vecs,2,get1Dobj,A,B))[1]
  w <- vecs[,idx]
  return(w)
}

##################################################
#         1D objective function gradient         #
##################################################
get1Dderiv <- function(w,A,B){
  p <- dim(A)[1]
  foo <- eigen((A + B), symmetric = TRUE)
  if(p == 1) B.int <- foo$vec %*% 1/foo$val %*% t(foo$vec)
  else B.int <- foo$vec %*% diag(1/foo$val) %*% t(foo$vec)
  dF <- c(2/(t(w)%*%A%*%w))*A%*%w + c(2/(t(w)%*%B.int%*%w))*B.int%*%w - c(4/(t(w)%*%w))*w
  return(dF)
}

##################################################
# 1D manifold algorithm for u-dim envelope       #
# the algorithm based on M and inv(M+U)          #
# and Polak-Ribiere conjugate gradient (PRCG)    #
##################################################
manifold1Dplus <- function(M,U,u){
  p <- dim(M)[1]
  Mnew <- M
  Unew <- U
  G <- matrix(0,p,u)
  G0 <- diag(1,p)
  for(i in 1:u){  # used to be A = solve(Mnew + Unew), B = Mnew
    ans <- optim(get1Dini(Mnew,Unew),get1Dobj,get1Dderiv,
                 A=Mnew,B=Unew,method="CG",
                 control=list(maxit=500,type=2))
    w <- c(ans$par)
    gk <- c(1/sqrt(sum(w^2)))*w
    if(p == 1) G[,i] <- G0 * gk
    else G[,i] <- G0%*%gk
    G0 <- qr.Q(qr(G[,1:i]),complete=T)
    G0 <- G0[,(i+1):p]
    Mnew <- t(G0)%*%M%*%G0
    Unew <- t(G0)%*%U%*%G0
  }
  return(G)
}

