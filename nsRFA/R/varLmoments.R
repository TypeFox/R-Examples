varLmoments <- function (x, matrix=TRUE) {
 y <- sort(x)
 n <- length(y)
 nn <- rep(n-1, n)
 pp <- seq(0, n-1)
 p1 <- pp/nn
 p2 <- p1*(pp-1)/(nn-1)
 p3 <- p2*(pp-2)/(nn-2)
 b0 <- sum(y)/n
 b1 <- sum(p1*y)/n
 b2 <- sum(p2*y)/n
 b3 <- sum(p3*y)/n
 l1 <- b0
 l2 <- 2*b1 - b0
 l3 <- 6*b2 - 6*b1 + b0
 l4 <- 20*b3 - 30*b2 + 12*b1 - b0
 tau <- l2/l1
 tau3 <- l3/l2
 tau4 <- l4/l2
 Y1 <- y %*% t(rep(1,n))
 Y <- Y1*t(Y1)
 Q <- seq(1,n) %*% t(rep(1,n))
 P <- t(Q)
 varb0 <- b0^2 - 1/n/(n-1) *2* sum(Y*lower.tri(Y))

 W11 <- 1/n/(n-1)/(n-2)/(n-3) *2*((P-1)*(Q-3))
 V11 <- W11*lower.tri(W11)*Y
 varb1 <- b1^2 - sum(V11)

 W10 <- 1/n/(n-1)/(n-2)*((Q-2)+(P-1))
 V10 <- W10*lower.tri(W10)*Y
 covb0b1 <- b0*b1-sum(V10)

 W20 <- 1/n/(n-1)/(n-2)/(n-3) * ((Q-2)*(Q-3)+(P-1)*(P-2))
 V20 <- W20*lower.tri(W20)*Y
 covb0b2 <- b0*b2-sum(V20)

 W21 <- 1/n/(n-1)/(n-2)/(n-3)/(n-4)*((P-1)*(Q-3)*(Q-4)+(P-1)*(P-2)*(Q-4))
 V21 <- W21*lower.tri(W21)*Y
 covb1b2 <- b1*b2-sum(V21)

 W22 <- 1/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)*2*((P-1)*(P-2)*(Q-4)*(Q-5))
 V22 <- W22*lower.tri(W22)*Y
 varb2 <- b2*b2-sum(V22)

 W30 <- 1/n/(n-1)/(n-2)/(n-3)/(n-4) * ((Q-2)*(Q-3)*(Q-4)+(P-1)*(P-2)*(P-3))
 V30 <- W30*lower.tri(W30)*Y
 covb0b3 <- b0*b3-sum(V30)

 W31 <- 1/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)*((P-1)*(Q-3)*(Q-4)*(Q-5)+(P-1)*(P-2)*(P-3)*(Q-5))
 V31 <- W31*lower.tri(W31)*Y
 covb1b3 <- b1*b3-sum(V31)

 W32 <- 1/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)/(n-6)*((P-1)*(P-2)*(P-3)*(Q-5)*(Q-6)+(P-1)*(P-2)*(Q-4)*(Q-5)*(Q-6))
 V32 <- W32*lower.tri(W32)*Y
 covb2b3 <- b2*b3-sum(V32)

 W33 <- 1/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)/(n-6)/(n-7)*2*((P-1)*(P-2)*(P-3)*(Q-5)*(Q-6)*(Q-7))
 V33 <- W33*lower.tri(W33)*Y
 varb3 <- b3*b3-sum(V33)
 if (n > 7) {
  T <- matrix(c(varb0,covb0b1,covb0b2,covb0b3,
                covb0b1,varb1,covb1b2,covb1b3,
                covb0b2,covb1b2,varb2,covb2b3,
                covb0b3,covb1b3,covb2b3,varb3),nrow=4,ncol=4,byrow=TRUE)
  C <- matrix(c(1,0,0,0,
                -1,2,0,0,
                1,-6,6,0,
                -1,12,-30,20),nrow=4,ncol=4,byrow=TRUE)
  varL <- C %*% T %*% t(C)
  dimnames(varL) <- list(c("l1","l2","l3","l4"),c("l1","l2","l3","l4"))
  if (matrix==FALSE) {
   varl1 <- varL[1,1]
   varl2 <- varL[2,2]
   varl3 <- varL[3,3]
   varl4 <- varL[4,4]
   covl1l2 <- varL[1,2]
   covl2l3 <- varL[2,3]
   covl2l4 <- varL[2,4]
   varlcv <- tau^2*(varl1/l1^2+varl2/l2^2-2*covl1l2/l1/l2)
   varlca <- tau3^2*(varl2/l2^2+varl3/l3^2-2*covl2l3/l2/l3)
   varlkur <- tau4^2*(varl2/l2^2+varl4/l4^2-2*covl2l4/l2/l4)
   varL <- c(varl1,varl2,varl3,varl4,varlcv,varlca,varlkur)
   names(varL) <- c("var.l1","var.l2","var.l3","var.l4","var.lcv","var.lca","var.lkur")
  }
 }
 else if (n > 5) {
  T <- matrix(c(varb0,covb0b1,covb0b2,
                covb0b1,varb1,covb1b2,
                covb0b2,covb1b2,varb2),nrow=3,ncol=3,byrow=TRUE)
  C <- matrix(c(1,0,0,
                -1,2,0,
                1,-6,6),nrow=3,ncol=3,byrow=TRUE)
  varL <- C %*% T %*% t(C)
  dimnames(varL) <- list(c("l1","l2","l3"),c("l1","l2","l3"))
  if (matrix==FALSE) {
   varl1 <- varL[1,1]
   varl2 <- varL[2,2]
   varl3 <- varL[3,3]
   covl1l2 <- varL[1,2]
   covl2l3 <- varL[2,3]
   varlcv <- tau^2*(varl1/l1^2+varl2/l2^2-2*covl1l2/l1/l2)
   varlca <- tau3^2*(varl2/l2^2+varl3/l3^2-2*covl2l3/l2/l3)
   varL <- c(varl1,varl2,varl3,varlcv,varlca)
   names(varL) <- c("var.l1","var.l2","var.l3","var.lcv","var.lca")
  }
 }
 else if (n > 3) {
  T <- matrix(c(varb0,covb0b1,
                covb0b1,varb1),nrow=2,ncol=2,byrow=TRUE)
  C <- matrix(c(1,0,
                -1,2),nrow=2,ncol=2,byrow=TRUE)
  varL <- C %*% T %*% t(C)
  dimnames(varL) <- list(c("l1","l2"),c("l1","l2"))
  if (matrix==FALSE) {
   varl1 <- varL[1,1]
   varl2 <- varL[2,2]
   covl1l2 <- varL[1,2]
   varlcv <- tau^2*(varl1/l1^2+varl2/l2^2-2*covl1l2/l1/l2)
   varL <- c(varl1,varl2,varlcv)
   names(varL) <- c("var.l1","var.l2","var.lcv")
  }
 }
 else {
  varL <- varb0
  names(varL) <- "var.l1"
 }
 return(varL)
}


# ----------------------------------------------------------------------------------------------- #

varLCV <- function (x) {
 y <- sort(x)
 n <- length(y)
 nn <- rep(n-1, n)
 pp <- seq(0, n-1)
 p1 <- pp/nn
 b0 <- sum(y)/n
 b1 <- sum(p1*y)/n
 l1 <- b0
 l2 <- 2*b1 - b0
 tau <- l2/l1
 Y1 <- y %*% t(rep(1,n))
 Y <- Y1*t(Y1)
 Q <- seq(1,n) %*% t(rep(1,n))
 P <- t(Q)
 varb0 <- b0^2 - 1/n/(n-1) *2* sum(Y*lower.tri(Y))

 W11 <- 1/n/(n-1)/(n-2)/(n-3) *2*((P-1)*(Q-3))
 V11 <- W11*lower.tri(W11)*Y
 varb1 <- b1^2 - sum(V11)

 W10 <- 1/n/(n-1)/(n-2)*((Q-2)+(P-1))
 V10 <- W10*lower.tri(W10)*Y
 covb0b1 <- b0*b1-sum(V10)

 T <- matrix(c(varb0,covb0b1,
               covb0b1,varb1),nrow=2,ncol=2,byrow=TRUE)
 C <- matrix(c(1,0,
               -1,2),nrow=2,ncol=2,byrow=TRUE)
 varL <- C %*% T %*% t(C)
 varl1 <- varL[1,1]
 varl2 <- varL[2,2]
 covl1l2 <- varL[1,2]
 varlcv <- tau^2*(varl1/l1^2+varl2/l2^2-2*covl1l2/l1/l2)
 names(varlcv) <- "var.lcv"
 return(varlcv)
}


# ----------------------------------------------------------------------------------------------- #

varLCA <- function (x) {
 y <- sort(x)
 n <- length(y)
 nn <- rep(n-1, n)
 pp <- seq(0, n-1)
 p1 <- pp/nn
 p2 <- p1*(pp-1)/(nn-1)
 b0 <- sum(y)/n
 b1 <- sum(p1*y)/n
 b2 <- sum(p2*y)/n
 l1 <- b0
 l2 <- 2*b1 - b0
 l3 <- 6*b2 - 6*b1 + b0
 tau <- l2/l1
 tau3 <- l3/l2
 Y1 <- y %*% t(rep(1,n))
 Y <- Y1*t(Y1)
 Q <- seq(1,n) %*% t(rep(1,n))
 P <- t(Q)
 varb0 <- b0^2 - 1/n/(n-1) *2* sum(Y*lower.tri(Y))

 W11 <- 1/n/(n-1)/(n-2)/(n-3) *2*((P-1)*(Q-3))
 V11 <- W11*lower.tri(W11)*Y
 varb1 <- b1^2 - sum(V11)

 W10 <- 1/n/(n-1)/(n-2)*((Q-2)+(P-1))
 V10 <- W10*lower.tri(W10)*Y
 covb0b1 <- b0*b1-sum(V10)

 W20 <- 1/n/(n-1)/(n-2)/(n-3) * ((Q-2)*(Q-3)+(P-1)*(P-2))
 V20 <- W20*lower.tri(W20)*Y
 covb0b2 <- b0*b2-sum(V20)

 W21 <- 1/n/(n-1)/(n-2)/(n-3)/(n-4)*((P-1)*(Q-3)*(Q-4)+(P-1)*(P-2)*(Q-4))
 V21 <- W21*lower.tri(W21)*Y
 covb1b2 <- b1*b2-sum(V21)

 W22 <- 1/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)*2*((P-1)*(P-2)*(Q-4)*(Q-5))
 V22 <- W22*lower.tri(W22)*Y
 varb2 <- b2*b2-sum(V22)

 T <- matrix(c(varb0,covb0b1,covb0b2,
               covb0b1,varb1,covb1b2,
               covb0b2,covb1b2,varb2),nrow=3,ncol=3,byrow=TRUE)
 C <- matrix(c(1,0,0,
               -1,2,0,
               1,-6,6),nrow=3,ncol=3,byrow=TRUE)
 varL <- C %*% T %*% t(C)
 varl2 <- varL[2,2]
 varl3 <- varL[3,3]
 covl2l3 <- varL[2,3]
 varlca <- tau3^2*(varl2/l2^2+varl3/l3^2-2*covl2l3/l2/l3)
 names(varlca) <- "var.lca"
 return(varlca)
}


# ----------------------------------------------------------------------------------------------- #

varLkur <- function (x) {
 L <- Lmoments(x)
 l2 <- L[2]
 l4 <- L[5]*L[2]
 tau4 <- L[5]
 varL <- varLmoments(x)
 varl2 <- varL[2,2]
 varl4 <- varL[4,4]
 covl2l4 <- varL[2,4]
 varlkur <- tau4^2*(varl2/l2^2+varl4/l4^2-2*covl2l4/l2/l4)
 names(varlkur) <- "var.lkur"
 return(varlkur)
}

