vif.lm <- function(x) {

 # VIF = variance inflation factors
 # INPUT
 # x = output of 'lm'

 xind <- x$model[2:x$rank]
 k <- dim(xind)[2]
 nomi <- names(xind)
 VIFs <- rep(NA,k)
 for(i in 1:k) {
  f <- formula(paste(nomi[i],"~ 1 +",paste(nomi[-i],collapse=" + "),collapse=""))
  VIFs[i] <- 1/(1 - summary(lm(f, xind))$r.squared)
 }
 names(VIFs) <- nomi
 return(VIFs)
}


# ------------------------------------------------------------------ #

R2.lm <- function(x) {

 # INPUT
 # x = output of 'lm'

 R2 <- c(summary(x)$r.squared,summary(x)$adj.r.squared)
 names(R2) <- c("R2","adjR2")

 return(R2)
}


# ----------------------------------------------------------------- #

prt.lm <- function(x) {

 # INPUT
 # x = output of 'lm'

 prt <- summary(x)$coefficients[-1,"Pr(>|t|)"]

 return(prt)
}


# ------------------------------------------------------------------ #

RMSE.lm <- function(x) {

 # INPUT
 # x = output of 'lm'

 res <- x$residuals

 RMSE <- sqrt(sum((res)^2)/length(res))

 return(RMSE)
}


# ------------------------------------------------------------------ #

MAE.lm <- function(x) {

 # INPUT
 # x = output of 'lm'

 res <- x$residuals

 MAE <- sum(abs(res))/length(res)

 return(MAE)
}


# -------------------------------------------------------------------- #

predinterval.lm  <- function(x,level=0.95) {

 # INPUT
 # x = output of 'lm'

 pred <- predict(x, interval="prediction", level=level, newdata=x$model)

 return(pred)
}


# -------------------------------------------------------------------- #

jackknife1.lm <- function(x) {

 # INPUT
 # x = output of 'lm'

 f <- x$term
 #f <- x$call[[2]]
 m <- x$model

 n <- dim(m)[1]

 pred <- rep(NA,n)
 names(pred) <- row.names(m)
 for (i in 1:n) pred[i] <- predict(lm(f,m[-i,]), newdata=m[i,-1, drop=FALSE])

 return(pred)
}


# ------------------------------------------------------------------ #

RMSEjk.lm <- function(x) {

 # INPUT
 # x = output of 'lm'

 orig <- x$model[,1]
 pred <- jackknife1.lm(x)
 res <- pred - orig

 RMSE <- sqrt(sum((res)^2)/length(res))

 return(RMSE)
}


# ------------------------------------------------------------------ #

MAEjk.lm <- function(x) {

 # INPUT
 # x = output of 'lm'

 orig <- x$model[,1]
 pred <- jackknife1.lm(x)
 res <- pred - orig

 MAE <- sum(abs(res))/length(res)

 return(MAE)
}


# ------------------------------------------------------------------ #

# mantel.lm <- function (x, Nperm=1000) {
#  
#  f <- x$term
#  m <- x$model
#  betas <- x$coefficients[-1]
#  k <- length(betas)
#  
#  betas0 <- matrix(NA,Nperm,k)
#  for (i in 1:Nperm) {
#   m0 <- data.frame(cbind(sample(m[,1]),m[,-1]))
#   names(m0) <- names(m)
#   x0 <- lm(f, m0)
#   betas0[i,] <- x0$coefficients[-1]
#  }
# 
#  Ps <- rep(NA,k)
#  for (i in 1:k) {
#   Ps[i] <- ecdf(betas0[,i])(betas[i])
#  }
#  Ps <- 1 - abs((sign(betas)>0) - Ps)
#  Ps <- round(Ps,4)
# 
#  output <- c(Ps)
#  names(output) <- paste("P.",names(m)[-1],sep="")
#  return(output)
# }
 
mantel.lm <- function (x, Nperm = 1000) {
    f <- x$term
    m <- x$model
    betas <- x$coefficients[-1]
    mY <- m[, 1]
    mX <- m[, -1, drop=FALSE]
    if (class(mY) != "dist") stop("mantel.lm: distance matrices of class \"dist\" must be tested.")
    M <- as.matrix(mY)
    nsiti <- dim(M)[1]
    k <- length(betas)
    betas0 <- matrix(NA, Nperm, k)
    for (i in 1:Nperm) {
        permutazione <- sample(1:nsiti)
        mYperm <- as.dist(M[permutazione,permutazione])
        m0 <- data.frame(cbind(as.numeric(mYperm), mX))
        names(m0) <- names(m)
        x0 <- lm(f, m0)
        betas0[i, ] <- x0$coefficients[-1]
    }
    Ps <- rep(NA, k)
    for (i in 1:k) {
        Ps[i] <- ecdf(betas0[, i])(betas[i])
    }
    Ps <- 1 - abs((sign(betas) > 0) - Ps)
    Ps <- round(Ps, 4)
    output <- c(Ps)
    names(output) <- paste("P.", names(m)[-1], sep = "")
    return(output)
}
 
