##
## Zivot-Andrews Test
##
ur.za <- function(y, model=c("intercept", "trend", "both"), lag=NULL){
  y <- na.omit(as.vector(y))
  n <- length(y)
  model <- match.arg(model)
  if(is.null(lag)) lag <- 0
  lag <- as.integer(lag)
  if(length(lag) > 1 || lag < 0){
    warning("\nPlease, specify maximal number of lags for differenced series as positive integer; lag=1 is now used.")
    lag <- 1}
  datmat <- matrix(NA, n, lag + 3)
  if(n < ncol(datmat) + 2){
    stop("\nInsufficient number of obeservations.")}
  idx <- 1:(n-1)
  trend <- seq(1, n)
  datmat[,1] <- y
  datmat[,2] <- c(NA, y)[1:n]
  datmat[,3] <- trend
  datmat <- as.data.frame(datmat)
  colnames(datmat)[1:3] <- c("y", "y.l1", "trend")
  if(lag > 0){
    for(i in 1:lag){
      datmat[ , i + 3] <- c(rep(NA, i + 1), diff(y))[1:n]
    }
  colnames(datmat) <- c("y", "y.l1", "trend", paste("y.dl", 1:lag, sep=""))
  }
  if(model=="intercept"){
    roll <- function(z){
      du <- c(rep(0, z), rep(1, (n-z)))
      rollmat <- cbind(datmat, du)
      roll.reg <- coef(summary(lm(rollmat)))
      (roll.reg[2,1]-1.0)/roll.reg[2,2]
    }
    roll.stat <- sapply(idx, roll)
    cval <- c(-5.34, -4.8, -4.58)
    bpoint <- which.min(roll.stat)
    du <- c(rep(0, bpoint), rep(1, (n-bpoint)))
    testmat <- cbind(datmat, du)
    test.reg <- lm(testmat) 
  }else if(model=="trend"){
    roll <- function(z){
      dt <- c(rep(0, z), 1:(n-z))
      rollmat <- cbind(datmat, dt)
      roll.reg <- coef(summary(lm(rollmat)))
      (roll.reg[2,1]-1.0)/roll.reg[2,2]
    }
    roll.stat <- sapply(idx, roll)
    cval <- c(-4.93, -4.42, -4.11)
    bpoint <- which.min(roll.stat)
    dt <- c(rep(0, bpoint), 1:(n-bpoint))
    testmat <- cbind(datmat, dt)
    test.reg <- lm(testmat) 
  }else if(model=="both"){
    test.reg <- lm(datmat)
    roll <- function(z){
      du <- c(rep(0, z), rep(1, (n-z)))
      dt <- c(rep(0, z), 1:(n-z))
      rollmat <- cbind(datmat, du, dt)
      roll.reg <- coef(summary(lm(rollmat)))
      (roll.reg[2,1]-1.0)/roll.reg[2,2]
    }
    roll.stat <- sapply(idx, roll)
    cval <- c(-5.57, -5.08, -4.82)
    bpoint <- which.min(roll.stat)
    du <- c(rep(0, bpoint), rep(1, (n-bpoint)))
    dt <- c(rep(0, bpoint), 1:(n-bpoint))
    testmat <- cbind(datmat, du, dt)
    test.reg <- lm(testmat) 
  }
  teststat <- roll.stat[bpoint]
  new("ur.za", y=y, model=model, lag=lag, teststat=teststat, cval=cval, bpoint=bpoint, tstats=roll.stat, res=test.reg$residuals, testreg=test.reg, test.name="Zivot-Andrews")
}
