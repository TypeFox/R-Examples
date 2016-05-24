vector.splithalf <-
function(x, set, typ="cor", sims=100, graph=TRUE, CI=.95, minval=-1.0, seed=2) {
  comb.data <- data.frame(cbind(x, set))
  comp.data <- subset(comb.data, complete.cases(comb.data))
  N <- nrow(comp.data)
  half.N <- N/2
  coded <- seq(1:N)
  LL <- (1 - CI) / 2
  UL <- 1 - LL

  if(seed!=F) {set.seed(seed)}

    if(typ=="cor") {
      storer <- rep(0, sims)
        for(i in 1:sims) {
          rand.assign <- sample(coded, N, FALSE)
          rvec.1 <- as.vector(cor(comp.data[rand.assign <= half.N,1], comp.data[rand.assign <= half.N,-1]))
          rvec.2 <- as.vector(cor(comp.data[rand.assign > half.N,1], comp.data[rand.assign > half.N,-1]))
          storer[i] <- (2*cor(rvec.1, rvec.2)) / (cor(rvec.1, rvec.2)+1)
        }
      Up.r <- mean(storer)
      if(graph==T) {
        op <- par(las=1, font.main=1)
        hist(storer, main="Histogram of Split-Half rs for Corrs", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
        abline(v=Up.r, col="red")
        abline(v=quantile(storer, LL), col="orange")
        abline(v=quantile(storer, UL), col="orange")
      }
      Up.r <- ifelse(Up.r >= minval, Up.r, minval)
      out <- cbind(N, Up.r, sd(storer), quantile(storer, LL), quantile(storer, UL))
      colnames(out) <- c("N", "Split-Half r", "SE", "LL", "UL")
      rownames(out) <- c("Corr")
    }

    if(typ=="betas") {
       storeXY <- rep(0, sims)
       storeYX <- rep(0, sims)
       Z <- scale2(comp.data)
        for(i in 1:sims) {
          rand.assign <- sample(coded, N, FALSE)
          bvec.1xy <- as.vector(cor(Z[rand.assign <= half.N,1], Z[rand.assign <= half.N,-1])) * (sd(Z[rand.assign <= half.N,-1]) / sd(Z[rand.assign <= half.N,1]))
          bvec.2xy <- as.vector(cor(Z[rand.assign > half.N,1], Z[rand.assign > half.N,-1])) * (sd(Z[rand.assign > half.N,-1]) / sd(Z[rand.assign > half.N,1]))
          storeXY[i] <- (2*cor(bvec.1xy, bvec.2xy)) / (cor(bvec.1xy, bvec.2xy)+1)
          bvec.1yx <- as.vector(cor(Z[rand.assign <= half.N,1], Z[rand.assign <= half.N,-1])) * (sd(Z[rand.assign <= half.N,1]) / sd(Z[rand.assign <= half.N,-1]))
          bvec.2yx <- as.vector(cor(Z[rand.assign > half.N,1], Z[rand.assign > half.N,-1])) * (sd(Z[rand.assign > half.N,1]) / sd(Z[rand.assign > half.N,-1]))
          storeYX[i] <- (2*cor(bvec.1yx, bvec.2yx)) / (cor(bvec.1yx, bvec.2yx)+1)
        }

      Up.rXY <- mean(storeXY)
      Up.rYX <- mean(storeYX)

      if(graph==T) {
        op <- par(mfrow=c(2,1), font.main=1, las=1)
        hist(storeXY, main="Histogram of Split-Half rs for X-Y Betas", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
        abline(v=Up.rXY, col="red")
        abline(v=quantile(storeXY, LL), col="orange")
        abline(v=quantile(storeXY, UL), col="orange")
        hist(storeYX, main="Histogram of Split-Half rs for Y-X Betas", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
        abline(v=Up.rYX, col="red")
        abline(v=quantile(storeYX, LL), col="orange")
        abline(v=quantile(storeYX, UL), col="orange")
      }

      Up.rXY <- ifelse(Up.rXY >= minval, Up.rXY, minval)
      Up.rYX <- ifelse(Up.rYX >= minval, Up.rYX, minval)
      outXY <- cbind(N, Up.rXY, sd(storeXY), quantile(storeXY, LL), quantile(storeXY, UL))
      outYX <- cbind(N, Up.rYX, sd(storeYX), quantile(storeYX, LL), quantile(storeYX, UL))
      out <- rbind(outXY, outYX)
      colnames(out) <- c("N", "Split-half r", "SE", "LL", "UL")
      rownames(out) <- c("X-Y betas", "Y-X betas")
    }

    if(typ=="XY") {
       Z <- scale2(comp.data)
       storeXY <- rep(0, sims)
        for(i in 1:sims) {
          rand.assign <- sample(coded, N, FALSE)
          bvec.1xy <- as.vector(cor(Z[rand.assign <= half.N,1], Z[rand.assign <= half.N,-1])) * (sd(Z[rand.assign <= half.N,-1]) / sd(Z[rand.assign <= half.N,1]))
          bvec.2xy <- as.vector(cor(Z[rand.assign > half.N,1], Z[rand.assign > half.N,-1])) * (sd(Z[rand.assign > half.N,-1]) / sd(Z[rand.assign > half.N,1]))
          storeXY[i] <- (2*cor(bvec.1xy, bvec.2xy)) / (cor(bvec.1xy, bvec.2xy)+1)
        }

      Up.rXY <- mean(storeXY)
      if(graph==T) {
        op <- par(font.main=1, las=1)
        hist(storeXY, main="Histogram of Split-Half rs for X-Y Betas", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
        abline(v=Up.rXY, col="red")
        abline(v=quantile(storeXY, LL), col="orange")
        abline(v=quantile(storeXY, UL), col="orange")
      }
      
      Up.rXY <- ifelse(Up.rXY >= minval, Up.rXY, minval)
      out <- cbind(N, Up.rXY, sd(storeXY), quantile(storeXY, LL), quantile(storeXY, UL))
      colnames(out) <- c("N", "Split-half r", "SE", "LL", "UL")
      rownames(out) <- c("X-Y betas")
    }

    if(typ=="YX") {
       storeYX <- rep(0, sims)
       Z <- scale2(comp.data)
        for(i in 1:sims) {
          rand.assign <- sample(coded, N, FALSE)
          bvec.1yx <- as.vector(cor(Z[rand.assign <= half.N,1], Z[rand.assign <= half.N,-1])) * (sd(Z[rand.assign <= half.N,1]) / sd(Z[rand.assign <= half.N,-1]))
          bvec.2yx <- as.vector(cor(Z[rand.assign > half.N,1], Z[rand.assign > half.N,-1])) * (sd(Z[rand.assign > half.N,1]) / sd(Z[rand.assign > half.N,-1]))
          storeYX[i] <- (2*cor(bvec.1yx, bvec.2yx)) / (cor(bvec.1yx, bvec.2yx)+1)
        }

      Up.rYX <- mean(storeYX)

      if(graph==T) {
        op <- par(font.main=1, las=1)
        hist(storeYX, main="Histogram of Split-Half rs for Y-X Betas", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
        abline(v=Up.rYX, col="red")
        abline(v=quantile(storeYX, LL), col="orange")
        abline(v=quantile(storeYX, UL), col="orange")
      }

      Up.rYX <- ifelse(Up.rYX >= minval, Up.rYX, minval)
      out <- cbind(N, Up.rYX, sd(storeYX), quantile(storeYX, LL), quantile(storeYX, UL))
      colnames(out) <- c("N", "Split-half r", "SE", "LL", "UL")
      rownames(out) <- c("Y-X betas")
    }

    if(typ=="all") {
      storer <- rep(0, sims)
      storeXY <- rep(0, sims)
      storeYX <- rep(0, sims)
      Z <- scale2(comp.data)
        for(i in 1:sims) {
          rand.assign <- sample(coded, N, FALSE)
          rvec.1 <- as.vector(cor(comp.data[rand.assign <= half.N,1], comp.data[rand.assign <= half.N,-1]))
          rvec.2 <- as.vector(cor(comp.data[rand.assign > half.N,1], comp.data[rand.assign > half.N,-1]))
          storer[i] <- (2*cor(rvec.1, rvec.2)) / (cor(rvec.1, rvec.2)+1)
          bvec.1xy <- as.vector(cor(Z[rand.assign <= half.N,1], Z[rand.assign <= half.N,-1])) * (sd(Z[rand.assign <= half.N,-1]) / sd(Z[rand.assign <= half.N,1]))
          bvec.2xy <- as.vector(cor(Z[rand.assign > half.N,1], Z[rand.assign > half.N,-1])) * (sd(Z[rand.assign > half.N,-1]) / sd(Z[rand.assign > half.N,1]))
          storeXY[i] <- (2*cor(bvec.1xy, bvec.2xy)) / (cor(bvec.1xy, bvec.2xy)+1)
          bvec.1yx <- as.vector(cor(Z[rand.assign <= half.N,1], Z[rand.assign <= half.N,-1])) * (sd(Z[rand.assign <= half.N,1]) / sd(Z[rand.assign <= half.N,-1]))
          bvec.2yx <- as.vector(cor(Z[rand.assign > half.N,1], Z[rand.assign > half.N,-1])) * (sd(Z[rand.assign > half.N,1]) / sd(Z[rand.assign > half.N,-1]))
          storeYX[i] <- (2*cor(bvec.1yx, bvec.2yx)) / (cor(bvec.1yx, bvec.2yx)+1)
        }

      Up.r <- mean(storer)
      Up.rXY <- mean(storeXY)
      Up.rYX <- mean(storeYX)

      if(graph==T) {
        op <- par(mfrow=c(2,2), las=1, font.main=1)
        hist(storer, main="Histogram of Split-Half rs for Corrs", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
        abline(v=Up.r, col="red")
        abline(v=quantile(storer, LL), col="orange")
        abline(v=quantile(storer, UL), col="orange")
        hist(storeXY, main="Histogram of Split-Half rs for X-Y Betas", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
        abline(v=Up.rXY, col="red")
        abline(v=quantile(storeXY, LL), col="orange")
        abline(v=quantile(storeXY, UL), col="orange")
        hist(storeYX, main="Histogram of Split-Half rs for Y-X Betas", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
        abline(v=Up.rYX, col="red")
        abline(v=quantile(storeYX, LL), col="orange")
        abline(v=quantile(storeYX, UL), col="orange")
      }

      Up.r <- ifelse(Up.r >= minval, Up.r, minval)
      Up.rXY <- ifelse(Up.rXY >= minval, Up.rXY, minval)
      Up.rYX <- ifelse(Up.rYX >= minval, Up.rYX, minval)
      outcorr <- cbind(N, Up.r, sd(storer), quantile(storer, LL), quantile(storer, UL))
      outXY <- cbind(N, Up.rXY, sd(storeXY), quantile(storeXY, LL), quantile(storeXY, UL))
      outYX <- cbind(N, Up.rYX, sd(storeYX), quantile(storeYX, LL), quantile(storeYX, UL))
      out <- rbind(outcorr, outXY, outYX)
      colnames(out) <- c("N", "Split-Half r", "SE", "LL", "UL")
      rownames(out) <- c("Corr", "X-Y betas", "Y-X betas")
    }
  return(out)
}
