genLogistik<-function (data, member, match="score",anchor = 1:ncol(data), type = "both", 
    criterion = "LRT") 
{
    R2 <- function(m, n) 1 - (exp(-m$null.deviance/2)/exp(-m$deviance/2))^(2/n)
    R2max <- function(m, n) 1 - (exp(-m$null.deviance/2))^(2/n)
    R2DIF <- function(m, n) R2(m, n)/R2max(m, n)
    dev <- R2full<-R2simple<- deltaR <- rep(NA, ncol(data))
    nGroup <- length(unique(member)) - 1
    mFull <- mSimple <- matrix(0, ncol(data), 2 + 2 * nGroup)
    if (type == "udif") {
        sigmaMat <- rep(NA, (2 + nGroup) * (2 + nGroup) * ncol(data))
        dim(sigmaMat) <- c(2 + nGroup, 2 + nGroup, ncol(data))
    }
    else {
        sigmaMat <- rep(NA, (2 + 2 * nGroup) * (2 + 2 * nGroup) * 
            ncol(data))
        dim(sigmaMat) <- c(2 + 2 * nGroup, 2 + 2 * nGroup, ncol(data))
    }
    for (item in 1:ncol(data)) {
if (match[1]=="score"){
        data2 <- data[, anchor]
        if (sum(anchor == item) == 0) 
            data2 <- cbind(data2, data[, item])
        SCORES <- rowSums(data2, na.rm = TRUE)
}
else SCORES<-match
        GROUP <- as.factor(member)
      ITEM<-data[,item]
    
        m0 <- switch(type, both = glm(ITEM ~ SCORES * GROUP, family = "binomial"), 
              udif = glm(ITEM ~ SCORES + GROUP, family = "binomial"), 
            nudif = glm(ITEM ~ SCORES * GROUP, family = "binomial"))
        m1 <- switch(type, both = glm(ITEM ~ SCORES, family = "binomial"),
              udif = glm(ITEM ~ SCORES, family = "binomial"),
              nudif = glm(ITEM ~ SCORES + GROUP, family = "binomial"))
        if (criterion == "LRT") {
            dev[item] <- deviance(m1) - deviance(m0)
            covMat <- summary(m0)$cov.scaled
            sigmaMat[, , item] <- covMat
        }
        else {
            if (criterion != "Wald") 
                stop("'criterion' must be either 'LRT' or Wald'", 
                  call. = FALSE)
            else {
                coeff <- as.numeric(coefficients(m0))
                covMat <- summary(m0)$cov.scaled
                sigmaMat[, , item] <- covMat
                if (type == "udif") {
                  C <- matrix(0, nGroup, length(coeff))
                  for (tt in 1:nGroup) C[tt, 2 + tt] <- 1
                }
                else {
                  if (type == "nudif") {
                    C <- matrix(0, nGroup, length(coeff))
                    for (tt in 1:nGroup) C[tt, 2 + nGroup + tt] <- 1
                  }
                  else {
                    C <- matrix(0, nGroup * 2, length(coeff))
                    for (tt in 1:(2 * nGroup)) C[tt, 2 + tt] <- 1
                  }
                }
                dev[item] <- t(C %*% coeff) %*% solve(C %*% covMat %*% 
                  t(C)) %*% C %*% coeff
            }
        }
        R2full[item] <- R2DIF(m0, nrow(data)) 
        R2simple[item] <- R2DIF(m1, nrow(data)) 
        deltaR[item] <- R2DIF(m0, nrow(data)) - R2DIF(m1, nrow(data))
        mFull[item, 1:length(m0$coefficients)] <- m0$coefficients
        mSimple[item, 1:length(m1$coefficients)] <- m1$coefficients
    }
    names <- c("(Intercept)", "SCORE")
GE <- sort(unique(member))
    for (i in 2:length(GE)) names <- c(names, paste("GROUP", 
        GE[i], sep = ""))
    for (i in 2:length(GE)) names <- c(names, paste("SCORE:GROUP", 
        GE[i], sep = ""))
    colnames(mFull) <- colnames(mSimple) <- names
    res <- list(stat = dev, R2M0=R2full,R2M1=R2simple,deltaR2 = deltaR, parM0 = mFull, 
        parM1 = mSimple, covMat = sigmaMat, criterion = criterion,
match=ifelse(match[1]=="score","score","matching variable"))
    return(res)
}
