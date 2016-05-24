##
##  PURPOSE:   Pseudo goodness-of-fit test for a normal mixture
##             * default method
##
##  AUTHOR:    Arnost Komarek
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   20/08/2009
##
##  FUNCTIONS: NMixPseudoGOF.default.R
##
## ==================================================================

## *************************************************************
## NMixPseudoGOF.default
## *************************************************************
NMixPseudoGOF.default <- function(x, scale, w, mu, Sigma, breaks, nbreaks=10, digits=3, ...)
{  
  ##### Dimension of the normal mixture
  ##### ===============================================
  if (is.numeric(x)) x <- data.frame(x=x)
  if (is.matrix(x)) x <- as.data.frame(x)
  if (!is.data.frame(x)) stop("x must be a data.frame")
  p <- ncol(x)
  if (p < 1) stop("length of breaks must be 1 or more")

  ##### Remove rows with NA's from data
  ##### ===============================================  
  row.NA <- apply(is.na(x), 1, sum) > 0
  x <- x[!row.NA,]
  
  ##### scale
  ##### ===============================================  
  if (missing(scale)) scale <- list(shift=rep(0, p), scale=rep(1, p))
  if (!is.list(scale)) stop("scale must be a list")
  if (length(scale) != 2) stop("scale must have 2 components")
  inscale <- names(scale)  
  iscale.shift <- match("shift", inscale, nomatch=NA)
  iscale.scale <- match("scale", inscale, nomatch=NA)
  if (is.na(iscale.shift)) stop("scale$shift is missing")
  if (length(scale$shift) == 1) scale$shift <- rep(scale$shift, p)
  if (length(scale$shift) != p) stop(paste("scale$shift must be a vector of length ", p, sep=""))    
  if (is.na(iscale.scale)) stop("scale$scale is missing")
  if (length(scale$scale) == 1) scale$scale <- rep(scale$scale, p)
  if (length(scale$scale) != p) stop(paste("scale$scale must be a vector of length ", p, sep=""))
  if (any(scale$scale <= 0)) stop("all elements of scale$scale must be positive")

  ##### Number of mixture components
  ##### ===============================================  
  K <- length(w)

  ##### Check mixture weights
  ##### ===============================================  
  if (any(w < 0) | any(w > 1)) stop("weights must lie between zero and 1")
  if (abs(sum(w) - 1) > 1e-5) warning("sum of weights differs from 1")

  ##### Check mixture means and variances
  ##### (make them numeric if p=1)
  ##### Adjust means and variances for scaling
  ##### Compute overall mean and variance of the mixture (to determine sensible break values if not given)
  ##### ====================================================================================================  
  if (p == 1){

    ### Check
    if (length(mu) != K) stop("incorrect mu")
    mu <- as.numeric(mu)
    if (is.list(Sigma)){
      if (any(sapply(Sigma, length) != 1)) stop("incorrect Sigma")
      Sigma <- unlist(Sigma)
    }
    if (length(Sigma) != K) stop("incorrect Sigma")

    ### Scale adjustment
    mu <- mu * scale$scale + scale$shift
    Sigma <- Sigma * scale$scale^2

    ### Overall mean and variance
    Emix <- sum(w * mu)
    Vmix <- sum(w * (Sigma + (mu - Emix)^2))    
  }else{

    ### Check
    if (K == 1){
      if (length(mu) != p) stop("incorrect mu")
      mu <- matrix(mu, nrow=K, ncol=p)
      if (!is.list(Sigma)) Sigma <- list(Sigma)    
    }
    if (nrow(mu) != K) stop(paste("mu must have ", K, " rows", sep=""))
    if (ncol(mu) != p) stop(paste("mu must have ", p, " columns", sep=""))
    if (any(!sapply(Sigma, is.matrix))) stop("all Sigma's must be matrices")
    if (any(sapply(Sigma, nrow) != p)) stop(paste("all Sigma's must have ", p, " rows", sep=""))
    if (any(sapply(Sigma, ncol) != p)) stop(paste("all Sigma's must have ", p, " columns", sep=""))

    ### Scale adjustment, overall mean and variance
    mu <- mu * matrix(rep(scale$scale, K), nrow=K, ncol=p, byrow=TRUE) + matrix(rep(scale$shift, K), nrow=K, ncol=p, byrow=TRUE)
    Emix <- apply(matrix(rep(w, p), nrow=K, ncol=p) * mu, 2, sum)
    Vmix <- matrix(0, nrow=p, ncol=p)
    for (k in 1:K){
      Sigma[[k]] <- diag(scale$scale) %*% Sigma[[k]] %*% diag(scale$scale)
      Vmix <- Vmix + w[k] * (Sigma[[k]] + matrix(mu[k,] - Emix, ncol=1) %*% matrix(mu[k,] - Emix, nrow=1)) 
    }
  }
  
  ##### Breaks
  ##### ===============================================  
  if (length(nbreaks) == 1) nbreaks <- rep(nbreaks, p)
  if (length(nbreaks) != p) stop(paste("nbreaks must be a vector of length ", p, sep=""))

  if (length(digits) == 1) digits <- rep(digits, p)
  if (length(digits) != p) stop(paste("digits must be a vector of length ", p, sep=""))

  csigma <- 2
  if (missing(breaks)){    
    breaks <- list()
    if (p == 1){
      rangeGrid <- Emix + c(-csigma, csigma)*sqrt(Vmix)
      if (nbreaks[1] == 1) breaks[[1]] <- Emix
      else                 breaks[[1]] <- round(seq(rangeGrid[1], rangeGrid[2], length=nbreaks), digits=digits)
    }else{     
      for (i in 1:p){
        rangeGrid <- Emix[i] + c(-csigma, csigma)*sqrt(Vmix[i, i])
        if (nbreaks[i] == 1) breaks[[i]] <- Emix[i]
        else                 breaks[[i]] <- round(seq(rangeGrid[1], rangeGrid[2], length=nbreaks[i]), digits=digits[i])
      }
    }
  }
  
  if (is.numeric(breaks)) breaks <- list(x1=breaks)  
  if (!is.list(breaks)) stop("breaks must be a list")
  if (length(breaks) != p) stop(paste("breaks must be a list of length ", p, sep=""))
  if (is.null(names(breaks))) names(breaks) <- paste("x", (1:p), sep="")  

  ##### Order breaks, add Inf, create labels of intervals, create mid-points
  ##### ===========================================================================  
  labels <- midpoints <- breaksInf <- list()
  for (i in 1:p){
    breaks[[i]] <- breaks[[i]][order(breaks[[i]])]
    if (length(breaks[[i]]) == 1){
      if (p == 1) midpoints[[i]] <- breaks[[i]] + c(-1, 1) * sqrt(Vmix)
      else        midpoints[[i]] <- breaks[[i]] + c(-1, 1) * sqrt(Vmix[i, i])
    }else{
      aver.dist <- mean(breaks[[i]][2:length(breaks[[i]])] - breaks[[i]][1:(length(breaks[[i]]) - 1)])
      mids <- 0.5 * (breaks[[i]][2:length(breaks[[i]])] + breaks[[i]][1:(length(breaks[[i]]) - 1)])
      midpoints[[i]] <- c(breaks[[i]][1] - aver.dist, mids, breaks[[i]][length(breaks[[i]])] + aver.dist)
    }  
    
    labLeft <- paste("(", c(-Inf, breaks[[i]]), sep="")
    labRight <- c(paste(breaks[[i]], "]", sep=""), "Inf)")
    labels[[i]] <- paste(labLeft, ", ", labRight, sep="")    
    breaksInf[[i]] <- c(breaks[[i]], Inf)      
  }  

  nBreak <- sapply(breaks, length)
  #if (any(nBreak <= 1)) stop("there must be at least two breaks in each margin")

  
  ##### For each observation, determine interval in which it lies
  ##### - store it in index matrix
  ##### ============================================================================  
  index <- matrix(0, nrow=nrow(x), ncol=p)
  for (i in 1:p){
    tmpx <- matrix(rep(x[,i], length(breaks[[i]])), nrow=length(x[,i]))
    tmpb <- matrix(rep(breaks[[i]], each=length(x[,i])), nrow=length(x[,i]))
    index[,i] <- apply(tmpx > tmpb, 1, sum) + 1                                   ## takes values 1, ..., nBreak[i] + 1
  }    

  ##### Compute pseudo chi-squared statistics for each univariate margin
  ##### ============================================================================  
  RESULT <- CHISQ <- list()

  F1 <- list()
  RESULT[[1]] <- list()
  NPARAM <- (K - 1) + K + K
  for (i in 1:p){
    Mean <- mu[,i]
    StdDev <- sqrt(sapply(Sigma, "[", i, i))
    F <- matrix(0, nrow=K, ncol=nBreak[i])
    for (k in 1:K){
      F[k,] <- w[k]*pnorm(breaks[[i]], mean=Mean[k], sd=StdDev[k])
    }
    F1[[i]] <- apply(F, 2, sum)
    PP <- c(F1[[i]], 1) - c(0, F1[[i]])                                      ## probabilities of bins
    
    EXPECT <- nrow(x) * PP                                                   ## expected counts
    
    TAB <- table(index[,i])
    OBSERVED <- rep(0, length(labels[[i]]))
    names(OBSERVED) <- paste(1:length(labels[[i]]))
    OBSERVED[names(TAB)] <- TAB                                              ## observed counts
    
    RESID <- (OBSERVED - EXPECT) / sqrt(EXPECT)
    if (i == 1) CHISQ[[1]] <- data.frame(ChiSq = sum(RESID^2), nBin = length(labels[[i]]), nParam = NPARAM)
    else        CHISQ[[1]] <- rbind(CHISQ[[1]], data.frame(ChiSq = sum(RESID^2), nBin = length(labels[[i]]), nParam = NPARAM))    
    RESULT[[1]][[i]] <- data.frame(Ind1 = 1:(nBreak[i] + 1), Mid1 = midpoints[[i]], Bin1 = labels[[i]], Observed = OBSERVED, Expected = EXPECT, Resid = RESID)
  }
  rownames(CHISQ[[1]]) <- names(RESULT[[1]]) <- paste(1:p)
    
  #layout(autolayout(p))
  #for (i in 1:p) plot(RESULT[[1]][[i]][,"Mid1"], abs(RESULT[[1]][[i]][,"Resid"]), type="h", col="red", lwd=2)

  ##### Compute pseudo chi-squared statistics for each bivariate margin
  ##### ====================================================================================  
  if (p > 1){
    RESULT[[2]] <- list()
    NPARAM <- (K - 1) + K * 2 + K * 3
    NAAM <- character()
    for (i in 1:(p - 1)){
      for (j in (i+1):p){
        NAAM <- c(NAAM, paste(i, "-", j, sep=""))
        nBin <- length(labels[[i]]) * length(labels[[j]])
        Bin1 <- rep(labels[[i]], length(labels[[j]]))
        Bin2 <- rep(labels[[j]], each=length(labels[[i]]))
        
        Break1 <- rep(breaks[[i]], nBreak[j])
        Break2 <- rep(breaks[[j]], each=nBreak[i])
        Grid <- cbind(Break1, Break2)
        
        ### Distribution function in verteces of rectangles
        ### ----------------------------------------------------
        FF <- rep(0, nrow(Grid))
        for (g in 1:nrow(Grid)) for (k in 1:K) FF[g] <- FF[g] + w[k] * mnormt::pmnorm(Grid[g,], mean=mu[k, c(i, j)], varcov=Sigma[[k]][c(i, j), c(i, j)])
        FF <- matrix(FF, nrow=nBreak[i], ncol=nBreak[j])

        ### Rectangle probabilities
        ### ----------------------------------------------------
        PP <- numeric(nBin)

            ## P(X1 in A, X2 <= Break2[1])
        PP[1] <- FF[1, 1]                                                                       ## P(X1<=Break1[1], X2<=Break2[1])
        if (nBreak[i] > 1) PP[2:nBreak[i]]   <- FF[2:nBreak[i], 1] - FF[1:(nBreak[i] - 1), 1]   ## First row of PP matrix with X2<=Break2[1] and finite X1
        PP[nBreak[i] + 1] <- F1[[j]][1] - FF[nBreak[i], 1]                                      ## P(X1>Break1[nBreak[i]], X2<=Break2[1])
        Filled <- nBreak[i] + 1

            ## P(X1 in A, X2 in B)
        if (nBreak[j] > 1){
          for (d2 in 2:nBreak[j]){
            PP[Filled + 1] <- FF[1, d2] - FF[1, d2 - 1]
            if (nBreak[i] > 1) PP[Filled + (2:nBreak[i])] <- FF[2:nBreak[i], d2] - FF[2:nBreak[i], d2 - 1] - FF[1:(nBreak[i] - 1), d2] + FF[1:(nBreak[i] - 1), d2 - 1]
            PP[Filled + nBreak[i] + 1] <- F1[[j]][d2] - FF[nBreak[i], d2] - F1[[j]][d2 - 1] + FF[nBreak[i], d2 - 1]
            Filled <- Filled + nBreak[i] + 1            
          }            
        }

            ## P(X1 in A, X2 > Break2[nBreak[j]]
        PP[Filled + 1] <- F1[[i]][1] - FF[1, nBreak[j]]
        if (nBreak[i] > 1) PP[Filled + (2:nBreak[i])] <- F1[[i]][2:nBreak[i]] - FF[2:nBreak[i], nBreak[j]] - F1[[i]][1:(nBreak[i] - 1)] + FF[1:(nBreak[i] - 1), nBreak[j]]
        PP[Filled + nBreak[i] + 1] <- 1 - F1[[i]][nBreak[i]] - F1[[j]][nBreak[j]] + FF[nBreak[i], nBreak[j]]

        ### Expected counts
        ### -----------------------------------------------------
        EXPECT <- nrow(x) * PP
        
        ### Observed counts
        ### -----------------------------------------------------
        TAB <- table(index[, i], index[, j])
        OBSERVED <- matrix(0, nrow=length(labels[[i]]), ncol=length(labels[[j]]))        
        rownames(OBSERVED) <- paste(1:length(labels[[i]]))
        colnames(OBSERVED) <- paste(1:length(labels[[j]]))        
        OBSERVED[rownames(TAB), colnames(TAB)] <- TAB
        OBSERVED <- as.numeric(OBSERVED)

        ### Pearson's residuals, put everything together
        ### -----------------------------------------------------
        RESID <- (OBSERVED - EXPECT) / sqrt(EXPECT)
        if (i == 1 & j == 2){
          CHISQ[[2]] <- data.frame(ChiSq = sum(RESID^2), nBin = nBin, nParam = NPARAM)
        }else{  
          CHISQ[[2]] <- rbind(CHISQ[[2]], data.frame(ChiSq = sum(RESID^2), nBin = nBin, nParam = NPARAM))
        }  
        RESULT[[2]][[paste(i, "-", j, sep="")]] <- data.frame(Ind1 = rep(1:(nBreak[i] + 1), nBreak[j] + 1),       Ind2 = rep(1:(nBreak[j] + 1), each=nBreak[i] + 1),
                                                              Mid1 = rep(midpoints[[i]], length(midpoints[[j]])), Mid2 = rep(midpoints[[j]], each=length(midpoints[[i]])),
                                                              Bin1 = rep(labels[[i]], length(labels[[j]])),       Bin2 = rep(labels[[j]], each=length(labels[[i]])),
                                                              Observed = OBSERVED, Expected = EXPECT, Resid = RESID)
        
        if (i == 1 & j == 3){
          COL <- rev(heat_hcl(33, c.=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))
          PP2 <- matrix(PP, nrow=length(labels[[i]]), ncol=length(labels[[j]]))
          RESID2 <- matrix(abs(RESID), nrow=length(labels[[i]]), ncol=length(labels[[j]]))

          layout(autolayout(3))
          par(bty="n")
          image(midpoints[[i]], midpoints[[j]], PP2, col=COL, xlab=paste("x", i, sep=""), ylab=paste("x", j, sep=""))
          image(breaks[[i]], breaks[[j]], FF, col=COL, xlab=paste("x", i, sep=""), ylab=paste("x", j, sep=""))

          IMBREAK <- seq(0.5, 5, by=0.5)
          IMCOL <- rev(heat_hcl(length(IMBREAK) + 1, c.=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))
          image(midpoints[[i]], midpoints[[j]], RESID2, col=IMCOL, xlab=paste("x", i, sep=""), ylab=paste("x", j, sep=""))          
        }
      }
    }

    rownames(CHISQ[[2]]) <- NAAM
  }  
  
  LTp <- p * (p + 1)/2

  RET <- list(result=RESULT, chisq=CHISQ)  
  return(RET)  
}
