ict.test <- function(y, treat, J = NA, alpha = 0.05, n.draws = 250000, gms = TRUE, pi.table = TRUE){

  if(class(y)=="matrix") design <- "modified"
  else design = "standard"
  
  if (design == "modified") {
    
    ## for the modified design, aggregate to a vector of y_i
    ## and create appropriate treatment indicator
    J <- ncol(y) - 1
    treat <- ifelse(!is.na(y[,J+1]), 1, 0)
    for (j in 1:(J+1))
      for (i in 1:nrow(y))
        y[i,j] <- ifelse(is.na(y[i,j]), 0, y[i,j])
    y <- apply(y, 1, sum)
    
  } else {
    if (is.na(J))
      stop("You must fill in the option J, the number of control items.")
  }
  
  condition.values <- sort(unique(treat))
  treatment.values <- condition.values[condition.values!=0]
  
  if(length(treatment.values) > 1) multi <- TRUE
  else multi <- FALSE
  
  y.all <- y
  treat.all <- treat
  
  bonferroni <- rep(NA, length(treatment.values))
  
  for (curr.treat in treatment.values) {
    
    y <- y.all[treat.all == 0 | treat.all == curr.treat]
    treat <- treat.all[treat.all == 0 | treat.all == curr.treat]
    treat <- treat > 0
    
    t.y1 <- pi.y1 <- rep(NA, J)
    for(j in 0:(J-1)) {
      pi.y1[j+1] <- mean(y[treat==0] <= j) - mean(y[treat==1] <= j)
      try(t.y1[j+1] <- t.test(y[treat==0] <= j, y[treat==1] <= j, alternative = "less")$p.value)
    }
    
    t.y0 <- pi.y0 <- rep(NA, J)
    for(j in 1:J) {
      pi.y0[j] <- mean(y[treat==1] <= j) - mean(y[treat==0] <= (j - 1))
      try(t.y0[j] <- t.test(y[treat==1] <= j, y[treat==0] <= (j - 1), alternative = "less")$p.value)
    }

    n <- length(y)
    
    y.comb <- c(0:(J-1), 1:J)
    t.comb <- c(rep(1, J), rep(0, J))
    
    cond.y1 <- pi.y1 == 0
    cond.y0 <- pi.y0 == 0
    
    if (gms == TRUE) {
      cond.y1 <- (pi.y1 == 0) | (sqrt(n) * pi.y1 / sqrt(var(pi.y1)) > sqrt(log(n)))
      cond.y0 <- (pi.y0 == 0) | (sqrt(n) * pi.y0 / sqrt(var(pi.y0)) > sqrt(log(n)))
    }
    
    rho.pi <- cov.pi <- sd.pi <- matrix(NA, ncol = length(y.comb), nrow = length(y.comb))
    
    sd <- rep(NA, length(y.comb))
    for(j in 1:length(y.comb)) {
      
      if(t.comb[j]==1) sd[j] <- sqrt(((mean(y[treat==1] <= y.comb[j])*(1-mean(y[treat==1] <= y.comb[j])))/sum(treat==1) + (mean(y[treat==0] <= y.comb[j])*(1-mean(y[treat==0] <= y.comb[j])))/sum(treat==0))) 
      
      if(t.comb[j]==0) sd[j] <- sqrt(((mean(y[treat==1] <= (y.comb[j]-1+1))*(1-mean(y[treat==1] <= (y.comb[j]-1+1))))/sum(treat==1) + (mean(y[treat==0] <= (y.comb[j]-1+0))*(1-mean(y[treat==0] <= (y.comb[j]-1+0))))/sum(treat==0)))
      
    }
    
    for(j in 1:length(y.comb)) {
      
      for(k in 1:length(y.comb)) {
        
        if(t.comb[j]==1 & t.comb[k]==1) {
          if(y.comb[j]==y.comb[k]) cov.pi[j,k] <- sd[j]^2
          else if(y.comb[j] < y.comb[k]) cov.pi[j,k] <- (mean(y[treat==1] <= y.comb[j])*(1-mean(y[treat==1] <= y.comb[k]))/sum(treat==1) + mean(y[treat==0] <= y.comb[j])*(1-mean(y[treat==0] <= y.comb[k]))/sum(treat==0)) 
          else cov.pi[j,k] <- 0
          
          if(y.comb[j] <= y.comb[k]) rho.pi[j,k] <- cov.pi[j,k] / (sd[j]*sd[k])
          else rho.pi[j,k] <- 0
        }
        
        if(t.comb[j]==0 & t.comb[k]==0) {
          
          if(y.comb[j]==y.comb[k]) cov.pi[j,k] <- sd[j]^2
          else if(y.comb[j] <= y.comb[k]) cov.pi[j,k] <- (mean(y[treat==1] <= (y.comb[j] -1 +1))*(1-mean(y[treat==1] <= (y.comb[k]-1+1)))/sum(treat==1) + mean(y[treat==0] <= (y.comb[j] - 1 + 0))*(1-mean(y[treat==0] <= (y.comb[k] - 1 + 0)))/sum(treat==0)) 				
          else cov.pi[j,k] <- 0
          
          if(y.comb[j] <= y.comb[k]) rho.pi[j,k] <- cov.pi[j,k] / (sd[j]*sd[k])
          else rho.pi[j,k] <- 0
        }
        
        if(t.comb[j]==0 & t.comb[k]==1) {
          if(y.comb[j] <= y.comb[k]) cov.pi[j,k] <- (-1)*(mean(y[treat==1] <= (y.comb[j]-1+1))*(1-mean(y[treat==1] <= y.comb[k]))/sum(treat==1) + mean(y[treat==0] <= (y.comb[j]-1+0))*(1-mean(y[treat==0] <= y.comb[k]))/sum(treat==0))
          else cov.pi[j,k] <- 0
          
          if(y.comb[j] <= y.comb[k]) rho.pi[j,k] <- cov.pi[j,k] / (sd[j]*sd[k])
          else rho.pi[j,k] <- 0
        }
        
        if(t.comb[j]==1 & t.comb[k]==0) {
          if(y.comb[j] < y.comb[k]) cov.pi[j,k] <- (-1) * ( mean(y[treat==1] <= y.comb[j])*(1-mean(y[treat==1] <= (y.comb[k]-1+1)))/sum(treat==1) + mean(y[treat==0] <= y.comb[j])*(1-mean(y[treat==0] <= (y.comb[k] - 1+ 0)))/sum(treat==0))
          else if(y.comb[j]==y.comb[k]) cov.pi[j,k] <- (-1) * (mean(y[treat==1] <= (y.comb[j]-1+1))*(1-mean(y[treat==1] <= y.comb[k]))/sum(treat==1) + mean(y[treat==0] <= (y.comb[j]-1+0))*(1-mean(y[treat==0] <= y.comb[j]))/sum(treat==0))
          else cov.pi[j,k] <- 0				
          
          if(y.comb[j] <= y.comb[k]) rho.pi[j,k] <- cov.pi[j,k] / (sd[j]*sd[k])
          else rho.pi[j,k] <- 0
        }
      }
    }
    
    for(i in 1:nrow(rho.pi)){
      for(j in 1:ncol(rho.pi)){
        if(y.comb[i]>y.comb[j]) rho.pi[i,j] <- rho.pi[j,i]
        if(y.comb[i]>y.comb[j]) cov.pi[i,j] <- cov.pi[j,i]
      }
    }
    
    if (length(pi.y1) > 0) {
      rho.pi.y1 <- rho.pi[1:length(pi.y1), 1:length(pi.y1)]
      cov.pi.y1 <- cov.pi[1:length(pi.y1), 1:length(pi.y1)]
    } 
    
    if (length(pi.y0) > 0) {
      rho.pi.y0 <- rho.pi[(length(pi.y1)+1):(length(pi.y1) + length(pi.y0)),
                          (length(pi.y1)+1):(length(pi.y1) + length(pi.y0))]
      cov.pi.y0 <- cov.pi[(length(pi.y1)+1):(length(pi.y1) + length(pi.y0)),
                          (length(pi.y1)+1):(length(pi.y1) + length(pi.y0))]
    }

    ## create pi and s.e. table for printing

    y.comb.tb <- c(0:J, 0:J)
    t.comb.tb <- c(rep(1, J+1), rep(0, J+1))
    
    sd.tb <- rep(NA, length(y.comb.tb))
    for(j in 1:length(y.comb.tb)) {
      
      if(t.comb.tb[j]==1) sd.tb[j] <- sqrt(((mean(y[treat==1] <= y.comb.tb[j])*(1-mean(y[treat==1] <= y.comb.tb[j])))/sum(treat==1) + (mean(y[treat==0] <= y.comb.tb[j])*(1-mean(y[treat==0] <= y.comb.tb[j])))/sum(treat==0))) 
      
      if(t.comb.tb[j]==0) sd.tb[j] <- sqrt(((mean(y[treat==1] <= (y.comb.tb[j]-1+1))*(1-mean(y[treat==1] <= (y.comb.tb[j]-1+1))))/sum(treat==1) + (mean(y[treat==0] <= (y.comb.tb[j]-1+0))*(1-mean(y[treat==0] <= (y.comb.tb[j]-1+0))))/sum(treat==0)))
      
    }
    
    pi.y1.tb <- rep(NA, J+1)
    for(j in 0:J) {
      pi.y1.tb[j+1] <- mean(y[treat==0] <= j) - mean(y[treat==1] <= j)
    }
    
    pi.y0.tb <- rep(NA, J+1)
    for(j in 0:J) {
      pi.y0.tb[j+1] <- mean(y[treat==1] <= j) - mean(y[treat==0] <= (j - 1))
    }
    
    tb <- round(rbind(cbind(pi.y1.tb, sd.tb[1:(J+1)]), cbind(pi.y0.tb, sd.tb[(J+2):((J+1)*2)])), 4)
    rownames(tb) <- c(paste("pi(y = ", 0:J, ", t = 1)", sep = ""),
                      paste("pi(y = ", 0:J, ", t = 0)", sep = ""))
    colnames(tb) <- c("est.", "s.e.")
    
    
    ## now reduce the number of tests based on pi = zero and GMS conditions

    pi.y1 <- pi.y1[cond.y1 == FALSE]
    pi.y0 <- pi.y0[cond.y0 == FALSE]
    
    t.y1 <- t.y1[cond.y1 == FALSE]
    t.y0 <- t.y0[cond.y0 == FALSE]
    
    rho.pi.y1 <- rho.pi.y1[cond.y1 == FALSE, cond.y1 == FALSE]
    cov.pi.y1 <- cov.pi.y1[cond.y1 == FALSE, cond.y1 == FALSE]
    rho.pi.y0 <- rho.pi.y0[cond.y0 == FALSE, cond.y0 == FALSE]
    cov.pi.y0 <- cov.pi.y0[cond.y0 == FALSE, cond.y0 == FALSE]
    
    ## begin test calculation
    
    ## calculate p value for sensitive item = 1
    
    if (length(pi.y1) > 1) {
      
      par.y1 <- rep(0, length(pi.y1))
      
      Dmat <- 2*ginv(cov.pi.y1)
      Amat <- diag(length(pi.y1))
      
      if (sum(pi.y1 < 0) > 0) {
        lambda <- solve.QP(Dmat, par.y1, Amat, bvec = -pi.y1)$value
      } else {
        lambda <- 0
      }
      
      w <- rep(NA, length(pi.y1)+1)
      
      rho.pi.y1.partial <- cor2pcor(rho.pi.y1)
      
      if (length(pi.y1)==2) {
        w[3] <- .5 * pi^(-1) * acos(rho.pi.y1[1,2])
        w[2] <- .5
        w[1] <- .5 - .5 * pi^(-1)  * acos(rho.pi.y1[1,2])
      } else if (length(pi.y1)==3) {      
        rho.pi.y1.partial.12.3 <- (rho.pi.y1[1,2] - rho.pi.y1[1,3] *
                                   rho.pi.y1[2,3])/(sqrt(1-rho.pi.y1[1,3]^2) * sqrt(1-rho.pi.y1[2,3]^2))
        rho.pi.y1.partial.13.2 <- (rho.pi.y1[1,3] - rho.pi.y1[1,2] *
                                   rho.pi.y1[3,2])/(sqrt(1-rho.pi.y1[1,2]^2) * sqrt(1-rho.pi.y1[3,2]^2))
        rho.pi.y1.partial.23.1 <- (rho.pi.y1[2,3] - rho.pi.y1[2,1] *
                                   rho.pi.y1[3,1])/(sqrt(1-rho.pi.y1[2,1]^2) * sqrt(1-rho.pi.y1[3,1]^2))
        w[1] <- .25 * pi^(-1) * (2 * pi - acos(rho.pi.y1[1,2]) - acos(rho.pi.y1[1,3]) - acos(rho.pi.y1[2,3]))
        w[2] <- .25 * pi^(-1) * (3 * pi - acos(rho.pi.y1.partial.12.3) -
                                 acos(rho.pi.y1.partial.13.2) - acos(rho.pi.y1.partial.23.1))
        w[3] <- .5 - w[1]
        w[4] <- .5 - w[2]
      } else if (length(pi.y1)==4) {
        w[4] <- .125 * pi^(-1) * (-4 * pi + acos(rho.pi.y1[4,3]) + acos(rho.pi.y1[3,2]) + acos(rho.pi.y1[4,2]))
        w[3] <- .25 * pi^(-2) * ( acos(rho.pi.y1[4,3]) * (pi - acos(rho.pi.y1.partial[2,1])))
        w[2] <- .125 * pi^(-1) * ( 8 * pi - acos(rho.pi.y1[4,3]) + acos(rho.pi.y1[3,2]) + acos(rho.pi.y1[4,2]))
        w[1] <- pmvnorm(mean = pi.y1, sigma = cov.pi.y1, lower = rep(0, length(pi.y1)))
        w[5] <- .5 - w[1] - w[3]
      } else if (length(pi.y1)>4) {
        draws <- mvrnorm(n = n.draws, mu = par.y1, Sigma = cov.pi.y1)
        pi.tilde <- matrix(NA, nrow = n.draws, ncol = length(pi.y1))
        
        for (i in 1:n.draws) {
          if (sum(draws[i,] < 0) > 1) {
            pi.tilde[i,] <- solve.QP(Dmat, par.y1, Amat, bvec = -draws[i,])$solution + draws[i,]
          } else {
            pi.tilde[i,] <- draws[i, ]
          }
        }
        
        pi.tilde.pos.count <-  apply(pi.tilde, 1, function(x) { sum(x > 0) })
        for(k in 0:J)
          w[k+1] <- mean(pi.tilde.pos.count==(J-k))  
      }
      
      p.y1 <- 0
      for(k in 0:length(pi.y1))
        p.y1 <- p.y1 + w[k+1] * pchisq(lambda, df = k, lower.tail = FALSE)
      
    } else if (length(pi.y1) == 1) {
      p.y1 <- t.y1
    }
    
    ## repeat for sensitive item = 0
    
    if (length(pi.y0) > 1) {
      
      par.y0 <- rep(0, length(pi.y0))
      
      Dmat <- 2*ginv(cov.pi.y0)
      Amat <- diag(length(pi.y0))
      
      if (sum(pi.y0 < 0) > 0) {
        lambda <- solve.QP(Dmat, par.y0, Amat, bvec = -pi.y0)$value
      } else {
        lambda <- 0
      }
      
      w <- rep(NA, length(pi.y0)+1)
      
      rho.pi.y0.partial <- cor2pcor(rho.pi.y0)
      
      if (length(pi.y0)==2) {
        w[3] <- .5 * pi^(-1) * acos(rho.pi.y0[1,2])
        w[2] <- .5
        w[1] <- .5 - .5 * pi^(-1)  * acos(rho.pi.y0[1,2])
      } else if (length(pi.y0)==3) {
        rho.pi.y0.partial.12.3 <- (rho.pi.y0[1,2] - rho.pi.y0[1,3] *
                                   rho.pi.y0[2,3])/(sqrt(1-rho.pi.y0[1,3]^2) * sqrt(1-rho.pi.y0[2,3]^2))
        rho.pi.y0.partial.13.2 <- (rho.pi.y0[1,3] - rho.pi.y0[1,2] *
                                   rho.pi.y0[3,2])/(sqrt(1-rho.pi.y0[1,2]^2) * sqrt(1-rho.pi.y0[3,2]^2))
        rho.pi.y0.partial.23.1 <- (rho.pi.y0[2,3] - rho.pi.y0[2,1] *
                                   rho.pi.y0[3,1])/(sqrt(1-rho.pi.y0[2,1]^2) * sqrt(1-rho.pi.y0[3,1]^2))
        w[1] <- .25 * pi^(-1) * (2 * pi - acos(rho.pi.y0[1,2]) - acos(rho.pi.y0[1,3]) - acos(rho.pi.y0[2,3]))
        w[2] <- .25 * pi^(-1) * (3 * pi - acos(rho.pi.y0.partial.12.3) -
                                 acos(rho.pi.y0.partial.13.2) - acos(rho.pi.y0.partial.23.1))
        w[3] <- .5 - w[1]
        w[4] <- .5 - w[2]
      } else if (length(pi.y0)==4) {
        w[4] <- .125 * pi^(-1) * (-4 * pi + acos(rho.pi.y0[4,3]) + acos(rho.pi.y0[3,2]) + acos(rho.pi.y0[4,2]))
        w[3] <- .25 * pi^(-2) * ( acos(rho.pi.y0[4,3]) * (pi - acos(rho.pi.y0.partial[2,1])))
        w[2] <- .125 * pi^(-1) * ( 8 * pi - acos(rho.pi.y0[4,3]) + acos(rho.pi.y0[3,2]) + acos(rho.pi.y0[4,2]))
        w[1] <- pmvnorm(mean = pi.y0, sigma = cov.pi.y0, lower = rep(0, length(pi.y0)))
        w[5] <- .5 - w[1] - w[3]
        
      } else if (length(pi.y0)>4) {
        
        draws <- mvrnorm(n = n.draws, mu = par.y0, Sigma = cov.pi.y0)
        pi.tilde <- matrix(NA, nrow = n.draws, ncol = length(pi.y0))
        
        for(i in 1:n.draws) {
          if (sum(draws[i,] < 0) > 1) {
            pi.tilde[i,] <- solve.QP(Dmat, par.y0, Amat, bvec = -draws[i,])$solution + draws[i,]
          } else {
            pi.tilde[i,] <- draws[i, ]
          }
        }
        
        pi.tilde.pos.count <-  apply(pi.tilde, 1, function(x) { sum(x > 0) })
        for(k in 0:length(pi.y0))
          w[k+1] <- mean(pi.tilde.pos.count==(J-k))
        
      }
      
      p.y0 <- 0
      for(k in 0:length(pi.y0))
        p.y0 <- p.y0 + w[k+1] * pchisq(lambda, df = k, lower.tail = FALSE)
      
    } else if (length(pi.y0) == 1) {
      p.y0 <- t.y0
    }
    
    if ((length(pi.y1) > 0) & (length(pi.y0) > 0))
      bonferroni[curr.treat] <- 2 * min(p.y1, p.y0)
    else if (length(pi.y1) > 0) 
      bonferroni[curr.treat] <- 2 * p.y1 
    else if (length(pi.y0) > 0)
      bonferroni[curr.treat] <- 2 * p.y0
    else
      bonferroni[curr.treat] <-1
    
  } ## end treatment value loop

  names(bonferroni) <- paste("Sensitive Item", treatment.values)

  if (pi.table == FALSE)
    return.object <- list(p = bonferroni)
  else
    return.object <- list(p = bonferroni, pi.table = tb)

  class(return.object) <- "ict.test"
  
  return.object
  
}

print.ict.test <- function(x, ...){
  
  cat("\nTest for List Experiment Design Effects\n\n")

  if (!is.null(x$pi.table)) {
    
    cat("Estimated population proportions \n")
    
    print(x$pi.table)
    
    cat("\n")
      
  }
     
  cat("Bonferroni-corrected p-value\nIf this value is below alpha, you reject the null hypothesis of no design effect. If it is above alpha, you fail to reject the null.\n\n")

  print(x$p)
  
  cat("\n")

  invisible(x)
  
}
