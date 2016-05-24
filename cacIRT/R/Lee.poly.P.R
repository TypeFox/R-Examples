Lee.poly.P <-function(cutscore, Pij, theta)	 {
  
  ut <- theta
  nn <- length(ut)
  nk <- dim(Pij)[2]	
  ni <- dim(Pij)[3]
  sc <- ni + 1
  nc <- length(cutscore)
  
  exp.TS <- rowSums(sapply(1:ni, function(i) Pij[,,i] %*% 0:(nk-1)))
  
  rec.mat <- gen.rec.raw(Pij, ut)
  
  esacc <- escon <- matrix(NA, nc, nn, dimnames = list(paste("cut at", cutscore), round(ut, 3)))
  
  for (j in 1:nc) {
    cuts <- c(0, cutscore[j], ((nk-1)*ni+1))
    categ <- cut(exp.TS, cuts, labels = FALSE, right=FALSE)
    bang <- ceiling(cuts)
    rec.s <- list(NA)
    for (i in 1:2) {
      rec.s[[i]] <- as.matrix(rec.mat[, (bang[i] + 1):bang[i + 1]])
    }
    for (i in 1:nn) {
      esacc[j, i] <- sum(rec.s[[categ[i]]][i, ])
    }
    escon[j, ] <- rowSums(rec.s[[1]])^2 + rowSums(rec.s[[2]])^2
  }
  if (nc > 1) {
    simul <- matrix(NA, nn, 2, dimnames = list(round(ut, 3), c("Accuracy", "Consistency")))
    cuts <- c(0, cutscore, ((nk-1)*ni+1))
    categ <- cut(exp.TS, cuts, labels = FALSE, right=FALSE)
    bang <- ceiling(cuts)
    rec.s <- list(NA)
    for (i in 1:(nc + 1)) {
      rec.s[[i]] <- as.matrix(rec.mat[, (cuts[i] + 1):cuts[i + 1]])
    }
    for (i in 1:nn) {
      simul[i, 1] <- sum(rec.s[[categ[i]]][i, ])
    }
    what <- matrix(0, nn, 1)
    for (i in 1:(nc + 1)) {
      what <- what + rowSums(rec.s[[i]])^2
    }
    simul[, 2] <- what
    ans <- (list(Marginal = rbind(cbind(Accuracy = rowMeans(esacc), 
                                        Consistency = rowMeans(escon)), Simultaneous = colMeans(simul)), 
                 Conditional = list(Accuracy = cbind(t(esacc), Simultaneous = simul[, 
                                                                                    1]), Consistency = cbind(t(escon), Simultaneous = simul[,2]))))
    ans
  }
  else ans <- (list(Marginal = cbind(Accuracy = rowMeans(esacc), 
                                     Consistency = rowMeans(escon)), Conditional = list(Accuracy = t(esacc), 
                                                                                        Consistency = t(escon))))
  ans
  
}