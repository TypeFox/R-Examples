Lee.pnr <- function(cutscore, pnr.out)	{
  ut <- pnr.out[[1]] #evalutation points (discrete total scores)
  we <- pnr.out[[2]] #density at each evaluation   
  Pij <- pnr.out[[3]] #probability array

  
  nn <- length(ut)
  nc <- length(cutscore)
  nk <- dim(Pij)[2]	
  ni <- dim(Pij)[3]
  highest <- (nk-1)*ni
  
  rec.mat <- gen.rec.raw(Pij, ut) 
  
  ut[ut<0] <- 1e-4
  ut[ut>highest] <- highest
  
  esacc <- escon <- matrix(NA, nc, nn, dimnames = list(paste("cut at", cutscore), round(ut, 3)))
  
  for (j in 1:nc) {
    cuts <- c(0, cutscore[j], highest + 1)
    categ <- cut(ut, cuts, labels = FALSE, right=FALSE)
    if(is.na(categ[1])) categ[1]<-1
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
    simul <- matrix(NA, nn, 2, dimnames = list(round(ut, 
                                                     3), c("Accuracy", "Consistency")))
    cuts <- c(0, cutscore, highest + 1)
    categ <- cut(ut, cuts, labels = FALSE, right=FALSE)
    bang <- ceiling(cuts)
    rec.s <- list(NA)
    for (i in 1:(nc + 1)) {
      rec.s[[i]] <- as.matrix(rec.mat[, (cuts[i] + 1):cuts[i+1]])
    }
    for (i in 1:nn) {
      simul[i, 1] <- sum(rec.s[[categ[i]]][i, ])
    }
    what <- matrix(0, nn, 1)
    for (i in 1:(nc + 1)) {
      what <- what + rowSums(rec.s[[i]])^2
    }
    simul[, 2] <- what
    ans <- list(Marginal = rbind(cbind(Accuracy = apply(esacc, 1, weighted.mean, we), 
                                        Consistency = apply(escon, 1, weighted.mean, we)), Simultaneous = apply(simul,2,weighted.mean, we)), 
                 Conditional = list(Accuracy = cbind(t(esacc), Simultaneous = simul[,1]), Consistency = cbind(t(escon), Simultaneous = simul[,2])))
    ans
  }
  else ans <- list(Marginal = rbind(cbind(Accuracy = apply(esacc, 1, weighted.mean, we), 
                                           Consistency = apply(escon, 1, weighted.mean, we))), 
                                     
                                     Conditional = list(Accuracy = t(esacc)), Consistency = t(escon))
  ans
}