varlog.lam <-
function(sight1,sight2){
    if(sight1$call$form != sight2$call$form){
      stop("Need same sightability model to calculate the covariance")
    }  
    if(is.null(sight1$sight$note)){
      varbeta <- vcov(sight1$sight)
      beta <- coef(sight1$sight)
    }else{
      varbeta <- sight1$sight$varbet
      beta <- sight1$sight$bet
    }  
    y1 <- sight1$odat$total
    y2 <- sight2$odat$total
    n1 <- nrow(sight1$odat)
    n2 <- nrow(sight2$odat)
    inv.srate1 <- 1/sight1$odat$samp.rates
    inv.srate2 <- 1/sight2$odat$samp.rates
    fo <- sight1$call$form
    class(fo) <- "formula"
    tempnm1 <- terms(fo, data = sight1$odat)
    tempnm2 <- attr(tempnm1, "term.labels")
    covars1 <- sight1$odat[, tempnm2]  
    covars2 <- sight2$odat[, tempnm2]  

  # Do as much of the matrix multiplication outside of loop as possible
    xdat1 <- as.matrix(cbind(rep(1, n1), covars1))
    xdat2 <- as.matrix(cbind(rep(1, n2), covars2))
    xb1 <- xdat1%*%beta  # X'beta
    xb2 <- xdat2%*%beta  # X'beta
    xbb <- kronecker(xb1, t(xb2), FUN = "+")  # X1+x2
    xvarbeta <- xdat1%*%varbeta%*%t(xdat2)  # X Sig X
    smat <- matrix(0, n1, n2) # holder for some of the terms in the expression for smat
    for(i in 1:n1){
      for(j in 1:n2){
         xtemp1 <- as.vector(xdat1[i, ], mode = "numeric")
         xtemp2 <- as.vector(xdat2[j, ], mode = "numeric")
         xtot <- t(xtemp1+xtemp2)
         smat[i,j] <- (xtot%*%varbeta%*%t(xtot))/2 
       }
     }    
  
    smat.cov <- exp(-xbb-smat)*(exp(xvarbeta)-1)
    y.p1 <- as.matrix(y1*inv.srate1, n1, 1)
    y.p2 <- as.matrix(y2*inv.srate2, n2, 1)
 
    cov.total <- t(y.p1)%*%smat.cov%*%y.p2
  
    var.tau1.tau2 <- matrix(c(sight2$est[2], cov.total, cov.total, sight1$est[2]), ncol = 2, byrow = TRUE) 
    dfs <- matrix(c(1/sight2$est[1], -1/sight1$est[1]), 1, 2)
    varloglam <- list(loglamda = log(sight2$est[1]/sight1$est[1]), varloglamda = dfs%*%var.tau1.tau2%*%t(dfs))
    return(varloglam = varloglam)
}
