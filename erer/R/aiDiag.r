aiDiag <- function(x, digits = 3, ...){
  if ((!inherits(x, "aiStaFit")) & (!inherits(x, "aiDynFit"))) {
    stop("Please provide an object from 'aiStaFit' or 'aiDynFit'.\n")
  }
  if (inherits(x, "aiStaFit")) dat <- x$y
  if (inherits(x, "aiDynFit")) dat <- x$daDyn     
  formu  <- x$formula
  nShare <- x$nShare   
  if (!is.null(x$est$rho)) {dat <- x$est$data; formu <- x$est$formula}
  
  test.BG <- matrix(0, nrow=nShare-1, ncol=2)
  for (i in 1:(nShare-1)) {
      BG <- bgtest(formu[[i]], order=1, order.by=NULL, 
        type=c("Chisq"), data=dat)
      test.BG[i,1] <- round(as.numeric(BG["statistic"]), digits)
      test.BG[i,2] <- round(as.numeric(BG["p.value"]), digits-1)
  }
  test.BP <- matrix(0, nrow=nShare-1, ncol=2)
  for (i in 1:(nShare-1)) {
      BP <- bptest(formu[[i]], varformula=NULL, studentize=TRUE, data=dat)
      test.BP[i,1] <- round(as.numeric(BP["statistic"]), digits)
      test.BP[i,2] <- round(as.numeric(BP["p.value"]), digits-1)
  }
  test.RESET <- matrix(0, nrow=nShare-1, ncol=2)
  for (i in 1:(nShare-1)) {
      RESET <- resettest(formu[[i]], power=2:3, type=c("fitted"), data=dat)
      test.RESET[i,1] <- round(as.numeric(RESET["statistic"]), digits)
      test.RESET[i,2] <- round(as.numeric(RESET["p.value"]), digits-1)
  }
  test.JB <- matrix(0, nrow=nShare-1, ncol=2)
  for (i in 1:(nShare-1)) {
      JB <- jarque.bera.test(na.omit(residuals(x$est$eq[[i]])))
      test.JB[i,1] <- round(as.numeric(JB["statistic"]), digits)
      test.JB[i,2] <- round(as.numeric(JB["p.value"]), digits-1)
  }

  test <- cbind(test.BG, test.BP, test.RESET, test.JB)
  colnames(test) <- c("BG.stat", "BG.p", "BP.stat", "BP.p", 
      "RESET.stat", "RESET.p", "JB.stat", "JB.p")
  dia <- cbind(Equation = x$share[-x$nOmit], format(test, trim = TRUE))
  dia <- data.frame(dia)
  return(dia)
}