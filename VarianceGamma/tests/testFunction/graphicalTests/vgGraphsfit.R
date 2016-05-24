test.vgGraphfitvgFitStart <- function () {
  maxrows <- nrow(params)
  vgCResults <- matrix(nrow=maxrows,ncol=3,dimnames=list(NULL,c("true","SL",
    "MoM")))
  sigmaResults <- matrix(nrow=maxrows,ncol=3,dimnames=list(NULL,c("true","SL",
    "MoM")))
  thetaResults <- matrix(nrow=maxrows,ncol=3,dimnames=list(NULL,c("true","SL",
    "MoM")))
  nuResults <- matrix(nrow=maxrows,ncol=3,dimnames=list(NULL,c("true","SL",
    "MoM")))

  rownum <- 1
  pdf("testvgFitStart.pdf",height=7,width=11,paper="a4")
  for (i in 1:nrow(params)){
      param <- params[i,]
      data.vector <- rvg(n,param = param)  
      paramStartSL <- try(as.numeric(vgFitStart(data.vector,
        startValues="SL",startMethodSL="Nelder-Mead")$paramStart),
        silent=FALSE)
      if (class(paramStartSL) == "try-error") paramStartSL <- rep(NA,4)

      paramStartMoM <- try(as.numeric(vgFitStart(data.vector,
        startValues="MoM",startMethodMoM="Nelder-Mead")$paramStart),
        silent=FALSE)
      if (class(paramStartMoM) == "try-error") paramStartMoM <- rep(NA,4)

      startValues <- list(paramStartSL,paramStartMoM)

      vgCResults[rownum,1] <- param[1]
      sigmaResults[rownum,1] <- param[2]
      thetaResults[rownum,1] <- param[3]
      nuResults[rownum,1] <- param[4]
      for (meth in 2:3){
        vgCResults[rownum,meth] <- startValues[[meth-1]][1]
        sigmaResults[rownum,meth] <- exp(startValues[[meth-1]][2])
        thetaResults[rownum,meth] <- startValues[[meth-1]][3]
        nuResults[rownum,meth] <- exp(startValues[[meth-1]][4])
      }

      rownum <- rownum+1
    }


pairs(vgCResults,main="vgC")
pairs(sigmaResults,main="sigma")
pairs(thetaResults,main="theta")
pairs(nuResults,main="nu")
}

test.vgGraphfitvgFit <- function () {
  maxrows <- nrow(params)
  rownum <- 1
  pdf("testVgFit.pdf",height=7,width=11)
  for (i in 1:nrow(params)){
    param <- params[i,]
    data.vector <- rvg(n, param = param)
    results <- vgFit(data.vector, startValues = "SL", method = method,
      hessian = hessian, plots = FALSE, printOut = FALSE)
    paramFit <- as.numeric(results$param)
    par(mfrow = c(1,2), mar=c(5,4,5,2)+0.1)
    curve(dvg(x, param = param), col = "red", type = "n", xlab = "sample",
      range(data.vector)[1], range(data.vector)[2])
    hist(data.vector, freq = FALSE, breaks = 15, main = "", add = TRUE)
    mtext(expression(bold("Fitted and Actual")), line = 3.75, cex = 1.15)
    mtext(bquote(paste(vgC==.(param[1]),", ",
                       sigma==.(param[2]),
                       theta==.(param[3]),", ",
                       nu==.(param[4]), sep = "")), line = 2.5, cex = 1.15)
    curve(dvg(x, param = param), add = TRUE, col = "red",
      range(data.vector)[1], range(data.vector)[2])
    curve(dvg(x, param = paramFit), add = TRUE, col = "blue",
      range(data.vector)[1], range(data.vector)[2])
    logHist(data.vector, breaks = 15, main = "", xlab = "sample")
    mtext(expression(bold("Fitted and Actual")), line = 3.75, cex = 1.15)
    mtext(bquote(paste(vgC==.(param[1]),", ",
                       sigma==.(param[2]),
                       theta==.(param[3]),", ",
                       nu==.(param[4]), sep = "")), line = 2.5, cex = 1.15)
    curve(log(dvg(x, param = param)), add = TRUE, col = "blue",
      range(data.vector)[1], range(data.vector)[2])
    curve(log(dvg(x, param = paramFit)), add = TRUE, col = "red",
      range(data.vector)[1], range(data.vector)[2])
    rownum <- rownum+1
  }
}
