test.vgGraphdpqrddvg <- function () {
  pdf("testddvg.pdf",height=7,width=11)
  par(mfrow = c(1,1))
  par(oma = c(2,0,4,0))
  par(mfrow = c(2,1))
  for (i in 1:maxrows) {
    param <- params[i,]
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
    if (nu < 2){
      maxDens <- dvg(vgMode(param = param), param = param, log = FALSE)
      vgRange <- vgCalcRange(param = param, tol = 10^(-2)*maxDens, density = TRUE)
    } else {
      if (theta < 0) {
        vgRange <- c(vgC-2,vgC+6)
      } else {
        if (theta >0 ) {
          vgRange <- c(vgC-6,vgC+2)
        } else {
          vgRange <- c(vgC-4,vgC+4)
        }
      }
    }
    curve(dvg(x, param = param, log = FALSE),col="red",type="l",xlab="x",
      vgRange[1],vgRange[2])
    mtext(expression(bold("Density of Variance Gamma")),
      line=3.5,cex=1.15)
    mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")), line=2.25,cex=1.15)
    abline(v = vgMode(param = param))

    curve(ddvg(x, param = param),col="red",type="l",xlab="x",
      vgRange[1],vgRange[2])
    abline(h = 0)
    mtext(expression(bold("Derivative of Density of Variance Gamma")),
      line=3.5,cex=1.15)
    mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")), line=2.25,cex=1.15)
    abline(v = vgMode(param = param))
  }
}

test.vgGraphvgdpqrBreaks <- function () {
  pdf("testvgBreaks.pdf",height=7,width=11)
  par(mfrow = c(1,1))
  par(oma = c(2,0,4,0))
  par(mfrow = c(2,1))
  for (i in 1:maxrows) {
    param <- params[i,]
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]

    if (nu < 2){
      maxDens <- dvg(vgMode(param = param), param = param, log = FALSE)
      vgRange <- vgCalcRange(param = param, tol = 10^(-2)*maxDens, density = TRUE)
    } else {
      if (theta < 0) {
        vgRange <- c(vgC-2,vgC+6)
      } else {
        if (theta >0 ) {
          vgRange <- c(vgC-6,vgC+2)
        } else {
          vgRange <- c(vgC-4,vgC+4)
        }
      }
    }
    curve(dvg(x, param = param, log = FALSE),col="red",type="l",xlab="x",
      vgRange[1],vgRange[2])
    mtext(expression(bold("Density of Variance Gamma")),
      line=3.5,cex=1.15)
    mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")), line=2.25,cex=1.15)
    bks <- unlist(vgBreaks(param = param))
    ##print(bks)
    abline(v = bks)
    abline(v = vgMode(param = param),col="blue")
    ##print(vgMode(param))
  }
}


### Create plots for qvg and pvg with quantile and probability lines
test.vgGraphdpqrpvgqvg <- function () {
  pdf("pvgqvg.pdf",height=7,width=11)      
  par(oma=c(2, 0, 4, 0))  
  for (i in 1:maxrows){
    param <- params[i,]
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
    vgURange <- qvg (p = 1 - 10^(-3),param = param)
    vgLRange <- qvg (p = 10^(-3),param = param)
    vgRange <- c(vgLRange,vgURange)
    ps <- seq(vgRange[1],vgRange[2],length=length(qs))
    
    par(mfrow=c(1,2))
    ## plot cdf and quantile function with extra lines
    curve(pvg(q= x,param = param),col="red",type="l",xlab="x",n=200,
                vgRange[1],vgRange[2])
    abline(h=qs,col="grey")
    abline(v=qvg(p = qs,param = param,nInterpol=numInt),col="grey")
    mtext(expression(bold("pvg with Quantiles")),
      line=3.5,cex=1.15)
    mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")),line=2.25,cex=1.15)
    curve(qvg(p = x,param = param,nInterpol=numInt),col="red",type="l",xlab="q",
      ylab="qvg(q,param)",n=200,0,1)
    abline(h=ps,col="grey")
    abline(v=pvg(q= ps,param = param),col="grey")
    mtext(expression(bold("qvg with Probabilities")),line=3.5,cex=1.15)
    mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")),line=2.25,cex=1.15)
  }  
}

## P-P and Q-Q plots
test.vgGraphdpqrppvgqqvg <- function () {
  params <- unique(params)
  maxrows <- dim(params)[1]
  qs <- c(0.001,0.01,0.025,0.05,0.1,0.2,0.4,0.5,0.6,0.8,0.9,0.95,0.975,0.99,
    0.999)
  qnames <- paste("q",qs,sep="")
  diffqnames <- paste("dq",qs,sep="")
  pnames <- paste("p",1:length(qs),sep="")
  diffpnames <- paste("dp",1:length(qs),sep="")
  
  pdf("Graphs/ppqqtest.pdf",height=7,width=11)
  par(oma=c(2, 0, 4, 0))

  for (i in 1:maxrows){
    param <- params[i,]
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
    dataVector <- rvg(n,param = param)
    
    par(mfrow=c(1,2))     
    ppvg(y = dataVector,param = param,main="",pch=".")
    mtext(expression(bold("Variance Gamma P-P Plot")),
          line=3.5,cex=1.15)
    mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")),line=2.25,cex=1.15)
    qqvg(y = dataVector,param = param,main="",pch=".")
    mtext(expression(bold("Variance Gamma Q-Q Plot")),
          line=3.5,cex=1.15)
    mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")),line=2.25,cex=1.15)
  }
}

test.vgGraphdpqrrvg <- function () {
  pdf("testrvghist.pdf",height=7, width=11)
  par(mfrow=c(1,2))
  par(oma=c(2, 0, 4, 0))
  for (i in 1:maxrows){
    param <- params[i,]
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
  par(mfrow=c(1,2))
  dataVector <- rvg(n, param = param)
  summary(dataVector)
  curve(dvg(x, param = param, log = FALSE), col="red", type="n", xlab="sample",
    ylab="Density", range(dataVector)[1], range(dataVector)[2])
  hist(dataVector, freq=FALSE, breaks=20 ,main="", add=TRUE)
  mtext(expression(bold("Test of rvg")), line=3.5, cex=1.15)
  mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")), line=2.25, cex=1.15)
  
  curve(dvg(x, param = param, log = FALSE), add=TRUE, col="red", 
    range(dataVector)[1], range(dataVector)[2])
  logHist(dataVector, breaks=20, main="", xlab="sample")
  mtext(expression(bold("Test of rvg")), line=3.5, cex=1.15)
  mtext(bquote(paste(vgC==.(vgC),", ", sigma==.(sigma),", ",
      theta==.(theta),", ", nu==.(nu),sep="")), line=2.25, cex=1.15)
  curve(log(dvg(x, param = param, log = FALSE)), add=TRUE, col="red", 
    range(dataVector)[1], range(dataVector)[2])
  }
}




