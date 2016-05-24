mlgarchSim <-
function(n, constant=c(0,0), arch=diag(c(0.1,0.05)),
  garch=diag(c(0.7,0.8)), xreg=NULL,
  backcast.values=list(lnsigma2=NULL, lnz2=NULL, xreg=NULL),
  innovations=NULL, innovations.vcov=diag(rep(1,length(constant))),
  check.stability=TRUE, verbose=FALSE)
{
  #check/change arguments:
  if(is.null(constant)){ constant <- 0 }
  if(is.null(arch)){ arch <- 0 }
  if(is.null(garch)){ garch <- 0 }
  c.code <- FALSE #necessary until c.code implemented

  #prepare:
  npluss1 <- n+1
  constant <- as.vector(constant)
  iRows <- NROW(constant)

  #xreg:
  if(is.null(xreg)){
    constantx <- matrix(constant,iRows,npluss1)
    Econstantx <- constant
  }else{
    if(!is.matrix(xreg)) stop("xreg is not a matrix")
    if(ncol(xreg)!=iRows) stop("ncol(xreg) must equal length(constant)")
    if(nrow(xreg)!=n) stop("nrow(xreg) must equal n")
    constantx <- constant + t(xreg)
    Econstantx <- as.numeric(rowMeans(constantx))
    constantx <- cbind(Econstantx, constantx) #backcast values (1 period)
  }

  #arch:
  if(length(arch) == 1 && arch==0){
    #arch <- list()
    #arch[[1]] <- matrix(0,iRows,iRows)
    arch <- matrix(0,iRows,iRows)
  }

  #garch:
  if(length(garch) == 1 && garch==0){
    garch <- matrix(0,iRows,iRows)
  }

  #make phi:
  phi <- arch + garch

  #check stability:
  if(check.stability){
    if(any(abs(eigen(phi)$values) >= 1)){
      mssg <- paste("The model may not be stable (one or more AR-roots is on or inside the unit circle)")
      print(mssg)
    }
  }

  #z series:
  if(is.null(innovations)){
    #for the future: check if diag(innovations.vcov)==1?
    innovations <- rmnorm(n, mean=rep(0,iRows),
      vcov=innovations.vcov)
  }
  mlnz2 <- log(innovations^2)
  Elnz2 <- colMeans(mlnz2)
  if(is.null(backcast.values$lnz2)){
    #to do: account for no lags
    mlnz2 <- rbind(as.numeric(Elnz2),mlnz2)
  }else{
    mlnz2 <- rbind(as.numeric(backcast.values$lnz2),mlnz2)
  }

  #make lnsigma2:
  lnsigma2 <- matrix(NA,npluss1,iRows)
  mI <- diag(rep(1,iRows))
  innovInit <- as.matrix(Econstantx) + arch %*% cbind(as.numeric(Elnz2))
  lnsigma2[1,] <- solve(mI-phi) %*% innovInit

  #recursion:
  if(c.code){
    #not available yet
  }else{
    for(i in 2:npluss1){
      #future code:
      #phiSum <- for(..) phi[[..]] etc.
      #archSum <- for(..) arch[[..]] etc.
      #lnsigma2[,] <- phiSum + archSum
      lnsigma2[i,] <- constantx[,i] + phi %*% lnsigma2[i-1,] + arch %*% mlnz2[i-1,]
    } #end for loop
  } #end if(c.code)

  #out:
  sigma <- exp(lnsigma2[-1,]/2)
  y <- sigma*innovations
  if(verbose){
    y <- cbind(y, sigma, innovations)
    colnames(y) <- c(paste("y",1:iRows,sep=""),
      paste("sd",1:iRows,sep=""),
      paste("innov",1:iRows,sep="") )
  }else{
    colnames(y) <- paste("y",1:iRows,sep="")
  }
  as.zoo(y)
}
