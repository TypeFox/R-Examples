#####################################################################################
#
#
#   FUNCTIONS TO LOAD AND PREPARE DATA FOR THE FIT
#   See reference manual for description
#

############################################
#
#    GUI
#
#
runUI <- function()
    runApp(
        system.file('gui',                                                    
         package='BBEST'),
         launch.browser=TRUE)
         
         
####################################
# FUNCTIONS FOR SEPARATE DATA BANKS
#

read.sqa <- function(file = stop("'file' must be specified")){
  
  sqa <- scan(file=file, what="list", sep="\n")
  N <- length(sqa)
  
  i.start <- 0
  nBanks <- 0
  for(i in 1:N){
    if(strsplit(sqa[i], split=" ")[[1]][1]=="#L"){
	    i.start[nBanks+1] <- i+1 
      nBanks <- nBanks + 1 
    }
  }
  ids <- 0
  for(i in 1:nBanks)
    ids[i] <- strsplit(sqa[i.start[i]-4], split = " ")[[1]][4]  
  
  dat <- read.table(file=file, header=FALSE, col.names=c("x", "y", "e1", "e2", "e3"))
  res <- list()
    
  bank <- 1
  i.start <- 1
  for(i in 2:length(dat$x)){
    if(dat$x[i] < dat$x[i-1]){
  	  res[[bank]] <- list()
      class(res[[bank]]) <- "data"
      res[[bank]]$x <- dat$x[i.start:(i-1)]
      res[[bank]]$y <- dat$y[i.start:(i-1)]
      res[[bank]]$SB <- rep(0, length(res[[bank]]$x))
      res[[bank]]$id <- ids[bank]
      bank <- bank + 1 
      i.start <- i	
	  }
  }
  res[[bank]] <- list()
  class(res[[bank]]) <- "data"
  res[[bank]]$x <- dat$x[i.start:length(dat$x)]
  res[[bank]]$y <- dat$y[i.start:length(dat$x)]
  res[[bank]]$SB <- rep(0, length(res[[bank]]$x))
  res[[bank]]$id <- ids[bank]
  
  for(i in 1:bank)
    res[[bank]]$y[ which( is.na(res[[bank]]$y) ) ] <- 0
  
  return(res)
}


###
prepare.banks.data <- function(data, n.banks=4, lambda_1, lambda_2, lambda_0, x_1, x_2, 
                         n.atoms, scatter.length, ADP, n.regions){
  for(i in 1:n.banks){
    cat("\n\n==================================\n\n")
    cat("Preparing bank # ", i, "\n")
    cat("\n")
    data[[i]] <- set.sigma(data[[i]], n.regions=n.regions)
    data[[i]] <- set.lambda(data[[i]], lambda=NA, lambda_1, lambda_2, lambda_0, x_1, x_2)
    data[[i]] <- set.SB(data[[i]], SB=NA, n.atoms, scatter.length, ADP, fit=FALSE)  
  }
  data
}

###
write.fix <- function(fit.results, file = stop("'file' must be specified")){

  N <- length(fit.results)
  if(!is.null(fit.results$fit.details)){
    fit.results <- list(fit.results)
  	N <- 1
  }
  options(warn=-1)
  for(i in 1:N){
    if(i==1) apnd <- FALSE else apnd <- TRUE 
    write(c(paste("#S ",i," Correction File for Bank ",fit.results[[i]]$fit.details$id,sep=""), "#L Q MULT ADD"), file=file, append=apnd)
	res <- cbind(fit.results[[i]]$x, rep(1,length(fit.results[[i]]$x)), -fit.results[[i]]$curves$bkg)
    write.table(res, file=file, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")  
  }
  options(warn=0)
}



####################################
# GENERAL PURPOSE FUNCTIONS
#


read.sqb <- function(file = stop("'file' must be specified")){
  sqb <- scan(file=file, what="list", sep="\n")
  N <- length(sqb)
  
  for(i in 1:N){
    if(strsplit(sqb[i], split=" ")[[1]][1]=="#L"){
	  i.start <- i+1 
      break
	}
  }

# _!_ This is very poor implementation...
# _!_ but it works

  Tfile <- file()
  for(i in i.start:N)
    write(sqb[i], file=Tfile, append=TRUE)
  dat <- read.table(Tfile)
  unlink(Tfile)
  dat[which(is.na(dat[,2])),2] <- 0
  return(list(x=dat[,1], y=dat[,2], sigma=rep(0, length(dat[,1])), lambda=rep(0, length(dat[,1])), SB=rep(0, length(dat[,1])) ))
}
###
set.data <- function(x, y, sigma=NA, lambda=NA, SB=NA){
  data <- list()
  data$x <- x
  data$y <- y
  if(length(sigma)==1)
    data$sigma <- rep(sigma, length(x))
  else
    data$sigma <- sigma
	
  data$lambda <- lambda
  data$SB <- SB
  class(data) <- "data"	  
  return(data)  
}

###
read.data <- function(file = stop("'file' must be specified"), ...){
  data <- read.table(file=file, header=TRUE,...)
  if(is.null(data$sigma))
    data$sigma <- NA
  if(is.null(data$lambda))
    data$lambda <- NA
  if(is.null(data$SB))
    data$SB <- 0
  if(is.null(data$x))
    colnames(data)[1] <- "x"
  if(is.null(data$y))
    colnames(data)[2] <- "y"
  data$y[which(is.na(data$y))]   <- 0
  class(data) <- "data"	   
  return(data)  
}

###
trim.data <- function(data, x.min, x.max){
  ind.min <- which(abs(data$x-x.min)==min(abs(data$x-x.min)))
  ind.max <- which(abs(data$x-x.max)==min(abs(data$x-x.max)))
  cut <- ind.min:ind.max
  dat <- list()
  class(dat) <- "data"	  
 
  dat$x <- data$x[cut]  
  dat$y <- data$y[cut]  
  if(!is.null(data$SB)) dat$SB <- data$SB[cut]  else dat$SB <- rep(0, length(dat$x))
  if(!is.null(data$sigma)) dat$sigma <- data$sigma[cut]  else dat$sigma <- rep(NA, length(dat$x))
  if(!is.null(data$lambda)) dat$lambda <- data$lambda[cut]  else dat$lambda <- rep(NA, length(dat$x))
  if(!is.null(data$smoothed)) dat$smoothed <- data$smoothed[cut]  else dat$smoothed <- rep(NA, length(dat$x))
  if(!is.null(data$id)) dat$id <- data$id  
  	  
  return(dat)  
}


###
set.sigma <- function(data, sigma=NA, x.bkg.only=NA, n.regions=10, thresh.scale=1){
  y.smoothed <- NA 
  k <- thresh.scale

  if(is.na(sigma[1])){
    if(is.na(x.bkg.only[1])){
      if(length(k)==1)
        k <- rep(k, n.regions)
      n <- floor(length(data$x)/n.regions)
      x.bkg.i <- 1:n
      sigma <- 0
      y.smoothed <- 0
      for(i in 1:n.regions){
        cat("\n step ", i, " of ", n.regions, "\n\n")
        if(i==n.regions)
          x.bkg.i <- x.bkg.i[1]:length(data$x)
        y.sm <- wmtsa::wavShrink(data$y[x.bkg.i], wavelet="s8", shrink.fun="hard", thresh.fun="universal",  thresh.scale=k[i], xform="modwt")
        sig <- sqrt(mean((y.sm-data$y[x.bkg.i])^2))
        sigma <- c(sigma, rep(sig, length(x.bkg.i)))
        y.smoothed <- c(y.smoothed, y.sm)
          x.bkg.i <- x.bkg.i + n
      }
      y.smoothed <- y.smoothed[-1]
      sigma <- sigma[-1]	  
    }
    else{
      x.min.i <- which(abs(data$x-x.bkg.only[1])==min(abs(data$x-x.bkg.only[1])))
      x.max.i <- which(abs(data$x-x.bkg.only[2])==min(abs(data$x-x.bkg.only[2])))
      x.bkg.i <- x.min.i:x.max.i
      y.smoothed <- wmtsa::wavShrink(data$y[x.bkg.i], wavelet="s8", shrink.fun="hard", thresh.fun="universal",  thresh.scale=k, xform="modwt")
      sigma <- rep(sqrt(mean((y.smoothed-data$y[x.bkg.i])^2)), length(data$x))
      y.smoothed <- c(data$y[1:(x.min.i-1)], y.smoothed, data$y[(x.max.i+1):length(data$x)])
    }
  }
  else{
    if(length(sigma)==1) sigma <- rep(sigma, length(data$x))
  }

  data$sigma <- sigma
  data$smoothed <- y.smoothed
  return(data)
}


###
set.lambda <- function(data, lambda=NA, lambda_1=NA, lambda_2=NA, lambda_0=NA, x_1=NA, x_2=NA){
  if(is.na(lambda[1])){
    if(is.na(x_1)) x_1 <- min(data$x)
    if(is.na(x_2)) x_2 <- max(data$x)
  
    lambda <- rep(lambda_0, length(data$x))
 
    a0 <- (lambda_1 - lambda_2) / (x_1-x_2)
    b0 <- (lambda_1*x_2 - lambda_2*x_1) / (x_2-x_1)
    ind.max <- which(abs(data$x-x_2)==min(abs(data$x-x_2)))
    ind.min <- which(abs(data$x-x_1)==min(abs(data$x-x_1)))
    lambda[ind.min:ind.max] <- a0 * data$x[ind.min:ind.max] + b0
  }
  lambda[which(lambda<=0)] <- 1e-6
  data$lambda <- lambda
  return(data)
}

###
set.SB <- function(data, SB=NA, n.atoms=NA, scatter.length=NA, ADP=NA, fit=FALSE, oneADP=TRUE, ADP.lim=c(0, 0.05)){
  if(is.na(SB[1])){
    if(is.na(n.atoms) || is.na(scatter.length))
      stop("Please provide SB or parameters n.atoms and scatter.length\n")
    if(is.na(ADP) && !fit)
      stop("Please provide ADP or set fit=TRUE\n")
    if(fit==TRUE){
      data$fitADP <- list(n.atoms=n.atoms, scatter.length=scatter.length, oneADP=oneADP, ADP.lim=ADP.lim)
	  SB <- rep(0, length(data$x))
	}
	else{
	  data$fitADP <- NULL
	  if(length(ADP)==1) ADP <- rep(ADP, length(n.atoms))
      N_total <- sum(n.atoms)
      f.av2 <- (sum(n.atoms*scatter.length)/N_total)^2
      f2.av <- sum(n.atoms*scatter.length^2)/N_total
      expADP <- 0
      for(j in 1:length(data$x))
        expADP[j] <- sum(n.atoms*scatter.length^2*exp(-ADP*data$x[j]^2)/N_total/f2.av)
   
      L <- (f.av2-f2.av)/f.av2
      SB <- 1-expADP*(1-L)
	}
  }
  data$SB <- SB
  return(data)
}


###
set.Gr <- function(data, r1=seq(0, 1, 0.005), r2=NA, rho.0, 
                   type1="gaussianNoise", type2=NA, sigma.f=NA, l=NA){

  K.DI <- list()
  KG.inv <- matrix.FT1 <- matrix.FT2 <- sigma.r <- bkg.r <- ff <- D<- NA
  K.DI$inv <- K.DI$det <- NA
  
  if(is.na(type1))
    cat("No constraints on G(r) behaviour included. \n")
  else if(type1=="gaussianNoise"){
	# noise in r-space  
    matrix.FT1 <- sineFT.matrix(Q=data$x, r=r1)
    delta <- c(diff(data$x)[1], diff(data$x))
    sigma.r <- 0
    cat("Calculating r-space noise... \n")
    for(j in 1:length(r1)){
      sigma.r[j] <- sum((2/pi*delta*data$x*sin(data$x*r1[j])*data$sigma)^2)
      sigma.r[j] <- sqrt(sigma.r[j])
    }
    # avoid dividing by zero  
    if(sigma.r[1]==0)
      sigma.r[1] <- sigma.r[2]
    cat("Calculating FT of the experimental data... \n")
    bkg.r <- sineFT(f.Q=data$y-1, Q=data$x, r=r1) + 4 * pi * rho.0 * r1
                                        # SB should be excluded from bkg estimation 
										# by setting term SB! 
										# data$y=S(Q)=F(Q)+1!
  }
  else if(type1=="correlatedNoise"){
    matrix.FT1 <- sineFT.matrix(Q=data$x, r=r1)
	cat("Calculating noise covariance matrix in r-space... \n")
    KG <- noise.cov.matrix.r(r=r1, Q=data$x, sigma=data$sigma)
    diag(KG) <- diag(KG) + abs(min(eigen(KG)$values)) * 1e4   # avoid singularity
    KG.inv <- solve(KG)  
    cat("Calculating FT of the experimental data... \n")
    bkg.r <- sineFT(f.Q=data$y-1, Q=data$x, r=r1) + 4 * pi * rho.0 * r1
  }
  else
    stop("Wrong type of low-r Gr contribution to likelihood. Should be either 'gaussianNoise' or 'correlatedNoise'\n")

  if(is.na(type2))
    cat("No constraints on bkg(r) behaviour included. \n")
  else if(type2=="gaussianProcess"){
    matrix.FT2 <- sineFT.matrix(Q=data$x, r=r2)
    K <- covMatrixSE(x=r2, sig=sigma.f, l=l)
    ff <- K$factor
    K <- K$cov
    K.DI <- covMatrix.DI(K)
  }
  else if(type2=="secondDeriv"){
    matrix.FT2 <- sineFT.matrix(Q=data$x, r=r2)
    D <- DMatrix(knots.x=r2)$matrix
  }  
  else
    stop("Wrong type of low-r2 condition. Should be either 'gaussianProcess' or 'secondDeriv'\n")
  
  cat("...done! \n")  
  
  data$Gr <- list(type1=type1, type2=type2, sigma.r=sigma.r, bkg.r=bkg.r, 
    matrix.FT1=matrix.FT1, matrix.FT2=matrix.FT2, KG.inv=KG.inv, D=D, 
    covMatrix=list(inv=K.DI$inv, det=K.DI$det, factor=ff), rho.0=rho.0, r1=r1, r2=r2)

  return(data)			  
}

###
write.fit.results <- function(fit.results, file = stop("'file' must be specified")){

  x <- fit.results$x
  y <- fit.results$curves$y - fit.results$curves$bkg
  SB <- fit.results$curves$SB
  bkg <- fit.results$curves$bkg
  if(length(fit.results$uncrt)>1) 
    stdev <- fit.results$uncrt$stdev
  else
    stdev <- rep(NA, length(x))
  scale <- fit.results$scale
  f <- fit.results$fit.details$scatter.length
  N <- fit.results$fit.details$n.atoms
  ADP <- fit.results$ADP
  if(is.null(f)) f <- NA
  if(is.null(N)) N <- NA
  if(is.null(ADP)) ADP <- NA
  m <- cbind(f, N, ADP)
  exp.data <- (y-SB)/scale + bkg + SB
  res <- cbind(x, y, stdev, SB, bkg, exp.data)
  
  knots.x <- fit.results$knots$x
  knots.y <- fit.results$knots$y
  knots <- cbind(knots.x, knots.y)
  knots <- format(knots,digits=6)
  
  options(warn=-1)
  write(c("# scale factor:", scale), file=file, append=FALSE)
  cat("\n", file=file, append=TRUE)
  
  write(c("# Atomic Displacement Parameters:"), file=file, append=TRUE)
  write.table(m, file=file, append=TRUE, col.names=c("f","N","ADP"), row.names=FALSE, quote=FALSE)
  cat("\n", file=file, append=TRUE)
  
  write(c("# knots positions:"), file=file, append=TRUE)
  write.table(knots, file=file, append=TRUE, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  cat("\n", file=file, append=TRUE)
  
  cat("############################################################## \n", file=file, append=TRUE)
  cat("# fit results \n", file=file, append=TRUE)
  cat("# columns: x; (scaled) corrected y; standard deviation in y due to noise and bkg uncertainty; coherent baseline; estimated background; raw data \n", file=file, append=TRUE)
  write.table(res, file=file, append=TRUE, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  options(warn=0)
}


sqa.split <- function(file = stop("'file' must be specified")){
  sqa <- scan(file=file, what="list", sep="\n")
  N <- length(sqa)
  
  i.start <- 0
  nBanks <- 0
  for(i in 1:N){
    if(strsplit(sqa[i], split=" ")[[1]][1]=="#L"){
	    i.start[nBanks+1] <- i+1 
      nBanks <- nBanks + 1 
    }
  }
  i.start[nBanks+1] <- length(sqa)+5
  
  name <- 0
  for(i in 1:nBanks){
    name <- strsplit(file, '[.]')[[1]]
    name <- paste(name[-length(name)], collapse = '.')
    name <- paste(name, "_b", i, ".sqa", sep="") 
    writeLines(sqa[ (i.start[i]-4):(i.start[i+1]-5)], con = name, sep = "\n", useBytes = FALSE)  
  }
  
}


fix.merge <- function(outfile, infile1, infile2, ...){
  files <- list(infile1, infile2, ...)
  N <- length(files)
  
  file_tmp <- scan(file=infile1, what="list", sep="\n")
  name_str_arr <- strsplit(file_tmp[1], split = " ")[[1]]
  name_str_arr[2] <- 1
  file_tmp[1] <- paste(name_str_arr, collapse = ' ') 
  
  writeLines(file_tmp, con = outfile, sep = "\n", useBytes = FALSE) 
  for(i in 2:N){
    file_tmp <- scan(file=files[[i]], what="list", sep="\n")
    name_str_arr <- strsplit(file_tmp[1], split = " ")[[1]]
    name_str_arr[2] <- i
    file_tmp[1] <- paste(name_str_arr, collapse = ' ') 
    write(file_tmp, file = outfile, sep = "\n", append=TRUE) 
  }
}
