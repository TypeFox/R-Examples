combos <- function(n)
  {
if(n < 2)
  {
  cat("\nn must be greater than 1.\n\n")
   }
else
  if(n > 11)
    {
  cat("\n CAUTION! Output size increases exponentially with n. \n")
  cat(" Result for n = ", n, "will be a ", 2^n-1, "*", n, " matrix.\n")
  cat(" Do you really want to proceed?\n")
  choice <- menu(c("Proceed","Quit"))
  if(choice == 1)
  combos1(n)
  else
  {}
}
else
  combos1(n)
}

combos1 <- function(n)
  {
  require(gtools)
  x <- cbind(combinations(n,1,1:n),array(0,dim=c(n,n-1)))
  for(i in 2:n)
   {
#    nc <- factorial(n)/(factorial(i)*factorial(n-i))
     nc <- dim(combinations(n,i,1:n))[1]
    x <- rbind(x,cbind(combinations(n,i,1:n),array(0,dim=c(nc,n-i))))
    }
  len <- dim(x)[1]
  x.index <- cbind(as.vector(1:len),as.vector(x))
  x.index <- cbind(x.index[,1][x.index[,2]>0],x.index[,2][x.index[,2]>0])
  x.bin <- array(0,dim=c(len,n))
  x.bin[x.index] <- 1
  list(ragged = x, binary = x.bin)
}
  
current.model <- function(y, current.comb, xcan, family = "gaussian", gof = "RMSPE")
  {
  comb.data <- data.frame(xcan[,current.comb])
  colnames(comb.data) <- colnames(xcan)[current.comb] 
  data <- data.frame(y,comb.data)
  depv <- names(data)[1]
  n.comb <- dim(comb.data)[2]
  xs <- vector("character",n.comb)
  for (i in 1:(n.comb - 1))
     xs[i] <- paste(names(comb.data)[i], "+", sep = "")
  xs[n.comb] <- names(comb.data)[n.comb]
  xss <- paste(xs, collapse = " ", sep = "")
  formu <- formula(paste(depv, "~", xss, sep = ""))
    if(gof == "RMSPE")
    gf <- sqrt(sum((glm(formu, data = data, family = family)$fitted.values - y)^2))
    if(gof == "logLik")
    gf <- as.vector(logLik(glm(formu, data = data, family = family)))
    if(gof == "Rsqu")
    gf <- summary(lm(formu, data = data))$r.squared
  gf
}

all.regs <- function(y, xcan, family = "gaussian", gof = "RMSPE",
                     print.vars = FALSE)
  {
  if(!is.vector(y) && dim(y)[2] != 1)
   cat("\ny must be a vector or a single column data frame")
  pcan <- dim(xcan)[2]
  n <- (2^pcan)-1
  combs <- combos1(pcan)$ragged
  if(gof != "RMSPE" && gof != "logLik" && gof != "Rsqu")
    {
     cat("\n gof (goodness of fit measure) must equal")
     cat("\n \"RMSPE\" (Root-mean-square \"prediction\" error")
     cat("\n \"logLik\" (Log-Likelihood) or")
     cat("\n \"Rsqu\" (R-squared)\n\n")
     }
    else
      {
    if(gof == "RMSPE")
      gfs <- sqrt(sum((glm(y ~ 1, family = family)$fitted.values - y)^2))
    if(gof == "logLik")
      gfs <- as.vector(logLik(glm(y ~ 1, family = family)))
    if(gof == "Rsqu")
      gfs <- 0
       }
  for(i in 1:n){
    if(i %% 500 == 0)
      cat(i,"regressions calculated:",n-i,"to go...\n")
    current.comb <- as.vector(combs[i,][combs[i,]>0])
    combn <- paste(names(data.frame(xcan)[current.comb]), "", collapse="")
    if(gof == "RMSPE")
    new.line <- current.model(y, current.comb, xcan,
                 family = family, gof = "RMSPE")
    if(gof == "logLik")
    new.line <- current.model(y, current.comb, xcan,
                 family = family, gof = "logLik")
    if(gof == "Rsqu")
    new.line <- current.model(y, current.comb,xcan, gof = "Rsqu")
    gfs <- c(gfs, new.line)
             }
  if(print.vars)
    {
    cat("regressions done: formatting results\n")
    var.names <- "Theta"
    for(i in 1:n)
      {
    current.comb <- as.vector(combs[i,][combs[i,]>0])
    combn <- paste(names(data.frame(xcan)[current.comb]), "", collapse="")
    new.line <- combn
    var.names <- c(var.names, new.line)
       }
    gfs <- data.frame("variable combination" = var.names, gof = gfs)
      }
    gfs
  }

partition <- function(gfs, pcan, var.names = NULL)
  {
      if(pcan >12)
        stop("Number of variables must be < 13 for current implementation",
             call. = FALSE)
       else 
       if(pcan >9)
           warning("hier.part produces a rounding error if number of variables >9
See documentation.",
                   call. = FALSE)
        {
      n <- 2^pcan
      if((is.vector(gfs) && length(gfs) != n) || (!is.vector(gfs) && dim(gfs)[1] != n))
        {
        cat("\nIncorrect number of goodness of fit measures.\n")
        cat("First element must be null model, last full model\n")
        cat("Total number of gof measures should = 2^pcan\n\n")
      }
      else
      if(is.vector(gfs))
        {
      theta <- gfs[1]
      fin <- gfs[2:n]
          }
      else
        {
      wgfs <- dim(gfs)[2]
#gfs should be a vector or an array with goodness of fit
#measures in last col.  wgfs selects last column
      theta <- gfs[1,wgfs]
      fin <- gfs[2:n,wgfs]
          }
      len <- length(fin)
      IJ <- vector("numeric", pcan * 2)
      storage.mode(pcan) <- "integer"
      storage.mode(len) <- "integer"
      storage.mode(theta) <- "double"       
      storage.mode(fin) <- "double"
      storage.mode(IJ) <- "double"
      IJ <- .C("hierpart",
                pcan,
                len,
                theta,
                fin,
                IJ = IJ,
                PACKAGE = "hier.part")$IJ
      IJ <- array(IJ, dim = c(pcan,2))
      IJ <- data.frame(t(data.frame(t(IJ), row.names = c("I","J"))),
                       row.names = var.names)
      IJ.perc <- IJ*100/sum(IJ)
      I <- data.frame(I=IJ[,1],row.names=var.names)
      I.perc <- I*100/sum(I)
      IJ <- cbind(IJ, "Total" = IJ$I +IJ$J)
  list(gfs = gfs, IJ = IJ, I.perc = I.perc)
     }
    }

hier.part <- function(y, xcan, family = "gaussian", gof = "RMSPE",
                      barplot = TRUE)
  {
    pcan <- dim(xcan)[2]
    if(pcan >12)
       stop("Number of variables must be < 13 for current implementation",
             call. = FALSE)
    else
      {
    gfs <- all.regs(y, xcan, family = family, gof = gof)
    HP <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
    if(barplot)
      {
       ymin <- min(c(0, floor(min(HP$I.perc)*0.1)*10))
       ymax <- ceiling(max(HP$I.perc)*0.1)*10
       barplot(t(HP$I.perc), col = c(1), ylim = c(ymin, ymax),
               ylab="% Independent effects (%I)")
      }
    list(gfs = gfs, IJ = HP$IJ, I.perc = HP$I.perc)
     }
    }         

rand.hp <- function(y, xcan, family = "gaussian", gof = "RMSPE",
                    num.reps = 100)
  {
    cat("\nPlease wait: running", num.reps, "randomizations \n")
    IJ <- hier.part(y, xcan, family = family, gof = gof,
                   barplot = FALSE)$IJ
    var.names = row.names(IJ)
    I <- data.frame(Obs = IJ[,1],row.names = var.names)
    I <- t(I)
    nvars <- dim(xcan)[2]
    npoints <- dim(xcan)[1]
    nreps <- num.reps + 1
    for(i in 1:num.reps)
      {
       for(j in 1:nvars)
        {
          o <- order(runif(npoints))
          xcan[,j] <-  xcan[,j][o]
         }
    Itemp <- hier.part(y, xcan, family = family, gof = gof,
                       barplot = FALSE)$IJ[,1]
    I <- rbind(I, t(Itemp))
     }
#RMSPE is smaller for better fits: logLik and Rsqu the opposite    
    if(gof == "RMSPE")
      I <- I * -1
    Z <- (I[1,1]-mean(I[2:nreps,1]))/sd(I[2:nreps,1])
    for(i in 2:nvars)
    Z <- c(Z,(I[1,i]-mean(I[2:nreps,i]))/sd(I[2:nreps,i]))
    sig <- vector("character", nvars)
    sig[Z >= rep(1.65,nvars)] <- "*"
    Iprobs <- data.frame(Obs = round(I[1,],2), Z.score = round(Z,2), sig95 = sig, row.names = var.names)
    if(gof == "RMSPE")
      I <- I * -1    
    list(Irands = I, Iprobs = Iprobs)
  }





