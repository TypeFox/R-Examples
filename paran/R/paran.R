`paran` <-
function(x=NA, iterations=0, centile=0, quietly=FALSE, status=TRUE, all=FALSE, cfa=FALSE, graph=FALSE, color=TRUE, col=c("black","red","blue"), lty=c(1,2,3), lwd=1, legend=TRUE, file="", width=640, height=640, grdevice="png", seed=0, mat=NA, n=NA) {

library(MASS)

# Confirm that x and mat were NOT both provided
  if ( !is.na(mat[[1]][1]) & !is.na(x[[1]][1]) ) {
    stop("\nYou must supply either x or mat but not both.")
    }

# Set number of variables
  if ( is.na(mat[[1]][1]) & !is.na(x[[1]][1]) ) {
    P <- length(as.matrix(x)[1,])    
    }
  if ( !is.na(mat[[1]][1]) & is.na(x[[1]][1]) ) {
    P <- length(mat[1,])    
    }

# Confirm correlation matrix, and not covariance matrix
  if (!is.na(mat[[1]][1])) {
    if (length(mat[1,]) != length(mat[,1])) {
      stop("\nThe matrix provided with the mat argument is not a correlation matrix.\n")
      }
    if ( is.element(FALSE , (diag(mat) == rep(1,P))) ) {
      stop("\nThe matrix provided with the mat argument is not a correlation matrix.\nParallel analysis is not compatible with the eigendecomposition of a covariance matrix.")
      }
    }

# Make the correlation matrix R from the data, or take
# from the supplied matrix
  if ( !is.na(x[[1]][1]) ) {
    N <- length(as.matrix(x)[,1])
    if ( !is.na(n) ) {
      warning("\nThe n argument is only for use with the mat argument. Ignoring n argument.")
      }
    R <- cor(x)
    cat("\nUsing eigendecomposition of correlation matrix.")
    }
  if ( !is.na(mat[[1]][1]) ) {
    if (!is.na(mat[[1]][1]) & is.na(n)) {
      stop("\nYou must also provide the sample size when using the matrix argument.")
      }
    N <- n
    R <- mat
    cat("\nUsing eigendecomposition of provided correlation matrix.")
    }

# quick validation of centile as an integer value
  centile <- round(centile)
  if (centile > 99 | centile < 0) {
    stop("\nYou must specify a centile value between 1 and 99.\n(Specifying centile 0 will use the mean.)")
    }

# Perform pca or cfa
  if (cfa == FALSE) {
    eigenvalues <- eigen(R, only.values = TRUE, EISPACK = FALSE)[[1]]
    }
  if (cfa == TRUE) {
    C <- R - ginv(diag(diag(ginv(R))))
    eigenvalues <- eigen(C, only.values = TRUE, EISPACK = FALSE)[[1]]
    }
  
# Get the eigenvalues .  .  .
  Ev <- eigenvalues

# note which model
  model <- "component"
  models <- "components"
  Model <- "Component"
  Models <- "Components"
  if (cfa == TRUE) {
    model <- "factor"
    models <- "factors   "
    Model <- "Factor   "
    Models <- "Factors"
    }

# clean up iteration and determine value
   if (iterations<1) {
    iterations <- 30*P
    }
   if (iterations<0) {
    cat("\nInvalid number of iterations! Using default value of ",iterations,"\n",sep="")
    }

# prepare to save the results of each pca
#    N <- length(as.matrix(x[1]))
    if ( cfa == FALSE ) {
      SimEvs <- matrix(NA,iterations,P)
      }
    if ( cfa == TRUE ) {
      SimEvs <- matrix(NA,iterations,P)
      }

# Let the user know the program is working if neccesary
  if (status==TRUE) {
    if (iterations >= 10) {
      cat("\nComputing: ")
      }
    }

  for (k in 1:iterations) {
   
# Yet _more_ letting the user know the program is working!
    if (status == TRUE) {
      if (k %% (iterations/10) == 1 & iterations >= 10 & k > 1) {
        pct <- (k%/%(iterations/10))*10
        cat(pct,"%  ",sep="")
        }
      if (k == iterations) {
        cat("100%\n")
        }
      }

# initialize previously created random dataset.
    Sim <- matrix(NA,N,P)
      
# Create the random dataset.
    # for normally distributed simulations
    if (seed != 0) {
      set.seed(seed*k)
      }
    Sim <- matrix(rnorm(N*P),N,P)

# Extract principal components or factors from the random dataset
# (which is the same size and dimension as the user dataset.)

    if (cfa == FALSE) {
      eigenvalues <- eigen(cor(Sim), only.values = TRUE, EISPACK = FALSE)[[1]]        
      }
    if (cfa == TRUE) {
      R <- cor(Sim)
      C <- R - ginv(diag(diag(ginv(R))))
      eigenvalues <- eigen(C, only.values = TRUE, EISPACK = FALSE)[[1]]
      }

# Get the eigenvalues .  .  .
    Evs <- eigenvalues

# Save eigenvalues
    SimEvs[k,] <- Evs

# end the for k loop
  }

# display if neccesary
  if (quietly == TRUE) {
    cat("\n")
    }
  if (quietly == FALSE) {
  
    cat("\n\nResults of Horn's Parallel Analysis for ",model," retention\n",sep="")

    if (iterations == 1) {
      if (centile == 0) {
        cat("1 iteration, using the mean estimate","\n",sep="")
        }
      if (centile != 0) {
        cat("1 iteration, using the ",centile," centile estimate","\n",sep="")
        }
      }

    if (iterations > 1) {
      if (centile == 0) {
        cat(iterations," iterations, using the mean estimate","\n",sep="")
        }
      if (centile != 0 & centile != 50) {
        cat(iterations," iterations, using the ",centile," centile estimate","\n",sep="")
        }
      if (centile == 50) {
        cat(iterations," iterations, using the ",centile," centile (median) estimate","\n",sep="")
        }    
      }

    cat("\n--------------------------------------------------","\n")
    cat(Model,"  Adjusted    Unadjusted    Estimated","\n")
    cat("            Eigenvalue  Eigenvalue    Bias","\n")
    cat("--------------------------------------------------","\n")
    }

  RndEv = c(1:P)*NA 

  if (centile > 0) {
    for (p in 1:P) {
      RndEv[[p]] <- quantile(SimEvs[,p],probs=centile/100)[[1]]
      }
    }
  if (centile==0) {
    for (p in 1:P) {
      RndEv[[p]] <- mean(SimEvs[,p])      }
    }

  if (Ev[[1]] < 1 | RndEv[[1]] < 1) { 
    if (quietly == FALSE) {
      cat("No components passed.","\n")
      cat("--------------------------------------------------","\n")
      stop
      }
    }

  Bias <- rep(0,P)
  AdjEv <- rep(1,P)
  for (p in 1:P) {
    if (cfa == TRUE) {
      Bias[p] <- RndEv[p]
      }
    if (cfa == FALSE) {
      Bias[p] <- RndEv[p] - 1
      }
    AdjEv[p] <- Ev[[p]] - Bias[p]
    }

  # calculate how many components or factors to return by counting those 
  # components or factors with adjusted eigenvalues greater than one for
  # PCA, or greater than zero for CFA.
  y <- NA
  for (x in 1:P) {
    y <- x
    if (cfa == TRUE) {
      if (AdjEv[x] <= 0) {
        y <- x - 1
        retained <- y
        break
        }
      }
    if (cfa == FALSE) {
      if (AdjEv[x] <= 1) {
        y <- x - 1
        retained <- y
        break
        }
      }
    }

  if ( all == TRUE ) {
    y <- P
    }

  for (x in 1:y) {
    if ( AdjEv[x] >=0 ) {
      AdjSpace = " "
      }
    if ( AdjEv[x] < 0 ) {
      AdjSpace = ""
      }
    if ( Ev[[x]] >= 0 ) {
      EvSpace = " "
      }
    if ( Ev[[x]] < 0 ) {
      EvSpace = ""
      }
    if ( Bias[x] >= 0 ) {
      BiasSpace = " "
      }
    if ( Bias[x] < 0 ) {
      BiasSpace = ""
      }

# Pad the rear of x in case of single-digits
    if ( x > 9 ) {
      xPad = ""
      }
    if ( x <= 9 ) {
      xPad = " "
      }

# Pad the front of AdjEv in case of eigenvalues > 10, 100, etc.
    AdjFPad = "   "
    if ( round(AdjEv[x]) >= 10 ) {
      AdjFPad = "  "
      }
    if ( round(AdjEv[x]) >= 100 ) {
      AdjFPad <- " "
      }

# Set the strtrim number SN
    SN <- 8
    if ( abs(AdjEv[x]) >= 10 ) {
      SN <- 9
      }
    if ( abs(AdjEv[x]) >= 100 ) {
      SN >= 10
      }
    if ( AdjEv[x] < 0 ) {
      SN <- SN + 1
      }

# Pad the front of Ev in case of eigenvalues > 10, 100, etc.
    EvFPad = "   "
    if ( round(Ev[[x]]) >= 10 ) {
      EvFPad = "  "
      }
    if ( round(Ev[[x]]) >= 100 ) {
      EvFPad = " "
      }

# Set the strtrim number SN
    EvSN <- 8
    if ( abs(Ev[[x]]) >= 10 ) {
      EvSN <- 9
      }
    if ( abs(Ev[[x]]) >= 100 ) {
      EvSN <- 10
      }
    if (abs(Ev[[x]]) >= .0000005) {
      EvZPad <- ""
      }
    if (abs(Ev[[x]]) < .0000005) {
      Ev[[x]] <- 0
      EvZPad <- ".000000"
      }

# Set the strtrim number SN
    BiasSN <- 8
    if ( Bias[x] >= 10 ) {
      BiasSN <- 9
      }
    if ( Bias[x] >= 100 ) {
      BiasSN >= 10
      }

    if (quietly == FALSE) {
      cat(x,xPad,"      ",AdjFPad,AdjSpace,strtrim(AdjEv[x],SN),EvFPad,EvSpace,strtrim(Ev[[x]],EvSN),EvZPad,"     ",BiasSpace,strtrim(Bias[x],BiasSN),"\n", sep="")
      }
    }
  if (quietly == FALSE) {
    cat("--------------------------------------------------","\n")
    if (cfa == TRUE) {
      cat("\nAdjusted eigenvalues > 0 indicate dimensions to retain.\n(",retained," ",models," retained)\n\n",sep="")
      }
    if (cfa == FALSE) {
      cat("\nAdjusted eigenvalues > 1 indicate dimensions to retain.\n(",retained," ",models," retained)\n\n",sep="")
      }
    }

# Graph it if needed
  if (graph == TRUE) {
    AdjEvCol = col[1]
    EvCol = col[2]
    RndEvCol = col[3]
    AdjEvLty = 1
    EvLty = 1
    RndEvLty = 1
    if (color == FALSE) {
      EvCol = "black"
      RndEvCol = "black"
      EvLty = lty[2]
      RndEvLty = lty[3]
      }
    if (cfa==FALSE) {
      par(yaxs='i', xaxs='i', lab=c(P,ceiling(max(AdjEv[1],Ev[1],RndEv[1])),2))
      plot.default(c(1:P), RndEv, type='o', main='Parallel Analysis', xlab='Components', ylab='Eigenvalues', pch=20, col=RndEvCol, lty=RndEvLty, lwd=lwd, xlim=c(.5,P+.5), ylim=c(min(AdjEv, Ev,RndEv)-.5,ceiling(max(AdjEv[[1]],Ev[[1]],RndEv[[1]]))))
      }
    if (cfa==TRUE) {
      par(xaxp=c(1,P,1))
      plot.default(c(1:P), RndEv, type='o', main='Parallel Analysis', xlab='Factors', ylab='Eigenvalues', pch=20, col=RndEvCol, lty=RndEvLty, lwd=lwd, xlim=c(.5,P+.5), ylim=c(min(AdjEv, Ev,RndEv)-.5,ceiling(max(AdjEv[[1]],Ev[[1]],RndEv[[1]]))))
      }
    if (cfa == TRUE) {
      abline(h=0, col='grey', lwd=.5)
      }
    if (cfa == FALSE) {
      abline(h=1, col='grey', lwd=.5)
      }
    points(c(1:P),AdjEv, type='o', col=AdjEvCol, lty=AdjEvLty, pch=21, bg='white', lwd=lwd)
    points(c(1:P),Ev, type='o', col=EvCol, lty=EvLty, pch=20, lwd=lwd)
    if (retained >= 1) {
      points(c(1:retained), AdjEv[1:retained], type='p', pch=19, col=AdjEvCol, lty=AdjEvLty, lwd=lwd)
      }

# Add a legend to help with interpretation (thanks to Ulrich Keller)
    if (legend==TRUE) {
      legend("topright", legend=c("Adjusted Ev (retained)", "Adjusted Ev (unretained)", "Unadjusted Ev", "Random Ev"), col=c(AdjEvCol, AdjEvCol, EvCol, RndEvCol), pch = c(19, 21, 20, 20), lty = c(AdjEvLty, AdjEvLty, EvLty, RndEvLty))
      }

# Save the graph it if they have requested it be saved using the graphic 
# device they have specified.
    if (file != "" & typeof(file) == "character") {
      dev.copy(device=grdevice, height=height, width=width, file=file)
      dev.off()
      }
    }

  invisible(list(Retained = retained, AdjEv = AdjEv, Ev = Ev, RndEv = RndEv, Bias = Bias, SimEvs = SimEvs))
}

