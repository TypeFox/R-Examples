###############################################################################
# DiagnosticPlots:                                                            #
#                                                                             #
# Define functions to produce diagnostic plots of various outputs generated   #
# from the functions in this package, including AlphaNP, AlphaMLE, JMLE, and  #
# Qrefine.                                                                    #
#                                                                             #
###############################################################################

#install.packages("R.oo")
#library(R.oo)

################################################################################

# Diagnostic Plots for Nonparametric Alpha estimation
# inputs: (1) x: the output from AlphaNP; 
#         (2) nperson: choice of examinee to be investigated
# outputs: bar plot of sorted loss function values for each candidate alpha

setMethodS3("plot", "AlphaNP", function(x, nperson, cex.main=1, cex.legend=0.8, ...)  
  {
    loss <- x$loss.matrix[, nperson]
    my.order <- order(loss)
    title <- paste(paste(paste(paste(paste("Loss Function of Examinee"), nperson), "with Method \'"), x$method, sep=""), "\'", sep="")    
    
    pattern <-AlphaPermute(dim(x$alpha.est)[2])
    npattern <- dim(pattern)[1]
    density.vec <- rep(-1, npattern)
    pattern.name <- NULL
    for (i in 1:npattern) 
    {
      pattern.name <- c(pattern.name, paste(pattern[i,], collapse=""))
      if (all(pattern[i,] == x$alpha.est[nperson,])) density.vec[i] <- 10
    }
    
    par()
    on.exit()
    barplot(loss[my.order], main=title, xlab="Alpha", ylab="Loss", names.arg=pattern.name[my.order], 
            density=density.vec[my.order], cex.main=cex.main, las=3)
    legend(x="topleft", legend=c("Estimated Alpha"), density=c(10), cex=cex.legend)
  }
)


################################################################################

# Diagnostic Plots for MLE Alpha estimation
# inputs: (1) x: the output from AlphaLE; 
#         (2) nperson: choice of examinee to be investigated
# outputs: bar plot of sorted negative log-likelihood values for each candidate alpha

setMethodS3("plot", "AlphaMLE", function(x, nperson, cex.main=1, ...) 
{
  loglike <- -x$loglike.matrix[, nperson]
  my.order <- order(loglike)
  title <- paste("Negative Log-likelihood Function of Examinee", nperson)
  
  pattern <-AlphaPermute(dim(x$alpha.est)[2])
  npattern <- dim(pattern)[1]
  density.vec <- rep(-1, npattern)
  pattern.name <- NULL
  for (i in 1:npattern) 
  {
    pattern.name <- c(pattern.name, paste(pattern[i,], collapse=""))
    if (all(pattern[i,] == x$alpha.est[nperson,])) density.vec[i] <- 10
  }
  
  par()
  on.exit()
  barplot(loglike[my.order], main=title, xlab="Alpha", ylab="-LogLike", names.arg=pattern.name[my.order], 
          density=density.vec[my.order], las=3, cex.main=cex.main)
}
)


################################################################################

# Diagnostic Plots for JMLE estimation
# inputs: (1) x: the output from JMLE; 
#         (2) nperson: choice of examinee to be investigated
# outputs: (1) bar plot of unsorted loss function values for each candidate alpha
#          (2) bar plot of unsorted negative log-likelihood values for each candidate alpha

setMethodS3("plot", "JMLE", function(x, nperson, cex.main=0.9, ...) 
{
  pattern <-AlphaPermute(dim(x$alpha.est)[2])
  npattern <- dim(pattern)[1]
  density.vec.NP <- density.vec.JMLE <- rep(-1, npattern)
  pattern.name <- NULL
  for (i in 1:npattern) 
  {
    pattern.name <- c(pattern.name, paste(pattern[i,], collapse=""))
    if (all(pattern[i,] == x$alpha.est.NP[nperson,])) density.vec.NP[i] <- 10
    if (all(pattern[i,] == x$alpha.est[nperson,])) density.vec.JMLE[i] <- 10
  }
  
  par()
  par(mfrow=c(2, 1))
  on.exit()
  
  # Loss function in nonparametric estimation
  
  title <- paste(paste("Loss Function of Examinee", nperson), "with Method \'Hamming\'")    
  barplot(x$NP.loss.matrix[, nperson], main=title, xlab="Alpha", ylab="Loss", 
          names.arg=pattern.name, density=density.vec.NP, las=3, cex.main=cex.main)
  
  # Negative log-likelihood in JMLE estimation
  
  title <- paste("Negative Log-likelihood Function of Examinee", nperson)    
  barplot(-x$loglike.matrix[, nperson], main=title, xlab="Alpha", ylab="-LogLike", 
          names.arg=pattern.name, density=density.vec.JMLE, las=3, cex.main=cex.main) 
}
)


################################################################################

# Diagnostic Plots for Q refinement function
# inputs: (1) x: the output from Qrefine; 
# outputs: panel plots of each refined item showing the RSS of all candidate q-vectors

setMethodS3("plot", "Qrefine", function(x, filename="Qrefine.plot.png", cex.main=1, cex.lab=1, cex.axis=1, cex.legend=1, ...) 
{
  pattern <-AlphaPermute(dim(x$modified.Q)[2])
  npattern <- dim(pattern)[1]
  
  nmodified <- length(unique(x$modified.entries[, 1]))
  #ncol <- floor(sqrt(nmodified))
  #nrow <- ceiling(nmodified / ncol)
  #my.fac <- factor(1:(nrow * ncol + 1))
  nrow <- nmodified + 1
  my.fac <- factor(1:nrow)
  par()
  
  par(cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
  #par(mfrow=c((nrow + 1), ncol))
  par(mfrow=c(nrow, 1))
  on.exit()
  
  for (fac in my.fac)
  {
    i <- as.numeric(fac)
    
    if (i <= nmodified)
    {
      modified.item <- unique(x$modified.entries[,1])[i]
      RSS <- matrix(NA, dim(pattern)[1], 1)
      
      for (j in 1:npattern)
      {
        my.Q <- x$modified.Q
        my.Q[modified.item,] <- pattern[j,]
        Alpha.est.NP <- AlphaNP(x$Y, my.Q, x$gate, method="Hamming")
        est.ideal <- Alpha.est.NP$est.ideal
        diff <- x$Y[,i] - est.ideal[,i]
        RSS[j] <- sum(diff ^ 2)      
      }
      
      my.order <- order(RSS)
      
      title <- paste("RSS of Q-Vectors for Item", modified.item)
      
      pattern.name <- NULL
      density.vec <- rep(-1, npattern)
      angle.vec <- rep(45, npattern)
      
      for (k in 1:npattern)
      {
        pattern.name <- c(pattern.name, paste(pattern[k,], collapse=""))
        if (all(pattern[k,] == x$initial.Q[modified.item,])) density.vec[k] <- 10
        if (all(pattern[k,] == x$modified.Q[modified.item,])) 
        {
          density.vec[k] <- 20
          angle.vec[k] <-135
        }
      }    
      
      barplot(RSS[my.order], main=title, xlab="Q-Vector", ylab="RSS", names.arg=pattern.name[my.order], 
              las=3, angle=angle.vec[my.order], density=density.vec[my.order])
    }else
    {
      plot(1, ann=F, bty='n',type='n',xaxt='n',yaxt='n')
      if (i == nrow)
      {
        legend("topleft", legend=c("Modified", "Initial", "Other"), density=c(20, 20, -1), angle=c(135, 45, 45), horiz=T, cex=cex.legend)
      }
    }
  }  
  #dev.off()
}
)
