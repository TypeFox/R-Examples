# The .thresML is a minimum threshold below which the loglikelihood is not considered
# (i.e., if the calculated loglikelihood is lower to .thresML, it is set to .thresML for numerical reasons).
.thresML = -6000


BayesianMCMC <- function (xcont, xhist=NA, infhist=NA, suphist=NA, nbans=NA, seuil=NA,
                          nbpas=1000, nbchaines=3, confint=c(0.05, 0.95), dist="GEV",
                          apriori=function(...){1}, parameters0=NA, varparameters0=NA) {

 # This is the main function! Other functions are called in this one and are written below.
 # To know the meaning of the inputs and the outputs you can type help(BayesianMCMC).

 reject <- round(nbpas/10)   # number of initial steps that will be excluded from the posterior distribution
 nbpas <- nbpas + reject
 returnperiods <- 10^(seq(0.1, 4, by=.1))
 nonexceedF <- 1 - 1/returnperiods

 # Here, if initial parameters and parameter variances are not provided, a guess is made using the function .chooseparameters0 (see below)
 if (any(is.na(parameters0))) {
  parameters0_est <- parameters0
  parameters0 <- .chooseparameters0(xcont, xhist, infhist, suphist, dist)
  parameters0[!is.na(parameters0_est)] <- parameters0_est[!is.na(parameters0_est)]
 }
 if (any(is.na(varparameters0))) {
  varparameters0_est <- varparameters0
  varparameters0 <- (0.1*parameters0)^2
  varparameters0[!is.na(varparameters0_est)] <- varparameters0_est[!is.na(varparameters0_est)]
 }
 lpar <- length(parameters0)

 if(all(!is.na(c(nbans, seuil)))) {
   if (length(nbans) != length(seuil)) {
     stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): nbans and seuil should have the same length")
   }
   # Case only one threshold
   if (length(nbans)==1 & length(seuil)==1) {
     if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, seuil, nbans)))) {
       # Using censored information, historical flood known (Stedinger and Cohn, Naulet case b)
       long <- length(xhist)
     }
     if (all(is.na(c(xhist, suphist))) & all(!is.na(c(infhist, seuil, nbans)))) {
       # Using censored information, historical flood unknown (Stedinger and Cohn, Naulet case a)
       long <- length(infhist)
     }
     if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, seuil, nbans)))) {
       # Taking into account only flood estimation intervals
       long <- length(suphist)
     }
     seuil <- rep(seuil, long)
     nbans <- c((nbans), rep(0,(long-1)))
     xhist <- xhist
     infhist <- infhist
     suphist <- suphist
   }
 }

 # Here I initialize the arrays which will contain the results of the MCMC
 parameters <- array(data=NA, dim=c(nbpas, lpar, nbchaines))	# array 3D
 varparameters <- array(data=NA, dim=c(nbpas, lpar, nbchaines))    # array 3D
 vraisdist <- array(data=NA, dim=c(nbpas, nbchaines))
 #vraistest <- rep(NA, nbchaines)
 #nbsaut <- rep(NA, nbchaines)
 #propsaut <- rep(0, nbchaines)
 qq <- array(data=NA, dim=c(nbpas, length(nonexceedF), nbchaines))    # array 3D
 acceptance_binary <- array(data=NA,dim=c(nbpas, nbchaines)) # array 2D

 # Tests if some data are missing
 if(all(is.na(c(xcont, xhist, infhist, suphist, seuil)))) {
  stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): no input data")
 }
 if(all(is.na(seuil)) && (all(!is.na(xhist)) || all(!is.na(infhist)))) {
   stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): seuil data missing")
 }
 if(all(is.na(nbans)) && (all(!is.na(xhist)) || all(!is.na(infhist)))) {
   stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): nbans data missing")
 }

 # The following if-else selects the type of likelihood depending on the input provided
 if(all(is.na(c(xhist, infhist, suphist, nbans, seuil)))) {
  # Using only the systematic data
  funzionetest <- call(".lnvrais5", quote(parameters[1,,j]), quote(xcont), quote(dist))
  funzionecand <- call(".lnvrais5", quote(parameterscand), quote(xcont), quote(dist))
 }
 else if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, seuil, nbans)))) {
  # Using censored information, historical flood known (Stedinger and Cohn, Naulet case b)
  funzionetest <- call(".lnvrais1", quote(parameters[1,,j]), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais1", quote(parameterscand), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist))
 }
 else if (all(is.na(c(xhist, suphist))) & all(!is.na(c(infhist, seuil, nbans)))) {
  # Using censored information, historical flood unknown (Stedinger and Cohn, Naulet case a)
  funzionetest <- call(".lnvrais2", quote(parameters[1,,j]), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais2", quote(parameterscand), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
 }
 else if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, seuil, nbans)))) {
  # Taking into account only flood estimation intervals
  funzionetest <- call(".lnvrais4", quote(parameters[1,,j]), quote(xcont), quote(infhist), quote(suphist),
                                    quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais4", quote(parameterscand), quote(xcont), quote(infhist), quote(suphist),
                                    quote(nbans), quote(seuil), quote(dist))
 }
 else stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): inconsistency in input data")


 # Algorithm with repetitions
 # initialisation
 propsaut <- array(data=NA, dim=c(nbpas, nbchaines))
 propsaut_window <- array(data=NA, dim=c(nbpas, nbchaines))

 for (j in 1:nbchaines) {
  parameters[1,,j] <- .parameterscandMOD(parameters0, varparameters0, dist)   # first step (see .parameterscandMOD below)
  varparameters[1,,j] <- varparameters0
  vraistest <- eval(funzionetest) + log(apriori(parameters[1,,j]))   # it is a log-likelihhod
  vraisdist[1,j] <- vraistest
  qq[1,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[1,,j], dist=dist)   # first quantiles (see .quantilesMOD below)
  nbsaut <- 0
  acceptance_binary[1,j] <- 0
  for (i in 2:nbpas) {
   parameterscand <- .parameterscandMOD(parameters[i-1,,j], varparameters[i-1,,j], dist)   # other steps
   vraiscand <- eval(funzionecand) + log(apriori(parameterscand))   # it is a log-likelihhod
   valtest <- min((exp(vraiscand - vraistest)), 1)
   test <- runif(1)
   #if ((valtest > test) & (vraiscand > .thresML)) {
   if (valtest > test) {   # I move to the new set of parameters
    nbsaut <- nbsaut + 1
    acceptance_binary[i,j] <- 1
    parameters[i,,j] <- parameterscand
    vraistest <- vraiscand
   }
   else {   # I remain where I am
    parameters[i,,j] <- parameters[i-1,,j]
    acceptance_binary[i,j] <- 0
   }
   vraisdist[i,j] <- vraistest
   qq[i,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[i,,j], dist=dist)   # other quantiles

   # acceptance rate
   propsaut[i,j] <- nbsaut/i
   acceptance_window <- 1000
   first_window <- 1.5*acceptance_window + 1
   if (i>=first_window) {
    propsaut_window[i,j] <- mean(acceptance_binary[((i-acceptance_window+1):i),j])
   }
   else {
    propsaut_window[i,j] <- nbsaut/i
   }
   if (propsaut_window[i,j] < 0.33) {
    varparameters[i,,j] <- varparameters[i-1,,j]*(1+(propsaut_window[i,j]-0.34)/0.34/1000)
   }
   else if (propsaut_window[i,j] > 0.35) {
    varparameters[i,,j] <- varparameters[i-1,,j]*(1+(propsaut_window[i,j]-0.34)/0.34/1000)
   }
   else {
    varparameters[i,,j] <- varparameters[i-1,,j]
   }
  }
 }
 #removing the first iterations of the chains
 nbpas <- nbpas - reject
 qq <- qq[-c(1:reject),,]
 parameters <- parameters[-c(1:reject),,]
  dimnames(parameters) <- list(c(1:nbpas), paste("par", seq(1,lpar)), paste("chain", seq(1,nbchaines)))
 varparameters <- varparameters[-c(1:reject),,]
  dimnames(varparameters) <- list(c(1:nbpas), paste("varp", seq(1,lpar)), paste("chain", seq(1,nbchaines)))
 vraisdist <- vraisdist[-c(1:reject),]
  dimnames(vraisdist) <- list(c(1:nbpas), paste("chain", seq(1,nbchaines)))
 propsaut <- propsaut[-c(1:reject),]
  dimnames(propsaut) <- list(c(1:nbpas), paste("chain", seq(1,nbchaines)))
 propsaut_window <- propsaut_window[-c(1:reject),]
  dimnames(propsaut_window) <- list(c(1:nbpas), paste("chain", seq(1,nbchaines)))

 # The maximum likelihood is here performed simply taking the maximum of vraisdist
 logML <- max(vraisdist)
 dummy1 <- which.max(apply(vraisdist, 2, max))
 dummy2 <- apply(vraisdist, 2, which.max)[dummy1]
 parametersML <- parameters[dummy2,,dummy1]
 intervals <- apply(qq, 2, quantile, probs=confint, na.rm=TRUE)
 qqML <- qq[dummy2,,dummy1]
 output <- list(xcont=xcont, xhist=xhist, infhist=infhist, suphist=suphist, nbans=nbans, seuil=seuil,
                nbpas=nbpas, nbchaines=nbchaines, dist=dist, confint=confint, apriori=apriori,
                parameters=parameters, quantiles=qq, varparameters=varparameters,
                parameters0=parameters0, varparameters0=varparameters0,
                vraisdist=vraisdist, propsaut=propsaut,
                returnperiods=returnperiods, intervals=intervals,
                parametersML=parametersML, quantilesML=qqML, logML=logML)
 class(output) <- "BayesianMCMC"
 return(output)
}




BayesianMCMCcont <- function (x, nbpas=NA) {
 if(is.na(nbpas)) nbpas <- x$nbpas
 parameters0 <- x$parametersML
 varparameters0 <- rowMeans(apply(x$parameters, c(2,3), var))

 output <- BayesianMCMC (xcont=x$xcont, xhist=x$xhist, infhist=x$infhist, suphist=x$suphist,
               nbans=x$nbans, seuil=x$seuil, nbpas=nbpas, nbchaines=x$nbchaines,
               confint=x$confint, dist=x$dist,
               apriori=x$apriori,
               parameters0=parameters0, varparameters0=varparameters0)
 if(output$logML < x$logML) {
  output$logML <- x$logML
  output$parametersML <- x$parametersML
  output$quantilesML <- x$quantilesML
 }
 return(output)
}




# ----------------------------- #

print.BayesianMCMC <- function (x, ...) {
 dummy <- data.frame(cbind(x$quantilesML, t(x$intervals)), row.names=signif(x$returnperiods, 4))
 names(dummy)[1] <- "ML"
 print(dummy)
}


# ----------------------------- #

plot.BayesianMCMC <- function (x, which=1, ask=FALSE, ...) {
 if (ask) {
  op <- par(ask = TRUE)
  on.exit(par(op))
 }
 show <- rep(FALSE, 100)
 show[which] <- TRUE
 if (show[1]) .plotdiagnMCMC02(x, ...)
 if (show[2]) .plotdiagnMCMC01(x, ...)
 if (show[3]) .plotdiagnMCMC04(x, ...)
 if (show[4]) .plotdiagnMCMC05(x, ...)
 if (show[5]) .plotdiagnMCMC06(x, ...)   # a-priori distribution
}


# --------------------- #

.plotdiagnMCMC01 <- function(x, ...) {
 # Diagnostic plot of the parameters
 lpar <- dim(x$parameters)[2]
 #graphics.off()
 #x11()
 op <- par(mfrow=c(lpar, 3))
  for (j in 1:lpar) {
   limiti <- range(x$parameters[,j,])
   plot(x$parameters[,j,1], type="l", col=2, ylim=limiti, ylab=paste("par",j), xlab="")
   for (i in 2:x$nbchaines) {
    lines(x$parameters[,j,i], col=1+i)
   }
   abline(h=x$parametersML[j])
   #hist(x$parameters[,j,1], border=2, breaks=11, #seq(limiti[1], limiti[2], length=11),
   #     xlab=paste("par",j), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   #for (i in 2:x$nbchaines) {
   # hist(x$parameters[,j,i], border=1+i, breaks=11, add=TRUE)
   #}
   #abline(v=x$parametersML[j])
   ht <- hist(x$parameters[,j,], plot=FALSE)
   br <- ht$breaks
   hist(x$parameters[,j,1], border=2, breaks=br, xlab=paste("par",j), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   for (i in 2:x$nbchaines) {
    hist(x$parameters[,j,i], border=1+i, breaks=br, add=TRUE)
   }
   abline(v=x$parametersML[j])
   plot(x$varparameters[,j,1], type="l", col=2, ylim=range(x$varparameters[,j,]), ylab=paste("var par",j), xlab="")
   for (i in 2:x$nbchaines) {
    lines(x$varparameters[,j,i], col=1+i)
   }
  }
 par(op)
}

.plotdiagnMCMC02 <- function(x, ...) {
 # Plot of the frequency curve
 T <- c(1,1000)
 X <- c(0, 1.3*max(c(x$xcont, x$xhist, x$infhist, x$suphist), na.rm=TRUE))
 plot(T, X, type="n", log="x", ...)
 grid(equilogs=FALSE)
 ret <- .pointspos3 (x$xcont, x$xhist, x$infhist, x$suphist, x$nbans, x$seuil)
 lines(x$returnperiods, x$quantilesML)
 lines(x$returnperiods, x$intervals[1,], lty=2)
 lines(x$returnperiods, x$intervals[2,], lty=2)
 invisible(ret)
}

.plotdiagnMCMC03 <- function(x, ...) {
 Nsim=10000
 #if(all(is.na(c(x$xhist, x$infhist, x$suphist, x$seuil)))) {
 # # Calcul sur les seules données systèmatiques
 # funzione <- call(".lnvrais5", quote(x$parametersML), quote(xcont), quote(dist))
 #}
 #else if (all(is.na(c(x$infhist, x$suphist))) & all(!is.na(c(x$xhist, x$seuil, x$nbans)))) {
 # # Calcul avec info censurée mais débits historiques connus (Stedinger et Cohn, Naulet cas b)
 # funzione <- call(".lnvrais1", quote(x$parametersML), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist))
 #}
 #else if (all(is.na(c(x$xhist, x$suphist))) & all(!is.na(c(x$infhist, x$seuil, x$nbans)))) {
 # # Calcul avec info censurée mais débits historiques non connus (Stedinger et Cohn, Naulet cas a)
 # funzione <- call(".lnvrais2", quote(x$parametersML), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
 #}
 #else if (all(is.na(c(x$xhist))) & all(!is.na(c(x$infhist, x$suphist, x$seuil, x$nbans)))) {
 # # Calcul avec prise en compte des seuls intervalles d'estimation de débit
 # funzione <- call(".lnvrais4", quote(x$parametersML), quote(xcont), quote(infhist), quote(suphist),
 #                                   quote(nbans), quote(seuil), quote(dist))
 #}
}

.plotdiagnMCMC04 <- function(x, ...) {
 # Plot of the acceptance rate and the likelihood
 op <- par(mfrow=c(2, 2))
  limiti <- range(x$vraisdist)
  plot(x$vraisdist[,1], type="l", col=2, ylim=limiti, ylab=paste("ln likelihood"), xlab="")
   for (i in 2:x$nbchaines) {
    lines(x$vraisdist[,i], col=1+i)
   }
  ht <- hist(x$vraisdist, plot=FALSE)
  br <- ht$breaks
  hist(x$vraisdist[,1], border=2, breaks=br, xlab=paste("ln likelihood"), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   for (i in 2:x$nbchaines) {
    hist(x$vraisdist[,i], border=1+i, breaks=br, add=TRUE)
   }

  limiti <- range(x$propsaut)
  plot(x$propsaut[,1], type="l", col=2, ylim=limiti, ylab=paste("acceptance rate"), xlab="")
   for (i in 2:x$nbchaines) {
    lines(x$propsaut[,i], col=1+i)
   }
  ht <- hist(x$propsaut, plot=FALSE)
  br <- ht$breaks
  hist(x$propsaut[,1], border=2, breaks=br, xlab=paste("acceptance rate"), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   for (i in 2:x$nbchaines) {
    hist(x$propsaut[,i], border=1+i, breaks=br, add=TRUE)
   }
 par(op)
}



.plotdiagnMCMC05 <- function(x, ...) {
 # A posteriori distribution of the parameters
 if (length(x$parametersML) == 2) {
  plot(x$parameters[,1:2,1], col=2, pch=".")
  for (i in 2:x$nbchaines) {
   points(x$parameters[,1:2,i], col=1+i, pch=".")
  }
 }
 else if (length(x$parametersML) == 3) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  #layout(matrix(c(1,1,1,2,2,2,3,3,3,4,5,6), 6, 2, byrow=FALSE))
  layout(matrix(c(1,2,3,0), 2, 2, byrow=FALSE))
  plot(x$parameters[,c(1,2),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(1,2),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[1], x$parametersML[2], pch=19, cex=1.5, col=7)

  plot(x$parameters[,c(1,3),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(1,3),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[1], x$parametersML[3], pch=19, cex=1.5, col=7)

  plot(x$parameters[,c(3,2),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(3,2),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[3], x$parametersML[2], pch=19, cex=1.5, col=7)

  #  par(mar=c(1,1,1,1))
  #  plot(density(x$parameters[,1,1]), col=1, main="", axes=FALSE)
  #  for (i in 2:x$nbchaines) {
  #   lines(density(x$parameters[,1,i]), col=i)
  #  }
  #  mtext("par 1", 3, -1.5, adj=0.03, cex=.8)
  #  box()

  #  plot(density(x$parameters[,2,1]), col=1, main="", axes=FALSE)
  #  for (i in 2:x$nbchaines) {
  #   lines(density(x$parameters[,2,i]), col=i)
  #  }
  #  mtext("par 2", 3, -1.5, adj=0.03, cex=.8)
  #  box()

  #  plot(density(x$parameters[,3,1]), col=1, main="", axes=FALSE)
  #  for (i in 2:x$nbchaines) {
  #   lines(density(x$parameters[,3,i]), col=i)
  #  }
  #  mtext("par 3", 3, -1.5, adj=0.03, cex=.8)
  #  box()
  par(def.par)
 }
}



.plotdiagnMCMC06 <- function(x, ...) {
 # plot the a-priori distribution of the parameters
 if (length(x$parametersML) == 2) {
  plot(x$parameters[,1:2,1], type="n", ...)
  xr <- range(x$parameters[,1,1])
  yr <- range(x$parameters[,2,1])
  xx <- seq(xr[1], xr[2], length=50)
  yy <- seq(yr[1], yr[2], length=50)
  xy <- expand.grid(xx, yy)
  zz <- rep(NA, length=dim(xy)[1])
  test <- rep(1, length=dim(xy)[1])
  for (i in 1:dim(xy)[1]) {
   zz[i] <- x$apriori(as.numeric(xy[i,]))
  }
  if((identical(test,zz)) == TRUE){
     warning("Impossible to plot the a-priori distribution of the parameters : no a-priori information")
  }
  else contour(xx, yy, zz, add=TRUE, col="darkgray")
 }
 else if (length(x$parametersML) == 3) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  x1r <- range(x$parameters[,1,1])
  x2r <- range(x$parameters[,2,1])
  x3r <- range(x$parameters[,3,1])
  xx1 <- seq(x1r[1], x1r[2], length=20)
  xx2 <- seq(x2r[1], x2r[2], length=20)
  xx3 <- seq(x3r[1], x3r[2], length=20)
  xx123 <- expand.grid(xx1, xx2, xx3)
  N <- dim(xx123)[1]
  zz <- rep(NA, length=N)
  test <- rep(1, length=N)
  for (i in 1:N) {
   zz[i] <- x$apriori(as.numeric(xx123[i,]))
  }
  if((identical(test,zz)) == TRUE){
     warning("Impossible to plot the a-priori distribution of the parameters : no a-priori information")
  }
  else{
    layout(matrix(c(1,2,3,0), 2, 2, byrow=FALSE))
    plot(x$parameters[,c(1,2),1], type="n")
    contour(xx1, xx2, tapply(zz, xx123[,c(1,2)], sum), add=TRUE, col="darkgray")

    plot(x$parameters[,c(1,3),1], type="n")
    contour(xx1, xx3, tapply(zz, xx123[,c(1,3)], sum), add=TRUE, col="darkgray")

    plot(x$parameters[,c(3,2),1], type="n")
    contour(xx3, xx2, tapply(zz, xx123[,c(3,2)], sum), add=TRUE, col="darkgray")
  }
  par(def.par)
 }
}





# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

.pointspos3 <- function (xcont, xhist, infhist, suphist, nbans, seuil, ...)
{
 if(all(is.na(c(xhist, infhist, suphist, seuil)))) {
  # Using only the systematic data (Cunnane plotting position)
  xcont <- sort(xcont, decreasing=TRUE)
  F <- c(1:length(xcont))
  F <- 1 - ((F - 0.4)/(length(xcont) + 1 - 2*0.4))
  T <- 1/(1 - F)
  x <- xcont
  #plot(T, x, log="x", type="n")
  #grid(equilogs = FALSE)
  points(T, x)
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, seuil, nbans)))) {
  # Using censored information, historical flood known (Stedinger and Cohn, Naulet case b)
  # Calculates the probabilities associated with the thresholds
  unseuil<-unique(seuil)
  unseuil <- sort(unseuil)
  unnbans<-rep(NA,length(unseuil))
  uncruessup<-rep(NA,length(unseuil))
  for (i in 1:length(unseuil)) {unnbans[i]<-sum(nbans[seuil==unseuil[i]])
                                uncruessup[i]<-length(xhist[seuil==unseuil[i] & xhist>=unseuil[i]])
                                }
  nbiter<-100
  toto<-.points1(uncruessup, unnbans, unseuil, nbiter)
  if (length(dim(toto))>0) {
    toto <- apply(toto,1,median)
  }
  #Tseuil <- 1/(1-toto)

  unseuil<-c(0,unseuil)
  toto<-c(0,toto)

  xhist[xhist==-1]<-NA
  xhist<-na.omit(xhist)
 
  
  
  # Calculates the probabilities of (recent and historical) floods
  xcrues<-c(xcont,xhist)
  allcrues <- matrix(data=NA,nrow=length(xcrues),ncol=2)
  isCont <- c(rep(TRUE,length(xcont)), rep(FALSE,length(xhist)))
  allcrues[,1] <- xcrues
  allcrues[,2] <- isCont
  allcrues <- data.frame(allcrues)
  colnames(allcrues)<-c("xcrues","isCont")
  allcrues <- allcrues[order(allcrues$xcrues),]
  xcrues<-allcrues$xcrues

  for (i in 1:length(unseuil)) {
    if (i<length(unseuil)) {
      xcruesseuil<-xcrues[xcrues>=unseuil[i] & xcrues<unseuil[i+1]]
      xxseuil <- 1 + length(xcruesseuil) - rank(xcruesseuil, ties.method="first")
      xxseuil <- toto[i+1]+(toto[i]-toto[i+1])*((xxseuil - 0.4)/(length(xcruesseuil) + 1 - 2*0.4))
      }
    else {
      xcruesseuil<-xcrues[xcrues>=unseuil[i]]
      xxseuil <- 1 + length(xcruesseuil) - rank(xcruesseuil, ties.method="first")
      xxseuil <- 1+(toto[i]-1)*((xxseuil - 0.4)/(length(xcruesseuil) + 1 - 2*0.4))
      }
    if (i==1) {xx<-xxseuil} else {xx<-c(xx,xxseuil)}
  }
  T <- 1/(1 - xx)
  
  

  # Plots the graph
  #points(T, xcrues, pch=3)
  xcont<-xcrues[allcrues$isCont == 1]
  Tcont<-rep(NA,length(xcont))
  Tcont<-T[allcrues$isCont == 1]
  #plot(T, xcrues, log="x", type="n") # creates an empty graph
  #grid(equilogs = FALSE)
  points(Tcont, xcont, pch=3)

  xhist<-xcrues[allcrues$isCont == 0]
  Txhist<-rep(NA,length(xhist))
  Txhist<-T[allcrues$isCont == 0]
  points(Txhist, xhist, pch=19)
  invisible(cbind(T,xcrues))
 }
 else if (all(is.na(c(xhist, suphist))) & all(!is.na(c(infhist, seuil, nbans)))) {
  # Using censored information, historical flood unknown (Stedinger and Cohn, Naulet case a)
  # Calculates the probabilities associated with the thresholds
  unseuil<-unique(seuil)
  unseuil <- sort(unseuil)
  unnbans<-rep(NA,length(unseuil))
  uncruessup<-rep(NA,length(unseuil))
  for (i in 1:length(unseuil)) {unnbans[i]<-sum(nbans[seuil==unseuil[i]])
                                uncruessup[i]<-length(infhist[seuil==unseuil[i] & infhist>=unseuil[i]])
                                }
  nbiter<-100
  toto<-.points1(uncruessup, unnbans, unseuil, nbiter)
  if (length(dim(toto))>0) {
    toto <- apply(toto,1,median)
  }
  #Tseuil <- 1/(1-toto)

  unseuil<-c(0,unseuil)
  toto<-c(0,toto)
  
  infhist[infhist==-1]<-NA
  infhist<-na.omit(infhist)

  # Calculates the probabilities of (recent and historical) floods
  xcrues <- c(xcont,infhist)
  allcrues <- matrix(data=NA,nrow=length(xcrues),ncol=2)
  isCont <- c(rep(TRUE,length(xcont)), rep(FALSE,length(infhist)))
  allcrues[,1] <- xcrues
  allcrues[,2] <- isCont
  allcrues <- data.frame(allcrues)
  colnames(allcrues)<-c("xcrues","isCont")
  allcrues <- allcrues[order(allcrues$xcrues),]
  xcrues<-allcrues$xcrues

  for (i in 1:length(unseuil)) {
    if (i<length(unseuil)) {
      xcruesseuil<-xcrues[xcrues>=unseuil[i] & xcrues<unseuil[i+1]]
      xxseuil <- 1 + length(xcruesseuil) - rank(xcruesseuil, ties.method="first")
      xxseuil <- toto[i+1]+(toto[i]-toto[i+1])*((xxseuil - 0.4)/(length(xcruesseuil) + 1 - 2*0.4))
      }
    else {
      xcruesseuil<-xcrues[xcrues>=unseuil[i]]
      xxseuil <- 1 + length(xcruesseuil) - rank(xcruesseuil, ties.method="first")
      xxseuil <- 1+(toto[i]-1)*((xxseuil - 0.4)/(length(xcruesseuil) + 1 - 2*0.4))
      }
    if (i==1) {xx<-xxseuil} else {xx<-c(xx,xxseuil)}
  }
  T <- 1/(1 - xx)

  # Plots the graph
  #points(T, xcrues, pch=3)
  xcont<-xcrues[allcrues$isCont == 1]
  Tcont<-rep(NA,length(xcont))
  Tcont<-T[allcrues$isCont == 1]
  #plot(T, xcrues, log="x", type="n") # creates an empty graph
  #grid(equilogs = FALSE)
  points(Tcont, xcont, pch=3)

  infhist<-xcrues[allcrues$isCont == 0]
  Tinfhist<-rep(NA,length(infhist))
  Tinfhist<-T[allcrues$isCont == 0]
  points(Tinfhist, infhist, pch=19)
  if(length(infhist)!=0){
    segments(Tinfhist, infhist, Tinfhist, 10*max(xcrues), lty=3)
  }
  
  invisible(cbind(T,xcrues))
 }
 else if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, seuil, nbans)))) {
  # Taking into account only flood estimation intervals
  # Calculates the probabilities associated with the thresholds
  unseuil<-unique(seuil)
  unseuil <- sort(unseuil)
  unnbans<-rep(NA,length(unseuil))
  uncruessup<-rep(NA,length(unseuil))
  meanhist<-(infhist + suphist)/2
  for (i in 1:length(unseuil)) {unnbans[i]<-sum(nbans[seuil==unseuil[i]])
                                uncruessup[i]<-length(meanhist[seuil==unseuil[i] & meanhist>=unseuil[i]])
                                }
  nbiter<-100
  toto<-.points1(uncruessup, unnbans, unseuil, nbiter)
  if (length(dim(toto))>0) {
    toto <- apply(toto,1,median)
  }
  #Tseuil <- 1/(1-toto)

  unseuil<-c(0,unseuil)
  toto<-c(0,toto)
  
  infhist[infhist==-1]<-NA
  infhist<-na.omit(infhist)
  suphist[suphist==-1]<-NA
  suphist<-na.omit(suphist)
  meanhist[meanhist==-1]<-NA
  meanhist<-na.omit(meanhist)


  # Calculates the probabilities of (recent and historical) floods
  xcrues<-c(xcont,meanhist)
  allcrues <- matrix(data=NA,nrow=length(xcrues),ncol=4)
  isCont <- c(rep(TRUE,length(xcont)), rep(FALSE,length(meanhist)))
  allcrues[,1] <- xcrues
  allcrues[,2] <- isCont
  allcrues[,3] <- c(rep(NA,length(xcont)), infhist)
  allcrues[,4] <- c(rep(NA,length(xcont)), suphist)
  allcrues <- data.frame(allcrues)
  colnames(allcrues)<-c("xcrues","isCont","infhist","suphist")
  allcrues <- allcrues[order(allcrues$xcrues),]
  xcrues<-allcrues$xcrues

  for (i in 1:length(unseuil)) {
    if (i<length(unseuil)) {
      xcruesseuil<-xcrues[xcrues>=unseuil[i] & xcrues<unseuil[i+1]]
      xxseuil <- 1 + length(xcruesseuil) - rank(xcruesseuil, ties.method="first")
      xxseuil <- toto[i+1]+(toto[i]-toto[i+1])*((xxseuil - 0.4)/(length(xcruesseuil) + 1 - 2*0.4))
      }
    else {
      xcruesseuil<-xcrues[xcrues>=unseuil[i]]
      xxseuil <- 1 + length(xcruesseuil) - rank(xcruesseuil, ties.method="first")
      xxseuil <- 1+(toto[i]-1)*((xxseuil - 0.4)/(length(xcruesseuil) + 1 - 2*0.4))
      }
    if (i==1) {xx<-xxseuil} else {xx<-c(xx,xxseuil)}
  }
  T <- 1/(1 - xx)

  # Plots the graph
  #  points(T, xcrues, pch=3)
  xcont<-xcrues[allcrues$isCont == 1]
  Tcont<-rep(NA,length(xcont))
  Tcont<-T[allcrues$isCont == 1]
  #plot(T, xcrues, log="x", type="n") # creates an empty graph
  #grid(equilogs = FALSE)
  points(Tcont, xcont, pch=3)

  meanhist<-xcrues[allcrues$isCont == 0]
  Tmeanhist<-rep(NA,length(meanhist))
  Tmeanhist<-T[allcrues$isCont == 0]
  points(Tmeanhist, meanhist, pch=19)
  infhist<-allcrues$infhist[allcrues$isCont == 0]
  suphist<-allcrues$suphist[allcrues$isCont == 0]
  segments(Tmeanhist, infhist, Tmeanhist, suphist, lty=3)
  invisible(cbind(T,xcrues))
 }
 else stop(".pointspos3(xcont, xhist, infhist, suphist, nbans, seuil): inconsistency in input data")
}

#-----------------------------------------------------------------------------#
.points1 <- function (cruessup, nbans, seuil, nbiter){
 #-------- Calculates the empirical probabilities of thresholds ----------------------#
 #---for a number of years and a number of exceedances associated to each threshold ---#
 toto<-runif(sum(nbans)*nbiter,min=0,max=1)
 toto <- matrix(toto,ncol=sum(nbans))      # you generate a matrix of numbers between 0 and 1 drawn from a uniform distribution, standing for probabilities
 annee<-1
 posseuil<-matrix(data=NA,nrow=nbiter,ncol=length(nbans))
 for (i in 1:length(nbans)) {
   totoseuil<-toto[,annee:(annee+nbans[i]-1)]  # for each threshold you take from the matrix only as many columns as you have years associated to this threshold
   if(is.vector(totoseuil)){
     totoseuil<-sort(totoseuil)                 # you sort your probabilities
     if(cruessup[i]!=0){            # cruessup is as defined in .pointspos3 (the function using .points1)
      totoseuil<-totoseuil[nbans[i]+1-cruessup[i]]
     }
     else {
      totoseuil<-totoseuil[nbans[i]-cruessup[i]]
     }
   }  
   else{
     totoseuil<-apply(totoseuil,1,sort)  
     if(cruessup[i]!=0){
      totoseuil<-totoseuil[nbans[i]+1-cruessup[i],]
     }
     else {
      totoseuil<-totoseuil[nbans[i]-cruessup[i],]
     }
   }   
    
   posseuil[,i]<-totoseuil
   annee<-annee+nbans[i]
 }
 posseuil<-apply(posseuil,1,sort)
 return(posseuil)
}

# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.chooseparameters0 <- function (xcont, xhist, infhist, suphist, dist) {   # This function is used to select the initial parameters for the different distributions
  campionissimo <- sort(c(xcont, xhist, infhist, suphist))   # Unexpected value, but not completely false either
  ll <- as.numeric(Lmoments(campionissimo))
  if (dist=="GEV") {
   parameters0 <- unlist(par.GEV(ll[1], ll[2], ll[4]))
  }
  else if (dist=="NORM") {
   parameters0 <- c(mean(campionissimo), sd(campionissimo))
  }
  else if (dist=="EXP") {
   parameters0 <- unlist(par.exp(ll[1], ll[2]))
   if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont)
  }
  else if (dist=="GENLOGIS") {
   parameters0 <- unlist(par.genlogis(ll[1], ll[2], ll[4]))
  }
  else if (dist=="GENPAR") {
   parameters0 <- unlist(par.genpar(ll[1], ll[2], ll[4]))
   if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont)
  }
  else if ((dist=="GUMBEL")||(dist=="EV1")) {
   parameters0 <- unlist(par.gumb(ll[1], ll[2]))
  }
  else if (dist=="KAPPA") {
   parameters0 <- unlist(par.kappa(ll[1], ll[2], ll[4], ll[5]))
  }
  else if ((dist=="LOGNORM")||(dist=="LN3")) {
   parameters0 <- unlist(par.lognorm(ll[1], ll[2], ll[4]))
  }
  else if ((dist=="LN")||(dist=="LN2")) {
   parameters0 <- c(mean(log(campionissimo)), sd(log(campionissimo)))
  }
  else if ((dist=="P3")||(dist=="GAM")) {
   parameters0 <- unlist(par.gamma(ll[1], ll[2], ll[4])[1:3])
   if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont) # Positive skewness!!!
  }
  else stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): distribution unknown")
  return(parameters0)
}


# ---------------------------- #

.parameterscandMOD <- function (parameters, varparameters, dist) {
  # Perform a step using different distributions (normal or lognormal) depending on the parameter
  # (essentially if it must be positive or can be also negative)
  if (dist=="GEV") {
   # I know that the scale parameter is positive
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else if (dist=="NORM") {
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter1 <- rnorm(1, mean=parameters, sd=sqrt(varparameters))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter1, parameter2)
  }
  else if (dist=="EXP") {
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter1 <- rnorm(1, mean=parameters, sd=sqrt(varparameters))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter1, parameter2)
  }
  else if (dist=="GENLOGIS") {
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else if (dist=="GENPAR") {
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else if ((dist=="GUMBEL")||(dist=="EV1")) {
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter1 <- rnorm(1, mean=parameters, sd=sqrt(varparameters))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter1, parameter2)
  }
  else if (dist=="KAPPA") {
   parameterscand <- rnorm(rep(1, 4), mean=parameters, sd=sqrt(varparameters))
  }
  else if ((dist=="LOGNORM")||(dist=="LN3")) {
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else if ((dist=="LN")||(dist=="LN2")) {
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter1 <- rnorm(1, mean=parameters, sd=sqrt(varparameters))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter1, parameter2)
  }
  else if ((dist=="P3")||(dist=="GAM")) {
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): distribution unknown")
  return(parameterscand)
}



# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.quantilesMOD <- function (F, parameters, dist="GEV") {
 # Calculates the quantiles for a given distribution, a given parameter set and a given non-exceedance probability
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  qq <- invF.GEV(F, xi, alfa, k)
 }
 else if (dist=="NORM") {
  qq <- qnorm(F, mean=parameters[1], sd=parameters[2])
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  qq <- invF.exp(F, xi, alfa)
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  qq <- invF.genlogis(F, xi, alfa, k)
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  qq <- invF.genpar(F, xi, alfa, k)
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  qq <- invF.gumb(F, xi, alfa)
 }
 else if (dist=="KAPPA") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  h <- parameters[4]
  qq <- invF.kappa(F, xi, alfa, k, h)
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  qq <- invF.lognorm(F, xi, alfa, k)
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  #varlogLN <- log(1 + parameters[2]/parameters[1]^2)
  #meanlogML <- log(parameters[1]) - varlogLN/2
  #qq <- qlnorm(F, meanlog=meanlogML, sdlog=sqrt(varlogLN))
  qq <- qlnorm(F, meanlog=parameters[1], sdlog=parameters[2])
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  qq <- invF.gamma(F, xi, beta, alfa)
 }
 else stop(".quantilesMOD(F, parameters, dist): distribution unknown")

 return(qq)
}

# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.lnvrais5 <- function (parameters, xcont, dist="GEV") {
 # Using only the systematic data
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
   if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
    lnvrais <- sum(log(f.GEV(xcont, xi, alfa, k)))
   }
  else lnvrais <- .thresML
 }
 else if (dist=="NORM") {
  lnvrais <- sum(log(dnorm(xcont, mean=parameters[1], sd=parameters[2])))
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvrais <- sum(log(f.exp(xcont, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(f.genlogis(xcont, xi, alfa, k)))
  }
  else lnvrais <- .thresML
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0 & all(xcont > xi)) {
   lnvrais <- sum(log(f.genpar(xcont, xi, alfa, k)))
  }
  else lnvrais <- .thresML
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvrais <- sum(log(f.gumb(xcont, xi, alfa)))
 }
 #else if (dist=="KAPPA") {
 # xi <- parameters[1]
 # alfa <- parameters[2]
 # k <- parameters[3]
 # h <- parameters[4]
 # if (sum((k*(xcont-xi)/alfa) > 1)==0) {
 #  lnvrais <- sum(log(f.kappa(xcont, xi, alfa, k, h)))
 # }
 # else lnvrais <- .thresML
 #}
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvrais <- .thresML
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvrais <- sum(log(dlnorm(xcont, meanlog=parameters[1], sdlog=parameters[2])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  if (alfa > 0 && is.na(beta)!=TRUE) {
   lnvrais <- sum(log(f.gamma(xcont, xi, beta, alfa)))
  }
  else lnvrais <- .thresML
 }
 else stop(".lnvrais5(parameters, xcont, dist): distribution unknown")
 if (is.nan(lnvrais)) lnvrais <- .thresML
 if (lnvrais < .thresML) lnvrais <- .thresML

 return(lnvrais)
}





# ---------------- #

.lnvrais1 <- function (parameters, xcont, xhist, nbans, seuil, dist="GEV") {
 nbans <- .datachange(nbans)
 # Case of xhist=-1
 nbans[xhist==-1] <- nbans[xhist==-1] + 1
 xhist[xhist==-1]<-NA
 xhist <- na.omit(xhist)
 
 # Using censored information, historical flood known (Stedinger and Cohn, Naulet case b)
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
   if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
    lnvraiscont <- sum(log(f.GEV(xcont, xi, alfa, k)))
   }
  else lnvraiscont <- .thresML
   if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
    lnvraishist <- sum((nbans-1) * log(F.GEV(seuil, xi, alfa, k)))
   }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(xhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(f.GEV(xhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(dnorm(xcont, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- sum((nbans - 1) * log(pnorm(seuil, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- lnvraishist + sum(log(dnorm(xhist, mean=parameters[1], sd=parameters[2])))
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.exp(xcont, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.exp(seuil, xi, alfa))) + sum(log(f.exp(xhist, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.genlogis(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.genlogis(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(xhist-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(f.genlogis(xhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0 && all(xcont > xi)) {
   lnvraiscont <- sum(log(f.genpar(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0 && all(seuil > xi)) {
   lnvraishist <- sum((nbans - 1) * log(F.genpar(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(xhist-xi)/alfa) > 1)==0 && all(xhist > xi)) {
   lnvraishist <- lnvraishist + sum(log(f.genpar(xhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.gumb(xcont, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.gumb(seuil, xi, alfa)) + sum(log(f.gumb(xhist, xi, alfa))))
 }
 #else if (dist=="KAPPA") {
 # xi <- parameters[1]
 # alfa <- parameters[2]
 # k <- parameters[3]
 # h <- parameters[4]
 # lnvraiscont <- sum(log(f.kappa(xcont, xi, alfa, k, h)))
 # lnvraishist <- sum((nbans - 1) * log(F.kappa(seuil, xi, alfa, k, h)) + sum(log(1 - F.kappa(xhist, xi, alfa, k, h))))
 #}
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.lognorm(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(xhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(f.lognorm(xhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(dlnorm(xcont, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- sum((nbans - 1) * log(plnorm(seuil, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- lnvraishist + sum(log(dlnorm(xhist, meanlog=parameters[1], sdlog=parameters[2])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  if (alfa > 0 && is.na(beta)!=TRUE) {
   lnvraiscont <- sum(log(f.gamma(xcont, xi, beta, alfa)))
   lnvraishist <- sum((nbans - 1) * log(F.gamma(seuil, xi, beta, alfa))) +
                  sum(log(f.gamma(xhist, xi, beta, alfa)))
  }
  else {
   lnvraiscont <- .thresML
   lnvraishist <- .thresML
  }
 }
 else stop(".lnvrais1(parameters, xcont, xhist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresML
 if (lnvrais < .thresML) lnvrais <- .thresML

 return(lnvrais)
}






# ---------------- #

.lnvrais2 <- function (parameters, xcont, infhist, nbans, seuil, dist="GEV") {
 nbans <- .datachange(nbans)
 # Case of infhist=-1
 nbans[infhist==-1] <- nbans[infhist==-1] + 1
 infhist[infhist==-1]<-NA
 infhist <- na.omit(infhist)
 # Using censored information, historical flood unknown (Stedinger and Cohn, Naulet case a)
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.GEV(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans-1) * log(F.GEV(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(1 - F.GEV(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(dnorm(xcont, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- sum((nbans - 1) * log(pnorm(seuil, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- lnvraishist + sum(log(1 - pnorm(infhist, mean=parameters[1], sd=parameters[2])))
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.exp(xcont, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.exp(seuil, xi, alfa))) + sum(log(1 - F.exp(infhist, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.genlogis(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.genlogis(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(1 - F.genlogis(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0 && all(xcont > xi)) {
   lnvraiscont <- sum(log(f.genpar(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0 && all(seuil > xi)) {
   lnvraishist <- sum((nbans - 1) * log(F.genpar(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(infhist-xi)/alfa) > 1)==0 && all(infhist > xi)) {
   lnvraishist <- lnvraishist + sum(log(1 - F.genpar(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.gumb(xcont, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.gumb(seuil, xi, alfa))) + sum(log(1 - F.gumb(infhist, xi, alfa)))
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.lognorm(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(1 - F.lognorm(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(dlnorm(xcont, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- sum((nbans - 1) * log(plnorm(seuil, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- lnvraishist + sum(log(1 - plnorm(infhist, meanlog=parameters[1], sdlog=parameters[2])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  if (alfa > 0 && is.na(beta)!=TRUE) {
    lnvraiscont <- sum(log(f.gamma(xcont, xi, beta, alfa)))
    lnvraishist <- sum((nbans - 1) * log(F.gamma(seuil, xi, beta, alfa))) +
                   sum(log(1 - F.gamma(infhist, xi, beta, alfa)))
  }
  else {
   lnvraiscont <- .thresML
   lnvraishist <- .thresML
  }
 }
 else stop(".lnvrais2(parameters, xcont, infhist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresML
 if (lnvrais < .thresML) lnvrais <- .thresML

 return(lnvrais)
}

# ---------------- #

.lnvrais4 <- function (parameters, xcont, infhist, suphist, nbans, seuil, dist="GEV") {
 nbans <- .datachange(nbans)
 # Case of infhist=suphist=-1
 nbans[infhist==-1] <- nbans[infhist==-1] + 1
 infhist[infhist==-1]<-NA
 suphist[infhist==-1]<-NA
 infhist <- na.omit(infhist)
 suphist <- na.omit(suphist)
 
 # Taking into account only flood estimation intervals
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.GEV(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans-1) * log(F.GEV(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(F.GEV(suphist, xi, alfa, k) - F.GEV(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(dnorm(xcont, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- sum((nbans - 1) * log(pnorm(seuil, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- lnvraishist + sum(log(pnorm(suphist, mean=parameters[1], sd=parameters[2]) -
                                       pnorm(infhist, mean=parameters[1], sd=parameters[2])))
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.exp(xcont, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.exp(seuil, xi, alfa))) +
                 sum(log(F.exp(suphist, xi, alfa) - F.exp(infhist, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.genlogis(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.genlogis(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(F.genlogis(suphist, xi, alfa, k) - F.genlogis(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0 && all(xcont > xi)) {
   lnvraiscont <- sum(log(f.genpar(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0 && all(seuil > xi)) {
   lnvraishist <- sum((nbans - 1) * log(F.genpar(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(infhist-xi)/alfa) > 1)==0 && all(infhist > xi)) {
   lnvraishist <- lnvraishist + sum(log(F.genpar(suphist, xi, alfa, k) - F.genpar(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.gumb(xcont, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.gumb(seuil, xi, alfa))) +
                 sum(log(F.gumb(suphist, xi, alfa) - F.gumb(infhist, xi, alfa)))
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (is.na(alfa)!=TRUE && sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(seuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.lognorm(seuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
  if (is.na(alfa)!=TRUE && sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(F.lognorm(suphist, xi, alfa, k) - F.lognorm(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(dlnorm(xcont, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- sum((nbans - 1) * log(plnorm(seuil, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- lnvraishist + sum(log(plnorm(suphist, meanlog=parameters[1], sdlog=parameters[2]) -
                                       plnorm(infhist, meanlog=parameters[1], sdlog=parameters[2])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  if (alfa > 0 && is.na(beta)!=TRUE) {
   lnvraiscont <- sum(log(f.gamma(xcont, xi, beta, alfa)))
   lnvraishist <- sum((nbans - 1) * log(F.gamma(seuil, xi, beta, alfa))) +
                  sum(log(F.gamma(suphist, xi, beta, alfa) - F.gamma(infhist, xi, beta, alfa)))
  }
  else {
   lnvraiscont <- .thresML
   lnvraishist <- .thresML
  }
 }
 else stop(".lnvrais4(parameters, xcont, infhist, suphist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresML
 if (lnvrais < .thresML) lnvrais <- .thresML

 return(lnvrais)
}



# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.datachange <- function (nbans) {
  # This function changes the nbans vector to fit into the way the likelihhod formulas in BayesianMCMC are written

  nbans_new <- rep(NA, length(nbans))  # new nbans vector built
  numberofyears <- rep(NA, length(nbans[nbans!=0]))   # to store the nbans value associated to each threshold
  numberofzeros <- rep(NA, length(nbans[nbans!=0]))   # to count the number of zeros for each threshold

  n <- 0 # Initialization of the number of successive zeros n
  j <- 1 # Index of numberofyears vector

  for (i in 1:length(nbans)) {        # i: loop on the initial nbans vector

          if(nbans[i]!=0) {
              numberofyears[j] <- nbans[i]
          }
          else {
              n <- n+1     # If nbans[i]=0 we count one more zero
          }

          # stopping the decount
          if(i!=length(nbans)) { # If it is not the last value
              if(nbans[i+1]!=0) {   # and the following value is not a zero => we stop the decount of zeros and store it
                  numberofzeros[j] <- n
                  j <- j+1    # we go to the next index in numberofyears
                  n <- 0   # we reinitialize the decount at first non negative value

              }
              #if nbans[i+1]==0 the count continues
          }
          else { # last i, i.e. last number in nbans
              numberofzeros[j] <- n
          }
  }

  # numberofyear ans numberofzeros are built, now building nbans_new
  index<-1
  for (i in 1:length(numberofyears)) {
            nbans_new[index] <- numberofyears[i] - numberofzeros[i]
      if(numberofzeros[i]!=0) {
          for (j in 1:numberofzeros[i]) {
              nbans_new[index+j] <- 1
          }
      }
      index <- index + numberofzeros[i] + 1
  }

  return(nbans_new)
}
