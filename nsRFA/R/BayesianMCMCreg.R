# The .thresMLreg is a minimum threshold below which the loglikelihood is not considered
# (i.e., if the calculated loglikelihood is lower to .thresMLreg, it is set to .thresMLreg for numerical reasons).
.thresMLreg = -6000


BayesianMCMCreg <- function (xcont, scont, xhist=NA, infhist=NA, suphist=NA, shist=NA, nbans=NA, seuil=NA,
                          nbpas=1000, nbchaines=3, confint=c(0.05, 0.95), dist="GEV",
                          apriori=function(...){1}, parameters0=NA, varparameters0=NA) {

 # This is the main function! Other functions are called in this one and are written below.
 # To know the meaning of the inputs and the outputs you can type help(BayesianMCMCreg).

 reject <- round(nbpas/10)   # number of initial steps that will be excluded from the posterior distribution
 nbpas <- nbpas + reject
 returnperiods <- 10^(seq(0.1, 4, by=.1))
 nonexceedF <- 1 - 1/returnperiods

 # Here, if initial parameters and parameter variances are not provided, a guess is made using the function .chooseparameters0reg (see below)
 if (any(is.na(parameters0))) {
  parameters0_est <- parameters0
  parameters0 <- .chooseparameters0reg(xcont, scont, xhist, infhist, suphist, shist, dist)
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
     if (!is.na(scont) && (length(scont) == 1)){
        scont <- rep(scont, length(xcont))
     }
     if (!is.na(shist) && (length(shist) == 1)){
        shist <- rep(shist, long)
     }
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
 if(all(is.na(c(xcont, scont, xhist, infhist, suphist, shist, seuil)))) {
  stop("BayesianMCMCreg(xcont, scont, xhist, infhist, suphist, shist, nbpas, nbchaines, dist): no input data")
 }
 if(all(is.na(seuil)) && (all(!is.na(xhist)) || all(!is.na(infhist)))) {
   stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): seuil data missing")
 }
 if(all(is.na(nbans)) && (all(!is.na(xhist)) || all(!is.na(infhist)))) {
   stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): nbans data missing")
 }

 # The following if-else selects the type of likelihood depending on the input provided
 if(all(is.na(c(xhist, infhist, suphist, shist, nbans, seuil)))) {
  # Using only the systematic data
  funzionetest <- call(".lnvrais5reg", quote(parameters[1,,j]), quote(xcont), quote(scont), quote(dist))
  funzionecand <- call(".lnvrais5reg", quote(parameterscand), quote(xcont), quote(scont), quote(dist))
 }
 else if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, shist, seuil, nbans)))) {
  # Using censored information, historical flood known (Stedinger and Cohn, Naulet case b)
  funzionetest <- call(".lnvrais1reg", quote(parameters[1,,j]), quote(xcont), quote(scont), quote(xhist), quote(shist), quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais1reg", quote(parameterscand), quote(xcont), quote(scont), quote(xhist), quote(shist), quote(nbans), quote(seuil), quote(dist))
 }
 else if (all(is.na(c(xhist, suphist))) & all(!is.na(c(infhist, shist, seuil, nbans)))) {
  # Using censored information, historical flood unknown (Stedinger and Cohn, Naulet case a)
  funzionetest <- call(".lnvrais2reg", quote(parameters[1,,j]), quote(xcont), quote(scont), quote(infhist), quote(shist), quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais2reg", quote(parameterscand), quote(xcont), quote(scont), quote(infhist), quote(shist), quote(nbans), quote(seuil), quote(dist))
 }
 else if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, shist, seuil, nbans)))) {
  # Taking into account only flood estimation intervals
  funzionetest <- call(".lnvrais4reg", quote(parameters[1,,j]), quote(xcont), quote(scont), quote(infhist), quote(suphist), quote(shist),
                                    quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais4reg", quote(parameterscand), quote(xcont), quote(scont), quote(infhist), quote(suphist), quote(shist),
                                    quote(nbans), quote(seuil), quote(dist))
 }
 else stop("BayesianMCMCreg(xcont, scont, xhist, infhist, suphist, shist, nbpas, nbchaines, dist): inconsistency in input data")



 # Algorithm with repetitions
 # initialisation
 propsaut <- array(data=NA, dim=c(nbpas, nbchaines))
 propsaut_window <- array(data=NA, dim=c(nbpas, nbchaines))

 for (j in 1:nbchaines) {
  parameters[1,,j] <- .parameterscandMODreg(parameters0, varparameters0, dist)   # first step (see .parameterscandMODreg below)
  varparameters[1,,j] <- varparameters0
  vraistest <- eval(funzionetest) + log(apriori(parameters[1,,j]))   # it is a log-likelihhod
  vraisdist[1,j] <- vraistest
  qq[1,,j] <- .quantilesMODreg(F=nonexceedF, parameters=parameters[1,,j], dist=dist)   # first quantiles (see .quantilesMODreg below)
  nbsaut <- 0
  acceptance_binary[1,j] <- 0
  for (i in 2:nbpas) {
   parameterscand <- .parameterscandMODreg(parameters[i-1,,j], varparameters[i-1,,j], dist)   # other steps
   vraiscand <- eval(funzionecand) + log(apriori(parameterscand))   # it is a log-likelihhod
   valtest <- min((exp(vraiscand - vraistest)), 1)
   test <- runif(1)
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
   qq[i,,j] <- .quantilesMODreg(F=nonexceedF, parameters=parameters[i,,j], dist=dist)   # other quantiles

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
 output <- list(xcont=xcont, scont=scont, xhist=xhist, infhist=infhist, suphist=suphist, shist=shist, nbans=nbans, seuil=seuil,
                nbpas=nbpas, nbchaines=nbchaines, dist=dist, confint=confint, apriori=apriori,
                parameters=parameters, quantiles=qq, varparameters=varparameters,
                parameters0=parameters0, varparameters0=varparameters0,
                vraisdist=vraisdist, propsaut=propsaut,
                returnperiods=returnperiods, intervals=intervals,
                parametersML=parametersML, quantilesML=qqML, logML=logML)
 class(output) <- "BayesianMCMCreg"
 return(output)
}




BayesianMCMCregcont <- function (x, nbpas=NA) {
 if(is.na(nbpas)) nbpas <- x$nbpas
 parameters0 <- x$parametersML
 varparameters0 <- rowMeans(apply(x$parameters, c(2,3), var))

 output <- BayesianMCMCreg (xcont=x$xcont, scont=x$scont, xhist=x$xhist, infhist=x$infhist, suphist=x$suphist, shist=x$shist,
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

print.BayesianMCMCreg <- function (x, ...) {
 dummy <- data.frame(cbind(x$quantilesML, t(x$intervals)), row.names=signif(x$returnperiods, 4))
 names(dummy)[1] <- "ML"
 print(dummy)
}


# ----------------------------- #
plotBayesianMCMCreg_surf <- function(x, surf, ask=FALSE, ...) {
  if (ask) {
   op <- par(ask = TRUE)
   on.exit(par(op))
  }
  # Plots flow probabilities for different values of surface
  paramMCMC<-x$parameters[,,1]
  qq<-x$quantiles[,,1]
  T <- x$returnperiods
  F <- 1-1/T
  quantiles=array(data=NA,dim=c(x$nbpas,length(F)))
  quantilesML=rep(NA,times=length(F))
  intervals=array(data=NA,dim=c(2,length(F)))

  #Graph for each surface
  for (i in 1:length(surf)){
    for (j in 1:length(F)){
      q<- surf[i]^paramMCMC[,1]*qq[,j]
      q<-sort(q)
      quantiles[,j]<- q
    }
    intervals[1,]<- quantiles[x$nbpas*0.05,]
    intervals[2,]<- quantiles[x$nbpas*0.95,]
    quantilesML<- surf[i]^x$parametersML[1]*x$quantilesML
    X<-c(intervals[1,10], intervals[1,20], intervals[1,30])
    Y<-c(quantilesML[10], quantilesML[20], quantilesML[30])
    Z<-c(intervals[2,10], intervals[2,20], intervals[2,30])
    plot(T,quantilesML,log="x",main=paste("S=",surf[i],"(km^2)"),xlab="Return Period T (years)",ylab=expression(paste("Discharge(",m^3,"/s)" )),lty=1,lwd =1.6, col=4, type="l",bty = "n")
    axis(1,col=1,cex=0.6,lwd =1.8)
    axis(2,col=1,cex=0.6,lwd =1.8)
    lines(T, intervals[1,], lty=2, col=2,lwd =1)
    lines(T, intervals[2,], lty=2, col=2,lwd =1)
    legend("topleft", legend=c("ML adjustment","5-95% credibility","gauged data", "ungauged extremes"), lty=c(1,2,0,0), pch=c(NA,NA,3,21), col=c(4,2,1,1), cex=0.7, bty = "n",merge = TRUE)
    grid(equilogs=FALSE)
    if (all(is.na(c(x$xhist, x$infhist, x$suphist)))){
       qcont<-x$xcont/x$scont^x$parametersML[1]
       qqcont<- surf[i]^x$parametersML[1]*qcont
       qqhist<-c(NA)
       qqinfhist<-c(NA)
       qqsuphist<-c(NA)
     }
     else if (all(is.na(c(x$infhist, x$suphist))) & all(!is.na(c(x$xhist)))){
       qcont <- x$xcont/x$scont^x$parametersML[1]
       qqcont <- surf[i]^x$parametersML[1]*qcont
       qhist <- x$xhist/x$shist^x$parametersML[1]
       qqhist<- surf[i]^x$parametersML[1]*qhist
       qseuil <- x$seuil/x$shist^x$parametersML[1]
       qqseuil<- surf[i]^x$parametersML[1]*qseuil
       qqinfhist<-c(NA)
       qqsuphist<-c(NA)
     }
     else if (all(is.na(c(x$xhist, x$suphist))) & all(!is.na(c(x$infhist)))){
       qcont <- x$xcont/x$scont^x$parametersML[1]
       qqcont<- surf[i]^x$parametersML[1]*qcont
       qqhist<-c(NA)
       qqsuphist<-c(NA)
       qinfhist<-x$infhist/x$shist^x$parametersML[1]
       qqinfhist<-surf[i]^x$parametersML[1]*qinfhist
       qseuil <- x$seuil/x$shist^x$parametersML[1]
       qqseuil<- surf[i]^x$parametersML[1]*qseuil
     }
     else if (all(is.na(c(x$xhist))) & all(!is.na(c(x$infhist, x$suphist)))){
       qcont <- x$xcont/x$scont^x$parametersML[1]
       qqcont<- surf[i]^x$parametersML[1]*qcont
       qqhist<-c(NA)
       qinfhist<-x$infhist/x$shist^x$parametersML[1]
       qsuphist<-x$suphist/x$shist^x$parametersML[1]
       qqinfhist<-surf[i]^x$parametersML[1]*qinfhist
       qqsuphist<-surf[i]^x$parametersML[1]*qsuphist
       qseuil <- x$seuil/x$shist^x$parametersML[1]
       qqseuil<- surf[i]^x$parametersML[1]*qseuil
     }
     else stop("plotsurf(x, surf) : inconsistency in input data")

     graph <- .pointspos3reg(qqcont, qqhist, qqinfhist, qqsuphist, x$nbans, qqseuil)
     invisible(graph)
  }

  #Graph of the mean, sums up the surfaces
  for (i in 1:length(surf)) {
     mean_site<-mean(x$xcont[x$scont==surf[i]])
     if (i==1) mq<-mean_site
     else mq<-c(mq,mean_site)
  }
  S=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,15000)  #surface values to determine confidence intervals
  qm=array(data=NA,dim=c(x$nbpas,length(S)))
  Qm=array(data=NA,dim=c(3,length(S)))
  q_mean<-paramMCMC[,2] - paramMCMC[,3]*(gamma(1+paramMCMC[,4])-1)/paramMCMC[,4]
  for (i in 1:length(S)) {
      qm[,i]<- sort(q_mean*S[i]^paramMCMC[,1])
      Qm[1,i]=mean(qm[,i])
      Qm[2,i]=qm[x$nbpas*0.05,i]
      Qm[3,i]=qm[x$nbpas*0.95,i]
  }
  plot(S,Qm[1,],log="xy",type="l",lty=1, col=4, xlim=c(1,15000),ylim=c(1,15000), ylab=expression(paste("log[Q] (",m^3,"/s)" )),xlab=expression(paste("log[S] (",km^2,")" )),bty = "n")
  lines(S,Qm[2,],type="l",lty=2,col=2)
  lines(S,Qm[3,],type="l",lty=2,col=2)
  #mqModel<-x$parametersML[2] - x$parametersML[3]*(gamma(1+x$parametersML[4])-1)/x$parametersML[4]
  #mQModel<- mqModel*surf^x$parametersML[1]
  #points(surf, mQModel,pch=1,col=1,cex=1.5 )
  points(surf,mq,cex=1.2,pch=18,col=2)
  axis(1,col=1,cex=0.6,lwd =1.5)
  axis(2,col=1,cex=0.6,lwd =1.5)
  #legend("topleft", legend=c(expression(paste(hat(Q)[mean]," of model")), expression(paste("90% CI of ",hat(Q)[mean])), expression(paste(Q[mean], " of model")), expression(paste(Q[mean], " of gauged series"))), ncol=1, lty=c(1,2,NA,NA), pch=c(NA,NA,1,18), col=c(4,2,1,2), cex=0.8,bty = "n", merge = FALSE)
  legend("topleft", legend=c(expression(paste(hat(Q)[mean]," of model")), expression(paste("90% CI of ",hat(Q)[mean])), expression(paste(Q[mean], " of gauged series"))), ncol=1, lty=c(1,2,NA), pch=c(NA,NA,18), col=c(4,2,2), cex=0.8,bty = "n", merge = FALSE)
  grid(equilogs=FALSE)
}

# ----------------------------- #

plot.BayesianMCMCreg <- function (x, which=1, ask=FALSE, ...) {
 if (ask) {
  op <- par(ask = TRUE)
  on.exit(par(op))
 }
 show <- rep(FALSE, 100)
 show[which] <- TRUE
 if (show[1]) .plotdiagnMCMCreg02(x, ...)
 if (show[2]) .plotdiagnMCMCreg01(x, ...)
 if (show[3]) .plotdiagnMCMCreg04(x, ...)
 if (show[4]) .plotdiagnMCMCreg05(x, ...)
 if (show[5]) .plotdiagnMCMCreg06(x, ...)   # a-priori distribution
}


# --------------------- #

.plotdiagnMCMCreg01 <- function(x, ...) {
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

.plotdiagnMCMCreg02 <- function(x, ...) {
 # Plot of the frequency curve
 qcont <- x$xcont/x$scont^x$parametersML[1]
 qhist <- x$xhist/x$shist^x$parametersML[1]
 qinfhist <- x$infhist/x$shist^x$parametersML[1]
 qsuphist <- x$suphist/x$shist^x$parametersML[1]
 qseuil <- x$seuil/x$shist^x$parametersML[1]
 nbans <- x$nbans
 T <- c(1,1000)
 X <- c(0, 1.3*max(c(qcont, qhist, qinfhist, qsuphist), na.rm=TRUE))
 plot(T, X, type="n", log="x", ...)
 grid(equilogs=FALSE)
 ret <- .pointspos3reg (qcont, qhist, qinfhist, qsuphist, nbans, qseuil)
 lines(x$returnperiods, x$quantilesML)
 lines(x$returnperiods, x$intervals[1,], lty=2)
 lines(x$returnperiods, x$intervals[2,], lty=2)
 invisible(ret)
}

.plotdiagnMCMCreg03 <- function(x, ...) {
 Nsim=10000
 #if(all(is.na(c(x$xhist, x$infhist, x$suphist, x$seuil)))) {
 # # Calcul sur les seules données systèmatiques
 # funzione <- call(".lnvrais5reg", quote(x$parametersML), quote(xcont), quote(dist))
 #}
 #else if (all(is.na(c(x$infhist, x$suphist))) & all(!is.na(c(x$xhist, x$seuil, x$nbans)))) {
 # # Calcul avec info censurée mais débits historiques connus (Stedinger et Cohn, Naulet cas b)
 # funzione <- call(".lnvrais1reg", quote(x$parametersML), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist))
 #}
 #else if (all(is.na(c(x$xhist, x$suphist))) & all(!is.na(c(x$infhist, x$seuil, x$nbans)))) {
 # # Calcul avec info censurée mais débits historiques non connus (Stedinger et Cohn, Naulet cas a)
 # funzione <- call(".lnvrais2reg", quote(x$parametersML), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
 #}
 #else if (all(is.na(c(x$xhist))) & all(!is.na(c(x$infhist, x$suphist, x$seuil, x$nbans)))) {
 # # Calcul avec prise en compte des seuls intervalles d'estimation de débit
 # funzione <- call(".lnvrais4reg", quote(x$parametersML), quote(xcont), quote(infhist), quote(suphist),
 #                                   quote(nbans), quote(seuil), quote(dist))
 #}
}

.plotdiagnMCMCreg04 <- function(x, ...) {
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



.plotdiagnMCMCreg05 <- function(x, ...) {
 # A posteriori distribution of the parameters
 if (length(x$parametersML) == 3) {
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
 else if (length(x$parametersML) == 4) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow=FALSE))
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

  plot(x$parameters[,c(1,4),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(1,4),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[1], x$parametersML[4], pch=19, cex=1.5, col=7)

  plot(x$parameters[,c(2,3),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(2,3),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[2], x$parametersML[3], pch=19, cex=1.5, col=7)

  plot(x$parameters[,c(2,4),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(2,4),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[2], x$parametersML[4], pch=19, cex=1.5, col=7)

  plot(x$parameters[,c(3,4),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(3,4),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[3], x$parametersML[4], pch=19, cex=1.5, col=7)

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



.plotdiagnMCMCreg06 <- function(x, ...) {
 # plot the a-priori distribution of the parameters
 if (length(x$parametersML) == 3) {
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
    par(def.par)
  }
 }
 else if (length(x$parametersML) == 4) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  x1r <- range(x$parameters[,1,1])
  x2r <- range(x$parameters[,2,1])
  x3r <- range(x$parameters[,3,1])
  x4r <- range(x$parameters[,4,1])
  xx1 <- seq(x1r[1], x1r[2], length=20)
  xx2 <- seq(x2r[1], x2r[2], length=20)
  xx3 <- seq(x3r[1], x3r[2], length=20)
  xx4 <- seq(x3r[1], x3r[2], length=20)
  xx1234 <- expand.grid(xx1, xx2, xx3, xx4)
  N <- dim(xx1234)[1]
  zz <- rep(NA, length=N)
  test <- rep(1, length=N)
  for (i in 1:N) {
   zz[i] <- x$apriori(as.numeric(xx1234[i,]))
  }
  if((identical(test,zz)) == TRUE){
     warning("Impossible to plot the a-priori distribution of the parameters : no a-priori information")
  }
  else{
    layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow=FALSE))
    plot(x$parameters[,c(1,2),1], type="n")
    contour(xx1, xx2, tapply(zz, xx1234[,c(1,2)], sum), add=TRUE, col="darkgray")

    plot(x$parameters[,c(1,3),1], type="n")
    contour(xx1, xx3, tapply(zz, xx1234[,c(1,3)], sum), add=TRUE, col="darkgray")

    plot(x$parameters[,c(1,4),1], type="n")
    contour(xx1, xx4, tapply(zz, xx1234[,c(1,4)], sum), add=TRUE, col="darkgray")

    plot(x$parameters[,c(2,3),1], type="n")
    contour(xx2, xx3, tapply(zz, xx1234[,c(2,3)], sum), add=TRUE, col="darkgray")

    plot(x$parameters[,c(2,4),1], type="n")
    contour(xx2, xx4, tapply(zz, xx1234[,c(2,4)], sum), add=TRUE, col="darkgray")

    plot(x$parameters[,c(3,4),1], type="n")
    contour(xx3, xx4, tapply(zz, xx1234[,c(3,4)], sum), add=TRUE, col="darkgray")
  }
  par(def.par)
 }
}





# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

.pointspos3reg <- function (qcont, qhist, qinfhist, qsuphist, nbans, qseuil, ...)
{
 #where, as defined in plotdiagnMCMC02 :
   #qcont <- x$xcont/x$scont^x$parametersML[1]
   #qhist <- x$xhist/x$shist^x$parametersML[1]
   #qinfhist <- x$infhist/x$shist^x$parametersML[1]
   #qsuphist <- x$suphist/x$shist^x$parametersML[1]
   #qseuil <- x$seuil/x$shist^x$parametersML[1]
 if(all(is.na(c(qhist, qinfhist, qsuphist, qseuil)))) {
  # Using only the systematic data (Cunnane plotting position)
  qcont <- sort(qcont, decreasing=TRUE)
  F <- c(1:length(qcont))
  F <- 1 - ((F - 0.4)/(length(qcont) + 1 - 2*0.4))
  T <- 1/(1 - F)
  x <- qcont
  #plot(T, x, log="x", type="n")
  #grid(equilogs = FALSE)
  points(T, x)
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(qinfhist, qsuphist))) & all(!is.na(c(qhist, qseuil, nbans)))) {
  # Using censored information, historical flood known (Stedinger and Cohn, Naulet case b)
  # Calculates the probabilities associated with the thresholds
  unseuil<-unique(qseuil)
  unseuil <- sort(unseuil)
  unnbans<-rep(NA,length(unseuil))
  uncruessup<-rep(NA,length(unseuil))
  for (i in 1:length(unseuil)) {unnbans[i]<-sum(nbans[qseuil==unseuil[i]])
                                uncruessup[i]<-length(qhist[qseuil==unseuil[i] & qhist>=unseuil[i]])
                                }
  nbiter<-100
  toto<-.points1reg(uncruessup, unnbans, unseuil, nbiter)
  if (length(dim(toto))>0) {
    toto <- apply(toto,1,median)
  }
  #Tseuil <- 1/(1-toto)

  unseuil<-c(0,unseuil)
  toto<-c(0,toto)
  
  qhist[qhist<0]<-NA
  qhist<-na.omit(qhist)

  # Calculates the probabilities of (recent and historical) floods
  qcrues<-c(qcont,qhist)
  allcrues <- matrix(data=NA,nrow=length(qcrues),ncol=2)
  isCont <- c(rep(TRUE,length(qcont)), rep(FALSE,length(qhist)))
  allcrues[,1] <- qcrues
  allcrues[,2] <- isCont
  allcrues <- data.frame(allcrues)
  colnames(allcrues)<-c("qcrues","isCont")
  allcrues <- allcrues[order(allcrues$qcrues),]
  qcrues<-allcrues$qcrues

  for (i in 1:length(unseuil)) {
    if (i<length(unseuil)) {
      qcruesseuil<-qcrues[qcrues>=unseuil[i] & qcrues<unseuil[i+1]]
      xxseuil <- 1 + length(qcruesseuil) - rank(qcruesseuil, ties.method="first")
      xxseuil <- toto[i+1]+(toto[i]-toto[i+1])*((xxseuil - 0.4)/(length(qcruesseuil) + 1 - 2*0.4))
      }
    else {
      qcruesseuil<-qcrues[qcrues>=unseuil[i]]
      xxseuil <- 1 + length(qcruesseuil) - rank(qcruesseuil, ties.method="first")
      xxseuil <- 1+(toto[i]-1)*((xxseuil - 0.4)/(length(qcruesseuil) + 1 - 2*0.4))
      }
    if (i==1) {xx<-xxseuil} else {xx<-c(xx,xxseuil)}
  }
  T <- 1/(1 - xx)

  # Plots the graph
  #points(T, xcrues, pch=3)
  qcont<-qcrues[allcrues$isCont == 1]
  Tcont<-rep(NA,length(qcont))
  Tcont<-T[allcrues$isCont == 1]
  #plot(T, qcrues, log="x", type="n") # creates an empty graph
  #grid(equilogs = FALSE)
  points(Tcont, qcont, pch=3)

  qhist<-qcrues[allcrues$isCont == 0]
  Tqhist<-rep(NA,length(qhist))
  Tqhist<-T[allcrues$isCont == 0]
  points(Tqhist, qhist, pch=1)
  invisible(cbind(T,qcrues))
 }
 else if (all(is.na(c(qhist, qsuphist))) & all(!is.na(c(qinfhist, qseuil, nbans)))) {
  # Using censored information, historical flood unknown (Stedinger and Cohn, Naulet case a)
  # Calculates the probabilities associated with the thresholds
  unseuil<-unique(qseuil)
  unseuil <- sort(unseuil)
  unnbans<-rep(NA,length(unseuil))
  uncruessup<-rep(NA,length(unseuil))
  for (i in 1:length(unseuil)) {unnbans[i]<-sum(nbans[qseuil==unseuil[i]])
                                uncruessup[i]<-length(qinfhist[qseuil==unseuil[i] & qinfhist>=unseuil[i]])
                                }
  nbiter<-100
  toto<-.points1reg(uncruessup, unnbans, unseuil, nbiter)
  if (length(dim(toto))>0) {
    toto <- apply(toto,1,median)
  }
  #Tseuil <- 1/(1-toto)

  unseuil<-c(0,unseuil)
  toto<-c(0,toto)
  
  qinfhist[qinfhist<0]<-NA
  qinfhist<-na.omit(qinfhist)

  # Calculates the probabilities of (recent and historical) floods
  qcrues <- c(qcont,qinfhist)
  allcrues <- matrix(data=NA,nrow=length(qcrues),ncol=2)
  isCont <- c(rep(TRUE,length(qcont)), rep(FALSE,length(qinfhist)))
  allcrues[,1] <- qcrues
  allcrues[,2] <- isCont
  allcrues <- data.frame(allcrues)
  colnames(allcrues)<-c("qcrues","isCont")
  allcrues <- allcrues[order(allcrues$qcrues),]
  qcrues<-allcrues$qcrues

  for (i in 1:length(unseuil)) {
    if (i<length(unseuil)) {
      qcruesseuil<-qcrues[qcrues>=unseuil[i] & qcrues<unseuil[i+1]]
      xxseuil <- 1 + length(qcruesseuil) - rank(qcruesseuil, ties.method="first")
      xxseuil <- toto[i+1]+(toto[i]-toto[i+1])*((xxseuil - 0.4)/(length(qcruesseuil) + 1 - 2*0.4))
      }
    else {
      qcruesseuil<-qcrues[qcrues>=unseuil[i]]
      xxseuil <- 1 + length(qcruesseuil) - rank(qcruesseuil, ties.method="first")
      xxseuil <- 1+(toto[i]-1)*((xxseuil - 0.4)/(length(qcruesseuil) + 1 - 2*0.4))
      }
    if (i==1) {xx<-xxseuil} else {xx<-c(xx,xxseuil)}
  }
  T <- 1/(1 - xx)

  # Plots the graph
  #points(T, xcrues, pch=3)
  qcont<-qcrues[allcrues$isCont == 1]
  Tcont<-rep(NA,length(qcont))
  Tcont<-T[allcrues$isCont == 1]
  #plot(T, qcrues, log="x", type="n") # creates an empty graph
  #grid(equilogs = FALSE)
  points(Tcont, qcont, pch=3)

  qinfhist<-qcrues[allcrues$isCont == 0]
  Tinfhist<-rep(NA,length(qinfhist))
  Tinfhist<-T[allcrues$isCont == 0]
  points(Tinfhist, qinfhist, pch=19)
  if(length(qinfhist)!=0){
    segments(Tinfhist, qinfhist, Tinfhist, 10*max(qcrues), lty=3)
  }
  invisible(cbind(T,qcrues))
 }
 else if (all(is.na(c(qhist))) & all(!is.na(c(qinfhist, qsuphist, qseuil, nbans)))) {
  # Taking into account only flood estimation intervals
  # Calculates the probabilities associated with the thresholds
  unseuil<-unique(qseuil)
  unseuil <- sort(unseuil)
  unnbans<-rep(NA,length(unseuil))
  uncruessup<-rep(NA,length(unseuil))
  qmeanhist<-(qinfhist + qsuphist)/2
  for (i in 1:length(unseuil)) {unnbans[i]<-sum(nbans[qseuil==unseuil[i]])
                                uncruessup[i]<-length(qmeanhist[qseuil==unseuil[i] & qmeanhist>=unseuil[i]])
                                }
  nbiter<-100
  toto<-.points1reg(uncruessup, unnbans, unseuil, nbiter)
  if (length(dim(toto))>0) {
    toto <- apply(toto,1,median)
  }
  #Tseuil <- 1/(1-toto)

  unseuil<-c(0,unseuil)
  toto<-c(0,toto)
  
  qinfhist[qinfhist<0]<-NA
  qinfhist<-na.omit(qinfhist)
  qsuphist[qsuphist<0]<-NA
  qsuphist<-na.omit(qsuphist)
  qmeanhist[qmeanhist<0]<-NA
  qmeanhist<-na.omit(qmeanhist)

  # Calculates the probabilities of (recent and historical) floods
  qcrues<-c(qcont,qmeanhist)
  allcrues <- matrix(data=NA,nrow=length(qcrues),ncol=4)
  isCont <- c(rep(TRUE,length(qcont)), rep(FALSE,length(qmeanhist)))
  allcrues[,1] <- qcrues
  allcrues[,2] <- isCont
  allcrues[,3] <- c(rep(NA,length(qcont)), qinfhist)
  allcrues[,4] <- c(rep(NA,length(qcont)), qsuphist)
  allcrues <- data.frame(allcrues)
  colnames(allcrues)<-c("qcrues","isCont","qinfhist","qsuphist")
  allcrues <- allcrues[order(allcrues$qcrues),]
  qcrues<-allcrues$qcrues

  for (i in 1:length(unseuil)) {
    if (i<length(unseuil)) {
      qcruesseuil<-qcrues[qcrues>=unseuil[i] & qcrues<unseuil[i+1]]
      xxseuil <- 1 + length(qcruesseuil) - rank(qcruesseuil, ties.method="first")
      xxseuil <- toto[i+1]+(toto[i]-toto[i+1])*((xxseuil - 0.4)/(length(qcruesseuil) + 1 - 2*0.4))
      }
    else {
      qcruesseuil<-qcrues[qcrues>=unseuil[i]]
      xxseuil <- 1 + length(qcruesseuil) - rank(qcruesseuil, ties.method="first")
      xxseuil <- 1+(toto[i]-1)*((xxseuil - 0.4)/(length(qcruesseuil) + 1 - 2*0.4))
      }
    if (i==1) {xx<-xxseuil} else {xx<-c(xx,xxseuil)}
  }
  T <- 1/(1 - xx)

  # Plots the graph
  #  points(T, xcrues, pch=3)
  qcont<-qcrues[allcrues$isCont == 1]
  Tcont<-rep(NA,length(qcont))
  Tcont<-T[allcrues$isCont == 1]
  #plot(T, qcrues, log="x", type="n") # creates an empty graph
  #grid(equilogs = FALSE)
  points(Tcont, qcont, pch=3)

  qmeanhist<-qcrues[allcrues$isCont == 0]
  Tmeanhist<-rep(NA,length(qmeanhist))
  Tmeanhist<-T[allcrues$isCont == 0]
  points(Tmeanhist, qmeanhist, pch=19)
  qinfhist<-allcrues$qinfhist[allcrues$isCont == 0]
  qsuphist<-allcrues$qsuphist[allcrues$isCont == 0]
  segments(Tmeanhist, qinfhist, Tmeanhist, qsuphist, lty=3)
  invisible(cbind(T,qcrues))
 }
 else stop(".pointspos3reg(qcont, qhist, qinfhist, qsuphist, nbans, qseuil): inconsistency in input data")
}

#-----------------------------------------------------------------------------#
.points1reg <- .points1  # .points1 is in BayesianMCMC.R


# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.chooseparameters0reg <- function (xcont, scont, xhist, infhist, suphist, shist, dist) {   # This function is used to select the initial parameters for the different distributions
  qcont <- xcont/scont^0.7
  qhist <- xhist/shist^0.7
  qinfhist <- infhist/shist^0.7
  qsuphist <- suphist/shist^0.7
  campionissimo <- sort(c(qcont, qhist, qinfhist, qsuphist))
  ll <- as.numeric(Lmoments(campionissimo))
  if (dist=="GEV") {
   parameters0 <- c(0.7, unlist(par.GEV(ll[1], ll[2], ll[4])))
  }
  else if (dist=="NORM") {
   parameters0 <- c(0.7, mean(campionissimo), sd(campionissimo))
  }
  else if (dist=="EXP") {
   parameters0 <- c(0.7, unlist(par.exp(ll[1], ll[2])))
   if(min(qcont) < parameters0[2]) parameters0[2] <- min(qcont)
  }
  else if (dist=="GENLOGIS") {
   parameters0 <- c(0.7, unlist(par.genlogis(ll[1], ll[2], ll[4])))
  }
  else if (dist=="GENPAR") {
   parameters0 <- c(0.7, unlist(par.genpar(ll[1], ll[2], ll[4])))
   if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont)
  }
  else if ((dist=="GUMBEL")||(dist=="EV1")) {
   parameters0 <- c(0.7, unlist(par.gumb(ll[1], ll[2])))
  }
  else if (dist=="KAPPA") {
   parameters0 <- c(0.7, unlist(par.kappa(ll[1], ll[2], ll[4], ll[5])))
  }
  else if ((dist=="LOGNORM")||(dist=="LN3")) {
   parameters0 <- c(0.7, unlist(par.lognorm(ll[1], ll[2], ll[4])))
  }
  else if ((dist=="LN")||(dist=="LN2")) {
   parameters0 <- c(0.7, mean(log(campionissimo)), sd(log(campionissimo)))
  }
  else if ((dist=="P3")||(dist=="GAM")) {
   parameters0 <- c(0.7, unlist(par.gamma(ll[1], ll[2], ll[4])[1:3]))
   if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont) # Positive skewness!!!
  }
  else stop("BayesianMCMCreg(xcont, scont, xhist, infhist, suphist, shist, nbpas, nbchaines, dist): distribution unknown")
  return(parameters0)
}


# ---------------------------- #

.parameterscandMODreg <- function (parameters, varparameters, dist) {
  # Perform a step using different distributions (normal or lognormal) depending on the parameter
  # (essentially if it must be positive or can be also negative)
  if (dist=="GEV") {
   # I know that the scale parameter is positive
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters124 <- rnorm(rep(1,3), mean=parameters[c(1,2,4)], sd=sqrt(varparameters[c(1,2,4)]))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters124[1], parameters124[2], parameter3, parameters124[3])
  }
  else if (dist=="NORM") {
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters12 <- rnorm(rep(1,2), mean=parameters[c(1,2)], sd=sqrt(varparameters[c(1,2)]))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters12[1], parameters12[2], parameter3)
  }
  else if (dist=="EXP") {
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters12 <- rnorm(rep(1,2), mean=parameters[c(1,2)], sd=sqrt(varparameters[c(1,2)]))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters12[1], parameters12[2], parameter3)
  }
  else if (dist=="GENLOGIS") {
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters124 <- rnorm(rep(1,3), mean=parameters[c(1,2,4)], sd=sqrt(varparameters[c(1,2,4)]))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters124[1], parameters124[2], parameter3, parameters124[3])
  }
  else if (dist=="GENPAR") {
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters124 <- rnorm(rep(1,3), mean=parameters[c(1,2,4)], sd=sqrt(varparameters[c(1,2,4)]))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters124[1], parameters124[2], parameter3, parameters124[3])
  }
  else if ((dist=="GUMBEL")||(dist=="EV1")) {
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters12 <- rnorm(rep(1,2), mean=parameters, sd=sqrt(varparameters))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters12[1], parameters12[2], parameter3)
  }
  else if (dist=="KAPPA") {
   parameterscand <- rnorm(rep(1, 5), mean=parameters, sd=sqrt(varparameters))
  }
  else if ((dist=="LOGNORM")||(dist=="LN3")) {
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters124 <- rnorm(rep(1,3), mean=parameters[c(1,2,4)], sd=sqrt(varparameters[c(1,2,4)]))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters124[1], parameters124[2], parameter3, parameters124[3])
  }
  else if ((dist=="LN")||(dist=="LN2")) {
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters12 <- rnorm(rep(1,2), mean=parameters, sd=sqrt(varparameters))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters12[1], parameters12[2], parameter3)
  }
  else if ((dist=="P3")||(dist=="GAM")) {
   LNvarparameters <- log(1 + varparameters[3]/parameters[3]^2)
   LNparameters <- log(parameters[3]) - LNvarparameters/2
   parameters124 <- rnorm(rep(1,3), mean=parameters[c(1,2,4)], sd=sqrt(varparameters[c(1,2,4)]))
   parameter3 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameters124[1], parameters124[2], parameter3, parameters124[3])
  }
  else stop("BayesianMCMCreg(xcont, scont, xhist, infhist, suphist, shist, nbpas, nbchaines, dist): distribution unknown")
  return(parameterscand)
}



# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.quantilesMODreg <- function (F, parameters, dist="GEV") {
 # Calculates the quantiles for a given distribution, a given parameter set and a given non-exceedance probability
 if (dist=="GEV") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  qq <- invF.GEV(F, xi, alfa, k)
 }
 else if (dist=="NORM") {
  qq <- qnorm(F, mean=parameters[2], sd=parameters[3])
 }
 else if (dist=="EXP") {
  xi <- parameters[2]
  alfa <- parameters[3]
  qq <- invF.exp(F, xi, alfa)
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  qq <- invF.genlogis(F, xi, alfa, k)
 }
 else if (dist=="GENPAR") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  qq <- invF.genpar(F, xi, alfa, k)
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  qq <- invF.gumb(F, xi, alfa)
 }
 else if (dist=="KAPPA") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  h <- parameters[5]
  qq <- invF.kappa(F, xi, alfa, k, h)
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  qq <- invF.lognorm(F, xi, alfa, k)
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  qq <- qlnorm(F, meanlog=parameters[2], sdlog=parameters[3])
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[2]
  beta <- parameters[3]
  alfa <- parameters[4]
  qq <- invF.gamma(F, xi, beta, alfa)
 }
 else stop(".quantilesMODreg(F, parameters, dist): distribution unknown")

 return(qq)
}

# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.lnvrais5reg <- function (parameters, xcont, scont, dist="GEV") {
 # Using only the systematic data
 qcont<- xcont/scont^parameters[1]
 if (dist=="GEV") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(F.GEV(qcont*1.01, xi, alfa, k) - F.GEV(qcont*0.99, xi, alfa, k)))
  }
  else lnvrais <- .thresMLreg
 }
 else if (dist=="NORM") {
  lnvrais <- sum(log(pnorm(qcont*1.01, mean=parameters[2], sd=parameters[3]) - pnorm(qcont*0.99, mean=parameters[2], sd=parameters[3])))
 }
 else if (dist=="EXP") {
  xi <- parameters[2]
  alfa <- parameters[3]
  lnvrais <- sum(log(F.exp(qcont*1.01, xi, alfa) - F.exp(qcont*0.99, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(F.genlogis(qcont*1.01, xi, alfa, k) - F.genlogis(qcont*0.99, xi, alfa, k)))
  }
  else lnvrais <- .thresMLreg
 }
 else if (dist=="GENPAR") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0 && all(qcont > xi)) {
   lnvrais <- sum(log(F.genpar(qcont*1.01, xi, alfa, k) - F.genpar(qcont*0.99, xi, alfa, k)))
  }
  else lnvrais <- .thresMLreg
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  lnvrais <- sum(log(F.gumb(qcont*1.01, xi, alfa) - F.gumb(qcont*0.99, xi, alfa)))
 }
 #else if (dist=="KAPPA") {
 # xi <- parameters[2]
 # alfa <- parameters[3]
 # k <- parameters[4]
 # h <- parameters[5]
 # if (sum((k*(qcont-xi)/alfa) > 1)==0) {
 #  lnvrais <- sum(log(f.kappa(qcont, xi, alfa, k, h)))
 # }
 # else lnvrais <- .thresMLreg
 #}
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(F.lognorm(qcont*1.01, xi, alfa, k) - F.lognorm(qcont*0.99, xi, alfa, k)))
  }
  else lnvrais <- .thresMLreg
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvrais <- sum(log(plnorm(qcont*1.01, meanlog=parameters[2], sdlog=parameters[3]) - plnorm(qcont*0.99, meanlog=parameters[2], sdlog=parameters[3])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[2]
  beta <- parameters[3]
  alfa <- parameters[4]
  if (alfa > 0 && is.na(beta)!=TRUE) {
   lnvrais <- sum(log(F.gamma(qcont*1.01, xi, beta, alfa) - F.gamma(qcont*0.99, xi, beta, alfa)))
  }
  else lnvrais <- .thresMLreg
 }
 else stop(".lnvrais5reg(parameters, xcont, scont, dist): distribution unknown")
 if (is.nan(lnvrais)) lnvrais <- .thresMLreg
 if (lnvrais < .thresMLreg) lnvrais <- .thresMLreg

 return(lnvrais)
}





# ---------------- #

.lnvrais1reg <- function (parameters, xcont, scont, xhist, shist, nbans, seuil, dist="GEV") {
 nbans <- .datachange(nbans)   # .datachange is in BayesianMCMC.R

 # Using censored information, historical flood known (Stedinger and Cohn, Naulet case b)
 qcont<-xcont/scont^parameters[1]
 qhist<-xhist/shist^parameters[1]
 qseuil<-seuil/shist^parameters[1]
 
 # Case of xhist=-1
 nbans[xhist==-1] <- nbans[xhist==-1] + 1
 qhist[xhist==-1]<-NA
 qhist <- na.omit(qhist)
  
 if (dist=="GEV") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.GEV(qcont*1.01, xi, alfa, k) - F.GEV(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans-1) * log(F.GEV(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qhist*1.01-xi)/alfa) > 1)==0 && sum((k*(qhist*0.99-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(F.GEV(qhist*1.01, xi, alfa, k) - F.GEV(qhist*0.99, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(dnorm(qcont, mean=parameters[2], sd=parameters[3])))
  lnvraishist <- sum((nbans - 1) * log(pnorm(seuil, mean=parameters[2], sd=parameters[3])))
  lnvraishist <- lnvraishist + sum(log(pnorm(qhist*1.01, mean=parameters[2], sd=parameters[3]) - pnorm(qhist*0.99, mean=parameters[2], sd=parameters[3])))
 }
 else if (dist=="EXP") {
  xi <- parameters[2]
  alfa <- parameters[3]
  lnvraiscont <- sum(log(F.exp(qcont*1.01, xi, alfa) - F.exp(qcont*0.99, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.exp(qseuil, xi, alfa))) + sum(log(F.exp(qhist*1.01, xi, alfa) - F.exp(qhist*0.99, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.genlogis(qcont*1.01, xi, alfa, k) - F.genlogis(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.genlogis(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qhist*1.01-xi)/alfa) > 1)==0 && sum((k*(qhist*0.99-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(F.genlogis(qhist*1.01, xi, alfa, k) - F.genlogis(qhist*0.99, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if (dist=="GENPAR") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0 && all(qcont > xi)) {
   lnvraiscont <- sum(log(F.genpar(qcont*1.01, xi, alfa, k) - F.genpar(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0 && all(qseuil > xi)) {
   lnvraishist <- sum((nbans - 1) * log(F.genpar(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qhist*1.01-xi)/alfa) > 1)==0 && sum((k*(qhist*0.99-xi)/alfa) > 1)==0 && all(qhist > xi)) {
   lnvraishist <- lnvraishist + sum(log(F.genpar(qhist*1.01, xi, alfa, k) - F.genpar(qhist*0.99, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  lnvraiscont <- sum(log(F.gumb(qcont*1.01, xi, alfa) - F.gumb(qcont*0.99, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.gumb(qseuil, xi, alfa)) + sum(log(F.gumb(qhist*1.01, xi, alfa) - F.gumb(qhist*0.99, xi, alfa))))
 }
 #else if (dist=="KAPPA") {
 # xi <- parameters[2]
 # alfa <- parameters[3]
 # k <- parameters[4]
 # h <- parameters[5]
 # lnvraiscont <- sum(log(f.kappa(xcont, xi, alfa, k, h)))
 # lnvraishist <- sum((nbans - 1) * log(F.kappa(seuil, xi, alfa, k, h)) + sum(log(1 - F.kappa(xhist, xi, alfa, k, h))))
 #}
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.lognorm(qcont*1.01, xi, alfa, k) - F.lognorm(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.lognorm(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qhist*1.01-xi)/alfa) > 1)==0 && sum((k*(qhist*0.99-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(F.lognorm(qhist*1.01, xi, alfa, k) - F.lognorm(qhist*0.99, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(plnorm(qcont*1.01, meanlog=parameters[2], sdlog=parameters[3]) - plnorm(qcont*0.99, meanlog=parameters[2], sdlog=parameters[3])))
  lnvraishist <- sum((nbans - 1) * log(plnorm(qseuil, meanlog=parameters[2], sdlog=parameters[3])))
  lnvraishist <- lnvraishist + sum(log(plnorm(qhist*1.01, meanlog=parameters[2], sdlog=parameters[3]) - plnorm(qhist*0.99, meanlog=parameters[2], sdlog=parameters[3])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[2]
  beta <- parameters[3]
  alfa <- parameters[4]
  if (alfa > 0 && is.na(beta)!=TRUE) {
   lnvraiscont <- sum(log(F.gamma(qcont*1.01, xi, beta, alfa) - F.gamma(qcont*0.99, xi, beta, alfa)))
   lnvraishist <- sum((nbans - 1) * log(F.gamma(qseuil, xi, beta, alfa))) +
                  sum(log(F.gamma(qhist*1.01, xi, beta, alfa) - F.gamma(qhist*0.99, xi, beta, alfa)))
  }
  else {
   lnvraiscont <- .thresMLreg
   lnvraishist <- .thresMLreg
  }
 }
 else stop(".lnvrais1reg(parameters, xcont, scont, xhist, shist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresMLreg
 if (lnvrais < .thresMLreg) lnvrais <- .thresMLreg

 return(lnvrais)
}




# ---------------- #

.lnvrais2reg <- function (parameters, xcont, scont, infhist, shist, nbans, seuil, dist="GEV") {
 nbans <- .datachange(nbans)   # .datachange is in BayesianMCMC.R
 
 # Using censored information, historical flood unknown (Stedinger and Cohn, Naulet case a)
 qcont<-xcont/scont^parameters[1]
 qinfhist<-infhist/shist^parameters[1]
 qseuil<-seuil/shist^parameters[1]
 
 # Case of infhist=-1
 nbans[infhist==-1] <- nbans[infhist==-1] + 1
 infhist[infhist==-1]<-NA
 infhist <- na.omit(infhist)
 
 if (dist=="GEV") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.GEV(qcont*1.01, xi, alfa, k) - F.GEV(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans-1) * log(F.GEV(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qinfhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(1 - F.GEV(qinfhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(pnorm(qcont*1.01, mean=parameters[2], sd=parameters[3]) - pnorm(qcont*0.99, mean=parameters[2], sd=parameters[3])))
  lnvraishist <- sum((nbans - 1) * log(pnorm(qseuil, mean=parameters[2], sd=parameters[3])))
  lnvraishist <- lnvraishist + sum(log(1 - pnorm(qinfhist, mean=parameters[2], sd=parameters[3])))
 }
 else if (dist=="EXP") {
  xi <- parameters[2]
  alfa <- parameters[3]
  lnvraiscont <- sum(log(F.exp(qcont*1.01, xi, alfa) - F.exp(qcont*0.99, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.exp(qseuil, xi, alfa))) + sum(log(1 - F.exp(qinfhist, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.genlogis(qcont*1.01, xi, alfa, k) - F.genlogis(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.genlogis(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qinfhist-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(1 - F.genlogis(qinfhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if (dist=="GENPAR") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0 && all(qcont > xi)) {
   lnvraiscont <- sum(log(F.genpar(qcont*1.01, xi, alfa, k) - F.genpar(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0 && all(qseuil > xi)) {
   lnvraishist <- sum((nbans - 1) * log(F.genpar(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qinfhist-xi)/alfa) > 1)==0 && all(qinfhist > xi)) {
   lnvraishist <- lnvraishist + sum(log(1 - F.genpar(qinfhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  lnvraiscont <- sum(log(F.gumb(qcont*1.01, xi, alfa) - F.gumb(qcont*0.99, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.gumb(qseuil, xi, alfa))) + sum(log(1 - F.gumb(qinfhist, xi, alfa)))
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.lognorm(qcont*1.01, xi, alfa, k) - F.lognorm(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.lognorm(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qinfhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(1 - F.lognorm(qinfhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(plnorm(qcont*1.01, meanlog=parameters[2], sdlog=parameters[3]) - plnorm(qcont*0.99, meanlog=parameters[2], sdlog=parameters[3])))
  lnvraishist <- sum((nbans - 1) * log(plnorm(qseuil, meanlog=parameters[2], sdlog=parameters[3])))
  lnvraishist <- lnvraishist + sum(log(1 - plnorm(qinfhist, meanlog=parameters[2], sdlog=parameters[3])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[2]
  beta <- parameters[3]
  alfa <- parameters[4]
  if (alfa > 0 && is.na(beta)!=TRUE) {
   lnvraiscont <- sum(log(F.gamma(qcont*1.01, xi, beta, alfa) - F.gamma(qcont*0.99, xi, beta, alfa)))
   lnvraishist <- sum((nbans - 1) * log(F.gamma(qseuil, xi, beta, alfa))) +
                  sum(log(1 - F.gamma(qinfhist, xi, beta, alfa)))
  }
  else {
   lnvraiscont <- .thresMLreg
   lnvraishist <- .thresMLreg
  }
 }
 else stop(".lnvrais2reg(parameters, xcont, scont, infhist, shist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresMLreg
 if (lnvrais < .thresMLreg) lnvrais <- .thresMLreg

 return(lnvrais)
}

# ---------------- #

.lnvrais4reg <- function (parameters, xcont, scont, infhist, suphist, shist, nbans, seuil, dist="GEV") {
 nbans <- .datachange(nbans)   # .datachange is in BayesianMCMC.R
    
 # Taking into account only flood estimation intervals
 qcont<-xcont/scont^parameters[1]
 qinfhist<-infhist/shist^parameters[1]
 qsuphist<-suphist/shist^parameters[1]
 qseuil<-seuil/shist^parameters[1]
 
 # Case of infhist=suphist=-1
 nbans[infhist==-1] <- nbans[infhist==-1] + 1
 infhist[infhist==-1]<-NA
 suphist[infhist==-1]<-NA
 infhist <- na.omit(infhist)
 suphist <- na.omit(suphist)
 
 if (dist=="GEV") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.GEV(qcont*1.01, xi, alfa, k) - F.GEV(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans-1) * log(F.GEV(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qinfhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(F.GEV(qsuphist, xi, alfa, k) - F.GEV(qinfhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(pnorm(qcont*1.01, mean=parameters[2], sd=parameters[3]) - pnorm(qcont*0.99, mean=parameters[2], sd=parameters[3])))
  lnvraishist <- sum((nbans - 1) * log(pnorm(qseuil, mean=parameters[2], sd=parameters[3])))
  lnvraishist <- lnvraishist + sum(log(pnorm(qsuphist, mean=parameters[2], sd=parameters[3]) -
                                       pnorm(qinfhist, mean=parameters[2], sd=parameters[3])))
 }
 else if (dist=="EXP") {
  xi <- parameters[2]
  alfa <- parameters[3]
  lnvraiscont <- sum(log(F.exp(qcont*1.01, xi, alfa) - F.exp(qcont*0.99, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.exp(qseuil, xi, alfa))) +
                 sum(log(F.exp(qsuphist, xi, alfa) - F.exp(qinfhist, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.genlogis(qcont*1.01, xi, alfa, k) - F.genlogis(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.genlogis(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qinfhist-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(F.genlogis(qsuphist, xi, alfa, k) - F.genlogis(qinfhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if (dist=="GENPAR") {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0 && all(qcont > xi)) {
   lnvraiscont <- sum(log(F.genpar(qcont*1.01, xi, alfa, k) - F.genpar(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0 && all(qseuil > xi)) {
   lnvraishist <- sum((nbans - 1) * log(F.genpar(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qinfhist-xi)/alfa) > 1)==0 && all(qinfhist > xi)) {
   lnvraishist <- lnvraishist + sum(log(F.genpar(qsuphist, xi, alfa, k) - F.genpar(qinfhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  lnvraiscont <- sum(log(F.gumb(qcont*1.01, xi, alfa) - F.gumb(qcont*0.99, xi, alfa)))
  lnvraishist <- sum((nbans - 1) * log(F.gumb(qseuil, xi, alfa))) +
                 sum(log(F.gumb(qsuphist, xi, alfa) - F.gumb(qinfhist, xi, alfa)))
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[2]
  alfa <- parameters[3]
  k <- parameters[4]
  if (is.na(alfa)!=TRUE && sum((k*(qcont*1.01-xi)/alfa) > 1)==0 && sum((k*(qcont*0.99-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(F.lognorm(qcont*1.01, xi, alfa, k) - F.lognorm(qcont*0.99, xi, alfa, k)))
  }
  else lnvraiscont <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qseuil-xi)/alfa) > 1)==0) {
   lnvraishist <- sum((nbans - 1) * log(F.lognorm(qseuil, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
  if (is.na(alfa)!=TRUE && sum((k*(qinfhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(F.lognorm(qsuphist, xi, alfa, k) - F.lognorm(qinfhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresMLreg
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(plnorm(qcont*1.01, meanlog=parameters[2], sdlog=parameters[3]) - plnorm(qcont*0.99, meanlog=parameters[2], sdlog=parameters[3])))
  lnvraishist <- sum((nbans - 1) * log(plnorm(qseuil, meanlog=parameters[2], sdlog=parameters[3])))
  lnvraishist <- lnvraishist + sum(log(plnorm(qsuphist, meanlog=parameters[2], sdlog=parameters[3]) -
                                       plnorm(qinfhist, meanlog=parameters[2], sdlog=parameters[3])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[2]
  beta <- parameters[3]
  alfa <- parameters[4]
  if (alfa > 0 && is.na(beta)!=TRUE) {
   lnvraiscont <- sum(log(F.gamma(qcont*1.01, xi, beta, alfa) - F.gamma(qcont*0.99, xi, beta, alfa)))
   lnvraishist <- sum((nbans - 1) * log(F.gamma(qseuil, xi, beta, alfa))) +
                  sum(log(F.gamma(qsuphist, xi, beta, alfa) - F.gamma(qinfhist, xi, beta, alfa)))
  }
  else {
   lnvraiscont <- .thresMLreg
   lnvraishist <- .thresMLreg
  }
 }
 else stop(".lnvrais4reg(parameters, xcont, scont, infhist, suphist, shist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresMLreg
 if (lnvrais < .thresMLreg) lnvrais <- .thresMLreg

 return(lnvrais)
}

