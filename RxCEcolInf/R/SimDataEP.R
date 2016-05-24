##  Function to generate test data from the model, include exit poll
##  library(mvtnorm)
gendata.ep <- function(nprecincts = 175,
                       nrowcat = 3,
                       ncolcat = 3,
                       colcatnames = c("Dem", "Rep", "Abs"),
                       mu0 = c(-.6, -2.05, -1.7, -.2, -1.45, -1.45),
                       rowcatnames = c("bla", "whi", "his", "asi"),
                       alpha = c(.35, .45, .2, .1),
                       housing.seg = 1,
                       nprecincts.ep = 40,
                       samplefrac.ep = 1/14,
                       K0 = NULL,
                       nu0 = 12,
                       Psi0 = NULL,
                       lambda = 1000,
                       dispersion.low.lim = 1,
                       dispersion.up.lim = 1,
                       outfile=NULL,
                       his.agg.bias.vec = c(0,0), # aggbias, c(1.7, -3)
                       HerfInvexp = 3.5,
                       HerfNoInvexp = 3.5,
                       HerfReasexp = 2){
  
  ## ####################################################
  ## Explanation of input variables
  ## nprecicts is the number of precincts in the jurisdiction
  ## nrowcat is the number of rows in the precinct tables
  ## ncolcat is the number of columns in the precinct tables
  ## colcatnames is the vector of names of the precinct table's
  ##       columns; it's length must = ncolcat
  ## mu0 is the prior vector of mus; must have dimension
  ##       nrowcat * (ncolcat -1); the default of
  ##       c(-.6, -2.05, -1.7, -.2, -1.45, -2) gives:
  ##       bla:  lambda_bD ~ 80%, TO_b ~ 45%
  ##       whi:  lambda_wD ~ 20%, TO_w ~ 50%
  ##       his:  lambda_hD ~ 70%, TO-h ~ 40%
  ## rowcatnames is the vectorof names of the precinct table's
  ##       rows; it's length must be >= nrowcat; if >, the first
  ##       nrowcat elements will be used
  ## alpha is the vector of parameters for the dirichlet distribution
  ##       that governs the racial fraction of the precinct rows; it
  ##       determines the expected values only; housing.seg governs
  ##       the dispersion; length(alpha) must be >= nrowcat; if
  ##       >, the first nrowcat elements will be used
  ## housing.seg governs the dispersion of the racial fractions
  ##       of the precincts; the draw of the racial fractions
  ##       is from Diri(housing.seg * alpha)
  ## nprecincts.ep is the number of precincts in the exit
  ##       poll sample
  ## samplefrac.ep is the fraction of voters in each precinct
  ##       that will be sampled
  ## K0 is a length(mu) square matrix, the prior of the mu
  ##       vector distribution in each precinct; if left at NULL,
  ##       it will be set to I_{length(mu)} * .3
  ## nu0 is the degrees of freedom for the Inv-Wis prior from which
  ##       Sigma is drawn
  ## Psi0 is a length(mu) square matrix, the matrix parameter
  ##       of the Inv-Wis from which Sigma is drawn; if left
  ##       at NULL, it will be set to I_{length(mu)} * length(mu)
  ## lambda is the parameter of the Poisson distribution from which
  ##       the number of voters in each precinct will be drawn
  ## outfile is a filepath & filename string; if non-NULL, the output object
  ##       will be saved there; user should specify a name with
  ##       an .Rdata suffix
  ## his.agg.bias.vec is a 2-dimensional vector; if non-null, it will
  ##       cause the white elements of the location parameter for
  ##       the omega_i's to be a function of the fraction his
  ##       in the ith precinct; currently will not work unless
  ##       nrowcat = ncolcat = 3
  ##  HerfInvexp, HerfNoInvexp, HerfReasexp are all exponents
  ##       for Herf-weighted exit poll precinct sampling probs.
  ##       HerfInvexp and HerfReasexp use inverted Herf exponentiation
  ##       of the precinct racial compositions to generate sampling
  ##       weights.  HerfNoInvexp uses not-inverted Herf exponentiation.

  #print("gendata.ep function current as of 19 March 2009")
  
  ## #####################################################
  ## utility functions
  
  "rdirichlet" <- function (n, alpha.use){
    l <- length(alpha.use)
    x <- matrix(rgamma(l * n, alpha.use), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
  }
  
  "rwish" <-
    function(v, S) {
      if (!is.matrix(S))
        S <- matrix(S)
      if (nrow(S) != ncol(S)) {
        stop(message="S not square in rwish().\n")
      }
      if (v < nrow(S)) {
        stop(message="v is less than the dimension of S in rwish().\n")
      }
      p <- nrow(S)
      CC <- chol(S) 
      Z <- matrix(0, p, p)
      diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
      if(p > 1) {
        pseq <- 1:(p-1)
        Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
      }
      return(crossprod(Z %*% CC))
    }
  
  
  "riwish" <-
    function(v, S) {
      return(solve(rwish(v,solve(S))))
    }

  MakeEPObsDataForOnePrecinct <- function(precnum,
                                          dataprecnum.ep,
                                          Herfprob,
                                          Vtrprob,
                                          fnrowcat = nrowcat,
                                          fncolcat = ncolcat,
                                          frowcatnames = rowcatnames.use,
                                          fcolcatnames = colcatnames){
    obsdata.ep <- as.data.frame(matrix(NA, nrow = sum(dataprecnum.ep), ncol = 5))
    colnames(obsdata.ep) <- c("PrecNum", "PrecProb", "VtrProb",
                              "VtrRace", "VtrChoice")
    obsdata.ep[,1] <- precnum
    obsdata.ep[,2] <- Herfprob
    obsdata.ep[,3] <- Vtrprob
    rowloc <- 1
    for(ii in 1:fnrowcat){
      for(jj in 1:fncolcat){
        cellnum <- ((ii-1) * ncolcat) +  jj
        cellcount <- dataprecnum.ep[cellnum]
        #print(c(ii, jj, cellnum, cellcount, rowloc, rowloc + cellcount -1))
        if(cellcount == 0) next
        obsdata.ep[rowloc:(rowloc + cellcount - 1), 4] <- frowcatnames[ii]
        obsdata.ep[rowloc:(rowloc + cellcount - 1), 5] <- fcolcatnames[jj]
        rowloc <- rowloc + cellcount
      }
    }
    return(obsdata.ep)
  }


  
  GenerateEPSample <- function(Herf, Invert = TRUE,
                               comp.list = completedata.list,
                               Nprecs.ep = nprecincts.ep,
                               racecounts =  t(mydata[1:nrowcat,]),
                               sampfrac.ep = samplefrac.ep,
                               Nrow = nrowcat, Ncol = ncolcat,
                               Nprecs = nprecincts,
                               rnames.use = rowcatnames.use,
                               cnames = colcatnames){
     #  make the column names
    returnmat.ep <- matrix(0, nrow = Nprecs, ncol = Nrow * Ncol)
    colnames.ep <- rep(NA, ncol(returnmat.ep))
    counter <- 1
    for(ii in 1:Nrow){
      for(jj in 1:Ncol){
        colnames.ep[counter] <- paste("KK.", rnames.use[ii], ".",
                                      cnames[jj], sep = "")
        counter <- counter + 1
      }
    }
    colnames(returnmat.ep) <- colnames.ep
    
    XXmat.sq <- (racecounts/apply(racecounts, 1, sum))^Herf
    if(Invert){
      Herfs <- 1/apply(XXmat.sq, 1, sum)
    }
    if(!Invert){
      Herfs <- apply(XXmat.sq, 1, sum)
    }
    Herfprobs <- Herfs/sum(Herfs)
    sampprecincts.ep <- sort(sample((1:Nprecs),
                                    Nprecs.ep,
                                    replace = F,
                                    prob = Herfprobs))
  
    #  obtain the samples
    sampsize.ep.vec <- rep(NA, Nprecs.ep) # to store num in-sample vtrs in each prec
    precsampprob.ep.vec <- rep(NA, Nprecs.ep) # to store w/in prec sample probs
  
    for(ii in 1:Nprecs.ep){
      precnum <- sampprecincts.ep[ii]  #  iith in-sample precinct number
      truecountvec <- as.vector(t(comp.list[[precnum]]))  # precinct's true counts
      prectot <- sum(truecountvec)  #  obtain precincts NN_i
      nnsamp <- floor(prectot * sampfrac.ep)  #  number vtrs in-sample in iith in-sample precinct
      sampsize.ep.vec[ii] <- nnsamp # store num vtrs in-sample in iith in-sample precinct
      precsampprob.ep.vec[ii] <- nnsamp/prectot  #  store w/in prec sample prob

      cellvec <- rep(1:(length(truecountvec)), times = truecountvec)
      tempdraw <- sample(cellvec, nnsamp, replace = F)
      draw <- rep(0, length(truecountvec))
      draw[as.numeric(names(table(tempdraw)))] <- table(tempdraw)
      if(min(truecountvec-draw) < 0){
        stop("Something SERIOUSLY wrong with the exit poll drawing function.")
      }
      #  Old method, clumsy and only an approximation
      #draw <- rep(prectot, length(truecountvec))  #  to force the while loop through once
      #while(min(truecountvec - draw) < 0){
      #  draw<- rmultinom(1, nnsamp, truecountvec/prectot) # draw iith in-sample precinct's ep sample
      #}
      returnmat.ep[precnum,] <- draw  #  store draw in big matrix
    }
  
    #  make a dataframe of observed exit polling data only, for use w/ "survey" package
    #      do the first in-sample precinct
    precnumtemp <- sampprecincts.ep[1]
    ObsData.ep <- MakeEPObsDataForOnePrecinct(precnum = precnumtemp,
                                              dataprecnum.ep = returnmat.ep[precnumtemp,],
                                              Herfprob = Herfprobs[precnumtemp],
                                              Vtrprob = precsampprob.ep.vec[1])
    for(qq in 2:Nprecs.ep){
      precnumtemp <- sampprecincts.ep[qq]
      ObsData.ep <- rbind(ObsData.ep,
                          MakeEPObsDataForOnePrecinct(precnum = precnumtemp,
                                                      dataprecnum.ep = returnmat.ep[precnumtemp,],
                                                      Herfprob = Herfprobs[precnumtemp],
                                                      Vtrprob = precsampprob.ep.vec[qq]))
    }
    VtrId <- 1:(nrow(ObsData.ep))
    ObsData.ep <- cbind(ObsData.ep, VtrId)
    list(returnmat.ep = returnmat.ep,
         ObsData.ep = ObsData.ep,
         sampprecincts.ep = sampprecincts.ep)
  }


  ## #####################################################


  ## #####################################################
  ## Peel off values actually to be used, given nrowcat and ncolcat
  rowcatnames.use <- rowcatnames[1:nrowcat]
  if(is.null(K0)){
    K0 <- diag(.3, length(mu0))
  }
  if(is.null(Psi0)){
    Psi0 <- diag(length(mu0))*length(mu0)
  }
  alpha.use <- alpha[1:nrowcat] * housing.seg

  ## #####################################################

  
  
  ## #####################################################
  ## create storage matrices
  ## mydata is marginal total data
  mydata <- matrix(NA, nrowcat+ncolcat+1, nprecincts)
  rownames(mydata) <- c(rowcatnames.use, colcatnames, "nunits")
  colnames(mydata) <- paste("precinct", 1:nprecincts, sep="")

  ## omegamat is the matrix of omega_i values (i on the columns)
  omegamat <- matrix(NA, (ncolcat-1)*nrowcat, nprecincts)
  rownames(omegamat) <- paste(
                              rep(
                                  paste(colcatnames[1:(ncolcat-1)],
                                        colcatnames[ncolcat], sep="."),
                                  nrowcat),
                              rep(rowcatnames.use, rep(ncolcat-1, nrowcat)),
                              sep=".")
  
  
  ## completedata.list is a list of the complete data table
  completedata.list <- vector("list", nprecincts)
  ## #####################################################

  
  ## generate hyper-parameters (Level 3)
  mu <- rmvnorm(1, mean=mu0, sigma=K0)
  Sigma <- riwish(nu0, Psi0) ## may need to be changed depending on
                             ## parameterization of inv. Wishart

  ## add row and column names to parameters
  rownames(K0) <- colnames(K0) <- rownames(omegamat)
  colnames(mu) <- rownames(omegamat)
  rownames(Sigma) <- colnames(Sigma) <- rownames(omegamat)

  for (i in 1:nprecincts){
    
    ## generate number of units in table i
    nunits.i <- rpois(1, lambda * runif(1, dispersion.low.lim, dispersion.up.lim))

    ## generate the row totals
    rowtotal.probs <- rdirichlet(1, alpha.use)
    rowtotals <- rmultinom(1, nunits.i, rowtotal.probs)
    XX <- rowtotals/nunits.i
    names(XX) <- rowcatnames.use

    ## generate Level 2 parameters
    if(sum(abs(his.agg.bias.vec))){
      if(length(mu) != 6) stop("Error:  Cannot yet simulate with Hispanic bias except for 3 x 3 precinct tables.")
      perhis <- XX["his"]
      mu.his.agg.bias <- mu
      mu.his.agg.bias[3] <- mu.his.agg.bias[3] + his.agg.bias.vec[1] * perhis
      mu.his.agg.bias[4] <- mu.his.agg.bias[4] + his.agg.bias.vec[2] * perhis
      omega.i <- rmvnorm(1, mean = mu.his.agg.bias, sigma = Sigma)
    }
    if(sum(abs(his.agg.bias.vec)) == 0){
      omega.i <- rmvnorm(1, mean=mu, sigma=Sigma)
    }

    ## convert omega.i to probabilities (theta.i)
    omega.i.mat <- matrix(omega.i, ncolcat-1, nrowcat)
    omega.i.aug.mat <- rbind(omega.i.mat, 0)
    theta.i.mat <- matrix(NA, nrow(omega.i.aug.mat), ncol(omega.i.aug.mat))
    for (j in 1:ncol(theta.i.mat)){
      theta.i.mat[,j] <- exp(omega.i.aug.mat[,j]) /
        sum( exp( omega.i.aug.mat[,j]))
    }

    

    ## generate the interior of the table
    mytable <- matrix(NA, nrowcat, ncolcat)
    rownames(mytable) <- rowcatnames.use
    colnames(mytable) <- colcatnames
    for (j in 1:nrowcat){
      mytable[j,] <- rmultinom(1, rowtotals[j], theta.i.mat[,j])
    }
    coltotals <- apply(mytable, 2, sum)

    
    ## fill out storage matrices
    mydata[,i] <- c(rowtotals, coltotals, sum(rowtotals))
    omegamat[,i] <- omega.i
    completedata.list[[i]] <- mytable
    
  } ## end i in 1:nprecincts

  
  EPInv <- EPNoInv <- EPReas <- vector("list", 3)
  if(nprecincts.ep){
    EPInv <- GenerateEPSample(Herf = HerfInvexp)
    EPNoInv <- GenerateEPSample(Herf = HerfNoInvexp, Invert = FALSE)
    EPReas <- GenerateEPSample(Herf = HerfReasexp)
  }

  
  ##  put mydata in format corresponding to new Analyze.R formula interface
  mydata <- t(mydata)
  tempcolct <- ncol(mydata)
  mydata <- mydata[,1:(tempcolct-1)]


  
  ####  Create a vector of true values for comparison to put of simulations.
  ####       Along the way, use the names as output in the Analyze.R function.
  ####       Naming is the bigger challenge
#  order of output in current version of Analyze.R:
#     mus
  sim.check.vec <- drop(mu)
  namespos <- 1
  namesvec <- NULL
  for(ii in 1:nrowcat){
    for(jj in 1:(ncolcat-1)){
      newname <- paste("mu.", rowcatnames.use[ii], ".",
                       colcatnames[jj], ".", namespos,
                       sep = "")
      namesvec <- c(namesvec, newname)
      namespos <- namespos+1
    }
  }
  #print("Debug")
  
#     sds of SIGMA
  sim.check.vec <- c(sim.check.vec, sqrt(diag(Sigma)))
  namespos <- 1
  for(ii in 1:nrowcat){
    for(jj in 1:(ncolcat-1)){
      newname <- paste("sd.", rowcatnames.use[ii], ".",
                       colcatnames[jj], ".", namespos,
                       sep = "")
      namesvec <- c(namesvec, newname)
      namespos <- namespos+1
    }
  }
  
#     correlations of SIGMA, all first col, all second col, etc.
  Sigma.cor <- cov2cor(Sigma)
  dim.Sigma.cor <- (dim(Sigma.cor))[1]
  for(ii in 1:(dim.Sigma.cor-1)){
    for(jj in (ii+1):dim.Sigma.cor){
      sim.check.vec <- c(sim.check.vec, Sigma.cor[ii,jj])
      newname <- paste("corr.", ii, ".", jj, sep = "")
      namesvec <- c(namesvec, newname)
    }
  }
  
#     NNs, row major along precinct table
  NNtotals <- matrix(0, nrow = nrowcat, ncol = ncolcat)
  for(ii in 1:nprecincts){
    NNtotals <- NNtotals + completedata.list[[ii]]
  }
  sim.check.vec <- c(sim.check.vec, as.vector(t(NNtotals)))
  for(ii in 1:nrowcat){
    for(jj in 1:ncolcat){
      newname <- paste("NN.", rowcatnames.use[ii], ".",
                       colcatnames[jj], sep = "")
      namesvec <- c(namesvec, newname)
    }
  }
  
#     LAMBDAs, row major along precinct table
  for(ii in 1:nrowcat){
    partrowvec <- NNtotals[ii, -ncolcat]
    LAMBDArow <- partrowvec/sum(partrowvec)
    sim.check.vec <- c(sim.check.vec, LAMBDArow)
  }
  for(ii in 1:nrowcat){
    for(jj in 1:(ncolcat-1)){
      newname <- paste("LAMBDA.", rowcatnames.use[ii], ".",
                       colcatnames[jj],
                       sep = "")
      namesvec <- c(namesvec, newname)
      namespos <- namespos+1
    }
  }
  
  
#     TURNOUT, in obvious order
  for(ii in 1:nrowcat){
    partrowvec <- NNtotals[ii, -ncolcat]
    TURNOUTrow <- sum(partrowvec)/sum(NNtotals[ii,])
    sim.check.vec <- c(sim.check.vec, TURNOUTrow)
  }
  for(ii in 1:nrowcat){
    newname <- paste("TURNOUT.", rowcatnames.use[ii],
                     sep = "")
    namesvec <- c(namesvec, newname)
  }
  
#     GAMMAs, in obvious order
  colsums <- apply(NNtotals, 2, sum)
  TURNOUTtotal <- sum(colsums[-ncolcat])
  for(ii in 1:nrowcat){
    partrowvec <- NNtotals[ii, -ncolcat]
    TURNOUTrowtotal <- sum(partrowvec)
    GAMMArow <- TURNOUTrowtotal/TURNOUTtotal
    sim.check.vec <- c(sim.check.vec, GAMMArow)
  }
  for(ii in 1:nrowcat){
    newname <- paste("GAMMA.", rowcatnames.use[ii],
                     sep = "")
    namesvec <- c(namesvec, newname)
  }
  
  
#     BETAs, row major along precinct table, but skipping last column
  for(ii in 1:nrowcat){
    partrowvec <- NNtotals[ii, -ncolcat]
    BETArow <- partrowvec/sum(NNtotals[ii,])
    sim.check.vec <- c(sim.check.vec, BETArow)
  }
  for(ii in 1:nrowcat){
    for(jj in 1:(ncolcat-1)){
      newname <- paste("BETA.", rowcatnames.use[ii], ".",
                       colcatnames[jj],
                       sep = "")
      namesvec <- c(namesvec, newname)
      namespos <- namespos+1
    }
  }

#     Assign the names
  names(sim.check.vec) <- namesvec
  

  ## create output object
  GQ.sim.object <- list(GQdata=mydata,
                        EPInv = EPInv,
                        EPNoInv = EPNoInv,
                        EPReas = EPReas,
                        HerfInvexp = HerfInvexp,
                        HerfNoInvexp = HerfNoInvexp,
                        HerfReasexp = HerfReasexp,
                        omega.matrix=t(omegamat),
                        interior.tables=completedata.list,
                        mu0=mu0, K0=K0, nu0=nu0, Psi0=Psi0,
                        mu=mu, Sigma=Sigma,
                        Sigma.diag=diag(Sigma),
                        Sigma.cor=Sigma.cor,
                        nprecincts=nprecincts,
                        nrowcat=nrowcat, ncolcat=ncolcat,
                        rowcatnames=rowcatnames.use, colcatnames=colcatnames,
                        lambda=lambda, alpha=alpha.use,
                        sim.check.vec=sim.check.vec)
  
  if (!is.null(outfile)){
    save(GQ.sim.object, file=outfile)
  }
  
  return(GQ.sim.object)

  ## Explanation of GQ.sim.object elements
  ## GQdata is a nprecincts x (nrowcat + ncolcat) matrix, with
  ##      the observed data for the jurisdiction; the first
  ##      nrowcat columns have the precinct table row
  ##      totals, then the column totals
  ## EPInv, EPNoInv, and EPReas are lists of exit polling output.
  ##      Each list is of length 3 and has the following elements.
  ##      ####  "returnmat.ep" is an nprecincts x (nrowcat * ncolcat)
  ##      matrix; for each precinct not in the exit poll
  ##      sample, the corresponding row is a vector of 0s;
  ##      for each precinct in the exit poll sample, the
  ##      corresponding row is a vector of observed exit
  ##      polling counts following the precinct table format
  ##      vectorized ROW major; will be set to all 0s if
  ##      nprecincts.ep = 0
  ##      ####  "ObsData.ep" is a dataframe, one row for each voter
  ##          observed in the exit poll, and 6 columns; columns:
  ##             PrecNum has the precinct number of the observation/voter,
  ##             Precprob has the probability that the precinct was in-sample,
  ##             VtrProb has the probability that the voter was in-sample,
  ##             VtrRace has the voter's exit-poll-reported race,
  ##             VtrChoice has the voter's exit-poll-reported choice
  ##             VtrId has a unique number identifying each voter;
  ##          This matrix is all set to serve as a "data" argument
  ##            in a call to svydesign from the survey package
  ##     ####  "sampprecincts.ep" is a vector of the numbers of the precincts
  ##         in the exit poll sample
  ##  HerfInvexp, HerfNoInvexp, and HerfReasexp are the exponents
  ##     that generated the data in EPInv, EPNoInv, and EPReas,
  ##     as explained in the arguments
  ## omega.matrix is a dim(mu) x nprecincts matrix of each
  ##     precinct's omegas; note the screwy Byron-authored
  ##     dimensionality of the matrix
  ## interior.tables is a list of length nprecincts; each
  ##     element of the list is a nrowcat x ncolcat matrix
  ##     with the internal cell counts
  ## mu0, K0, nu0, and Psi0 are the input parameter values
  ## mu is the drawn mu vector
  ## Sigma is the drawn SIGMA matrix
  ## Sigma.diag is a vector of variances (NOT sds!)
  ## Sigma.cor is a matrix, diagionals = 1, off-diagionals
  ##     are the correlations in Sigma
  ## nprecicts, nrowcat, ncolcat, colcatnames, and lambda
  ##     are all the input values
  ## rowcatnames is a vector of the inputted precinct
  ##     table row names that were actually used
  ## alpha is the vector of inputed alpha values
  ##     that were were actually used
  ## sim.check.vec outputs a bunch of the true parameters
  ##     all in the order that the parameters are outputted
  ##     in the Analyze.R function, so that these can be
  ##     slotted into a coverage simulation matrix

}

if(0){
#  debugging
library(mvtnorm)
set.seed(2598810)

Checkep <- function(EPDraw){
  blacols <- c(1, 2, 3)
  whicols <- c(4, 5, 6)
  hiscols <- c(7, 8, 9)
  Demcols <- c(1, 4, 7)
  Repcols <- c(2, 5, 8)
  Abscols <- c(3, 6, 9)
  CheckEPMat <- matrix(NA, nrow = 175, ncol = 6)
  colnames(CheckEPMat) <- c("bla", "whi", "his",
                            "Dem", "Rep", "Abs")
  CheckEPMat[,"bla"] <- apply(EPDraw[,blacols], 1, sum)
  CheckEPMat[,"whi"] <- apply(EPDraw[,whicols], 1, sum)
  CheckEPMat[,"his"] <- apply(EPDraw[,hiscols], 1, sum)
  CheckEPMat[,"Dem"] <- apply(EPDraw[,Demcols], 1, sum)
  CheckEPMat[,"Rep"] <- apply(EPDraw[,Repcols], 1, sum)
  CheckEPMat[,"Abs"] <- apply(EPDraw[,Abscols], 1, sum)
  CheckEPMat
}


BadVal <- 0
for(ii in 1:100){
  test <- gendata.ep(his.agg.bias.vec = c(1.7, -3),
                   housing.seg = 1.3)
  
  TestEI <- as.matrix(test$GQdata, nrow = nrow(test$GQdata),
                      ncol = ncol(test$GQdata))

  TestEP <- Checkep(test$EPInv$returnmat.ep)
  
  Badness <- TestEI - TestEP
  if(min(Badness) < 0){
    BadVal <- 1
    print("Oh, Shit!!!")
  }
}
BadVal

}

