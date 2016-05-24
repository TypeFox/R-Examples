## Objects loaded at startup from data/MTM.RData
if(getRversion() >= "2.15.1") globalVariables(c(
                  'MTM', ## Markov Transition Matrices
                  'Ktm', ## Kt limits to choose a matrix from MTM
                  'Ktlim' ## Daily kt range of each matrix.
                  ))
                  
markovG0 <- function(G0dm, solD){
  timeIndex <- index(solD)
  Bo0d <- solD$Bo0d
  Bo0dm <- aggregate(Bo0d, by=as.yearmon, FUN=mean)
  ktm <- G0dm/Bo0dm

  ##Calcula con que matriz debe trabajar para cada mes
  whichMatrix <- findInterval(ktm, Ktm)

  Ktd <- state <- numeric(length(timeIndex))
  state[1] <- 1
  Ktd[1] <- ktm[state[1]]
  for (i in 2:length(timeIndex)){
    iMonth <- month(timeIndex)[i]
    colMonth <- whichMatrix[iMonth]
    rng <- Ktlim[, colMonth]
    classes <- seq(rng[1], rng[2], length=11)
    matMonth <- MTM[(10*colMonth-9):(10*colMonth),]
    ## http://www-rohan.sdsu.edu/~babailey/stat575/mcsim.r
    state[i] <- sample(1:10, size=1, prob=matMonth[state[i-1],])
    Ktd[i] <- runif(1, min=classes[state[i]], max=classes[state[i]+1])
  }
  G0dmMarkov <- aggregate(Ktd * Bo0d, as.yearmon(timeIndex), FUN=mean)
  fix <- na.locf(G0dm/G0dmMarkov, x=as.POSIXct, xout=timeIndex)
  G0d <- Ktd * Bo0d * fix
  G0d
  }
