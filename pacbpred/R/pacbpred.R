pacbpred <-
function(
    niter,
    burnin = floor(niter*2/3),
    Xtrain,
    Xtest,
    Y,
    K = 8,
    cst,
    sigma2,
    alpha = .1,
    delta
    )
  {
    findSeq <- function(x,A,s)
      {
        ind <- NA
        for(i in rev(s))
          {
            vec <- x==A[i:(i+p-1)]
            if(all(vec))
              {
                return(i)
              }
          }
        return(ind)
      }
    dictionary <- function(m, t)
      {
        if(m%%2 == 0)
          {
            return(sin(m*t*pi/4))
          } else {
            return(cos((m + 1)*t*pi/4))
          }
      }
    S <- function(m)
      {
        return(which(m!=0))
      }
    risk <- function(theta,mat)
      {
        if(is.vector(mat))
          {
            r <- t(Y-mat*theta)%*%(Y-mat*theta)/n
          } else {
            r <- t(Y-mat%*%theta)%*%(Y-mat%*%theta)/n
          }
        return(r)
      }
    modelWeight <- function(theta, model, r, post)
      {
        if(sum(abs(theta)) <= cst)
          {
            card1 <- length(S(model))
            card2 <- sum(model)
            part1 <- -log(factorial(card2))
                                        #- lchoose(p, card1) + card2*lcst + log(factorial(card2)) ##CHECK
            part2 <- - delta*r + post
            w <- exp(part1 + part2)
            return(w)
          } else {
            return(0)
          }
      }
    thetaFunc <- function(mat)
      {
        leastSquares <- as.vector(lsfit(mat,Y,intercept=FALSE)$coef)
        theta <- sqrt(sigma2)*rnorm(n = length(leastSquares)) + leastSquares
        return(list(theta = theta, mc = leastSquares))
      }
    ratioFunc <- function(modelO, modelP, rO, rP, postO, postP)
      {
        cardOreg <- length(S(modelO))
        cardOdev <- sum(modelO)
        cardPreg <- length(S(modelP))
        cardPdev <- sum(modelP)
        if(cardOdev > 170 || cardPdev > 170 || useless)
          {
            return(0)
          }
        if(move == 1)
          {
            res1 <- Qreg[3]/Qreg[1]
          }
        if(move == 0)
          {
            res1 <- 1
          }
        if(move == -1)
          {
            res1 <- Qreg[1]/Qreg[3]
          }
        ## res2 <- 1/(weightsB[found]/weightsA[proposedModelIndex]) inverse!
        res2 <- weightsB[found]/weightsA[proposedModelIndex]
        res3part1 <- -log(factorial(cardPdev)) + log(factorial(cardOdev))
                                        #lchoose(p, cardOreg) - lchoose(p, cardPreg) + (cardPdev - cardOdev)*lcst + log(factorial(cardPdev)) - log(factorial(cardOdev)) ##CHECK
        res3part2 <- delta*(rO - rP) + postP - postO
        res3 <- exp(res3part1 + res3part2)
        if(is.nan(res1*res2*res3))
          {
            res <- 1
          } else {
            res <- res1*res2*res3
          }
        return(res)
      }
    if(missing(sigma2))
      {
        sigma2 <- var(Y)
      }
    if(missing(cst))
      {
        cst <- p*(max(Y)-min(Y))
      }
    if(missing(delta))
      {
        delta <- nrow(Xtrain)/(4*sigma2)
      }
    useless <- warn <- FALSE
    bufferSize <- 1000
    saveModels = FALSE
    n <- nrow(Xtrain)
    p <- ncol(Xtrain)
    start <- sample(x = (1:p), size = 1)
    lcst <- log(alpha) - log(2*cst) + log(sqrt(sigma2)) + log(2*pi)/2
    dataFullMatrix <- matrix(nrow = n, ncol = p*K, data = 0)
    for(i in 1:n)
      {
        for(j in 1:p)
          {
            for(k in 1:K)
              {
                dataFullMatrix[i, (j-1)*K + k] <- dictionary(k, Xtrain[i,j])
              }
          }
      }
    for(j in 1:(p*K))
      {
        dataFullMatrix[,j] <- dataFullMatrix[,j] - mean(dataFullMatrix[,j])
      }
    sing <- c()
    V <- matrix(ncol = p*K, nrow = p*K, data = 0)
    for(j in 1:p)
      {
        q <- dataFullMatrix[, ((j - 1)*K + 1):((j - 1)*K + K)]
        sing <- c(sing,svd(q)$d)
        dataFullMatrix[, ((j - 1)*K + 1):((j - 1)*K + K)] <- svd(q)$u
        V[((j - 1)*K + 1):((j - 1)*K + K),((j - 1)*K + 1):((j - 1)*K + K)] <- svd(q)$v
      }
    Sing <- diag(sing)
    models.mcmc <- matrix(nrow = p, ncol = niter, data = 0)
    theta.mcmc <- matrix(nrow = p*K, ncol = niter, data = 0)
    proposedModel.mcmc <- matrix(nrow = p, ncol = niter, data = 0)
    proposedTheta.mcmc <- matrix(nrow = p*K, ncol = niter, data = 0)
    move.mcmc <- numeric(niter)
    ratio.mcmc <- numeric(niter)
    accept <- array(dim = niter)
    overloadBuffer <- 0
    models.mcmc[start,1] <- 1
    for(j in S(models.mcmc[,1]))
      {
        theta.mcmc[(j - 1)*K + 1, 1] <- runif(1)
      }
    proposedModel.mcmc[,1] <- models.mcmc[,1]
    theta.mcmc[,1] <- theta.mcmc[,1]
    move.mcmc[1] <- NA
    ratio.mcmc[1] <- NA
    accept[1] <- NA
    knownModels <- c()
    modelsIndex <- c()
    risques <- numeric(niter)

    cat('MCMC:',sep='')
    for(iiter in 2:niter)
      {
        risques[iiter-1] <- risk(theta.mcmc[,iiter-1],dataFullMatrix)
        if(useless)
          {
            warn <- TRUE
          }
        useless <- FALSE
        cat('... ',floor(100*iiter/niter),'%',sep='')
        if(length(knownModels)>bufferSize)
          {
            knownModels <- c()
            modelsIndex <- c()
            overloadBuffer <- overloadBuffer + 1
          }
        currentModel <- models.mcmc[,iiter-1]
        cardCurrent <- length(S(currentModel))
        Qreg <- Q <- c(.3,.4,.3)
        if(cardCurrent == p || sum(currentModel) > n)
          {
            Q <- c(.9, .1, 0)
          }
        if(cardCurrent == 1)
          {
            Q <- c(0, .5, .5)
          }
        move <- sample(x = c(-1, 0, 1), size = 1, prob = Q)
        move.mcmc[iiter] <- move
        if(move == 1)
          {
            candidatesA <- (1:p)[-S(currentModel)]
            nbNeighborsA <- length(candidatesA)*K
            neighborsA <- matrix(nrow = p, ncol = nbNeighborsA, data = 0)
            neighborThetaA <- matrix(nrow = p*K, ncol = nbNeighborsA, data = 0)
            weightsA <- numeric(nbNeighborsA)
            postA <- numeric(nbNeighborsA)
            rA <- numeric(nbNeighborsA)
            icounter <- 0
            for(j in candidatesA)
              {
                for(k in 1:K)
                  {
                    icounter <- icounter + 1
                    iModel <- c()
                    visitedModel <- currentModel
                    visitedModel[j] <- k
                    neighborsA[,icounter] <- visitedModel
                    for(t in S(visitedModel))
                      {
                        iModel <- c(iModel, ((t - 1)*K + 1):((t - 1)*K + visitedModel[t]))
                      }
                    designMatrixA <- dataFullMatrix[,iModel]
                    #dimension <- dim(designMatrixA)
                    #mat2 <- designMatrixA
                    #returned_data = .C("print_matrix",as.double(designMatrixA),n=as.integer(dimension),MM=as.double(mat2))
                    #mat2 <- returned_data$MM
                    #print(mat2==2*designMatrixA)
                    #print(returned_data$n)
                    if(saveModels)
                      {
                        set <- findSeq(visitedModel,knownModels,modelsIndex)
                      }
                    if(saveModels && !is.na(set))
                      {
                        mc <- knownModels[(set+p):(set+p+sum(visitedModel)-1)]
                        theta <- sqrt(sigma2)*rnorm(n = length(mc)) + mc
                        neighborThetaA[iModel,icounter] <- theta
                      } else {
                        thetaPseudoprior <- thetaFunc(designMatrixA)
                        theta <- thetaPseudoprior$theta
                        neighborThetaA[iModel,icounter] <- theta
                        mc <- thetaPseudoprior$mc
                        modelsIndex <- c(modelsIndex,length(knownModels)+1)
                        knownModels <- c(knownModels,visitedModel,mc)
                      }
                    postA[icounter] <- (t(theta - mc)%*%(theta - mc))/(2*sigma2)
                    rA[icounter] <- risk(theta,designMatrixA)
                    weightsA[icounter] <- modelWeight(theta = theta, model = visitedModel, r = rA[icounter], post = postA[icounter])
                  }
              }
            if(sum(weightsA)==0)
              {
                weightsA <- rep(1/nbNeighborsA,nbNeighborsA)
                useless <- TRUE
              } else {
                weightsA <- weightsA/sum(weightsA)
              }
            proposedModelIndex <- sample(x = 1:icounter, size = 1, prob = weightsA)
            proposedModel <- neighborsA[,proposedModelIndex]
            iproposedModel <- c()
            for(t in S(proposedModel))
              {
                iproposedModel <- c(iproposedModel, ((t - 1)*K + 1):((t - 1)*K + proposedModel[t]))
              }
            proposedModel.mcmc[, iiter] <- proposedModel
            proposedTheta <- neighborThetaA[iproposedModel, proposedModelIndex]
            proposedTheta.mcmc[iproposedModel, iiter] <- proposedTheta
            
            candidatesB <- S(proposedModel)
            nbNeighborsB <- length(candidatesB)
            neighborsB <- matrix(nrow = p, ncol = nbNeighborsB, data = 0)
            neighborThetaB <- matrix(nrow = p*K, ncol = nbNeighborsB, data = 0)
            weightsB <- numeric(nbNeighborsB)
            postB <- numeric(nbNeighborsB)
            rRetour <- numeric(nbNeighborsB)
            icounter <- 0
            for(i in candidatesB)
              {
                icounter <- icounter + 1
                iModel <- c()
                visitedModel <- proposedModel
                visitedModel[i] <- 0
                neighborsB[, icounter] <- visitedModel
                for(t in S(visitedModel))
                  {
                    iModel <- c(iModel, ((t - 1)*K + 1):((t - 1)*K + visitedModel[t]))
                  }
                if(sum(visitedModel == currentModel) == p)
                  {
                    found <- icounter
                    iFoundModel <- iModel
                  }
                designMatrixB <- dataFullMatrix[,iModel]
                if(saveModels)
                  {
                    set <- findSeq(visitedModel,knownModels,modelsIndex)
                  }
                if(saveModels && !is.na(set))
                  {
                    mc <- knownModels[(set+p):(set+p+sum(visitedModel)-1)]
                    theta <- sqrt(sigma2)*rnorm(n = length(mc)) + mc
                    neighborThetaB[iModel,icounter] <- theta
                  } else {
                    thetaPseudoprior <- thetaFunc(designMatrixB)
                    theta <- thetaPseudoprior$theta
                    neighborThetaB[iModel,icounter] <- theta
                    mc <- thetaPseudoprior$mc
                    modelsIndex <- c(modelsIndex,length(knownModels)+1)
                    knownModels <- c(knownModels,visitedModel,mc)
                  }
                postB[icounter] <- (t(theta - mc)%*%(theta - mc))/(2*sigma2)
                rRetour[icounter] <- risk(theta,designMatrixB)
                weightsB[icounter] <- modelWeight(theta = theta, model = visitedModel, r = rRetour[icounter], post = postB[icounter])
              }
            if(sum(weightsB)==0)
              {
                weightsB <- rep(1/nbNeighborsB,nbNeighborsB)
              } else {
                weightsB <- weightsB/sum(weightsB)
              }
            
            ratio <- ratioFunc(modelO = currentModel,
                             modelP = proposedModel,
                             rO = rRetour[found],
                             rP = rA[proposedModelIndex],
                             postO = postB[found],
                             postP = postA[proposedModelIndex])
            ratio <- min(1, ratio)
            ratio.mcmc[iiter] <- ratio
            choice <- rbinom(n = 1, size = 1, prob = ratio)
            if(choice == 1)
              {
                models.mcmc[, iiter] <- proposedModel
                theta.mcmc[iproposedModel, iiter] <- proposedTheta
                accept[iiter] <- TRUE
              } else {
                models.mcmc[, iiter] <- models.mcmc[, iiter - 1]
                theta.mcmc[iFoundModel, iiter] <- neighborThetaB[iFoundModel,found]
                accept[iiter] <- FALSE
              }
          }
        if(move == 0)
          {
            candidatesA <- S(currentModel)
            nbNeighborsA <- length(candidatesA)*K
            neighborsA <- matrix(nrow = p, ncol = nbNeighborsA, data = 0)
            neighborThetaA <- matrix(nrow = p*K, ncol = nbNeighborsA, data = 0)
            weightsA <- numeric(nbNeighborsA)
            weights <- weightsA
            postA <- numeric(nbNeighborsA)
            rA <- numeric(nbNeighborsA)
            icounter <- 0
            for(j in candidatesA)
              {
                for(k in 1:K)
                  {
                    icounter <- icounter + 1
                    iModel <- c()
                    visitedModel <- currentModel
                    visitedModel[j] <- k
                    neighborsA[, icounter] <- visitedModel
                    for(t in S(visitedModel))
                      {
                        iModel <- c(iModel, ((t - 1)*K + 1):((t - 1)*K + visitedModel[t]))
                      }
                    designMatrixA <- dataFullMatrix[,iModel]
                    if(saveModels)
                      {
                        set <- findSeq(visitedModel,knownModels,modelsIndex)
                      }
                    if(saveModels && !is.na(set))
                      {
                        mc <- knownModels[(set+p):(set+p+sum(visitedModel)-1)]
                        theta <- sqrt(sigma2)*rnorm(n = length(mc)) + mc
                        neighborThetaA[iModel,icounter] <- theta
                      } else {
                        thetaPseudoprior <- thetaFunc(designMatrixA)
                        theta <- thetaPseudoprior$theta
                        neighborThetaA[iModel,icounter] <- theta
                        mc <- thetaPseudoprior$mc
                        modelsIndex <- c(modelsIndex,length(knownModels)+1)
                        knownModels <- c(knownModels,visitedModel,mc)
                      }
                    postA[icounter] <- (t(theta - mc)%*%(theta - mc))/(2*sigma2)
                    rA[icounter] <- risk(theta,designMatrixA)
                    weightsA[icounter] <- modelWeight(theta = theta, model = visitedModel, r = rA[icounter], post = postA[icounter])
                    if(sum(visitedModel == currentModel) == p)
                      {
                        found <- icounter
                        weights[found] <- weightsA[icounter]
                        weightsA[icounter] <- 0
                        iFoundModel <- iModel
                      }
                  }
              }
            weights[-found] <- weightsA[-found]
            if(sum(weightsA)==0)
              {
                weightsA <- rep(1/nbNeighborsA,nbNeighborsA)
                useless <- TRUE
              } else {
                weightsA <- weightsA/sum(weightsA)
              }
            proposedModelIndex <- sample(x = 1:icounter, size = 1, prob = weightsA)
            proposedModel <- neighborsA[,proposedModelIndex]
            iproposedModel <- c()
            for(t in S(proposedModel))
              {
                iproposedModel <- c(iproposedModel, ((t - 1)*K + 1):((t - 1)*K + proposedModel[t]))
              }
            proposedModel.mcmc[, iiter] <- proposedModel
            proposedTheta <- neighborThetaA[iproposedModel, proposedModelIndex]
            proposedTheta.mcmc[iproposedModel, iiter] <- proposedTheta
            
            weightsB <- weights
            weightsB[proposedModelIndex] <- 0
            if(sum(weightsB)==0)
              {
                weightsB <- rep(1/nbNeighborsA,nbNeighborsA)
              } else {
                weightsB <- weightsB/sum(weightsB)
              }
            
            ratio <- ratioFunc(modelO = currentModel,
                             modelP = proposedModel,
                             rO = rA[found],
                             rP = rA[proposedModelIndex],
                             postO = postA[found],
                             postP = postA[proposedModelIndex])
            ratio <- min(1, ratio)
            ratio.mcmc[iiter] <- ratio
            choice <- rbinom(n = 1, size = 1, prob = ratio)
            if(choice == 1)
              {
                models.mcmc[, iiter] <- proposedModel
                theta.mcmc[iproposedModel, iiter] <- proposedTheta
                accept[iiter] <- TRUE
              } else {
                models.mcmc[, iiter] <- models.mcmc[, iiter - 1]
                theta.mcmc[iFoundModel, iiter] <- neighborThetaA[iFoundModel,found]
                accept[iiter] <- FALSE
              }
          }
        if(move == -1)
          {
            candidatesA <- S(currentModel)
            nbNeighborsA <- length(candidatesA)
            neighborsA <- matrix(nrow = p, ncol = nbNeighborsA, data = 0)
            neighborThetaA <- matrix(nrow = p*K, ncol = nbNeighborsA, data = 0)
            weightsA <- numeric(nbNeighborsA)
            postA <- numeric(nbNeighborsA)
            rA <- numeric(nbNeighborsA)
            icounter <- 0
            for(j in candidatesA)
              {
                icounter <- icounter + 1
                iModel <- c()
                visitedModel <- currentModel
                visitedModel[j] <- 0
                neighborsA[, icounter] <- visitedModel
                for(t in S(visitedModel))
                  {
                    iModel <- c(iModel, ((t - 1)*K + 1):((t - 1)*K + visitedModel[t]))
                  }
                designMatrixA <- dataFullMatrix[,iModel]
                if(saveModels)
                  {
                    set <- findSeq(visitedModel,knownModels,modelsIndex)
                  }
                if(saveModels && !is.na(set))
                  {
                    mc <- knownModels[(set+p):(set+p+sum(visitedModel)-1)]
                    theta <- sqrt(sigma2)*rnorm(n = length(mc)) + mc
                    neighborThetaA[iModel,icounter] <- theta
                  } else {
                    thetaPseudoprior <- thetaFunc(designMatrixA)
                    theta <- thetaPseudoprior$theta
                    neighborThetaA[iModel,icounter] <- theta
                    mc <- thetaPseudoprior$mc
                    modelsIndex <- c(modelsIndex,length(knownModels)+1)
                    knownModels <- c(knownModels,visitedModel,mc)
                  }
                postA[icounter] <- (t(theta - mc)%*%(theta - mc))/(2*sigma2)
                rA[icounter] <- risk(theta,designMatrixA)
                weightsA[icounter] <- modelWeight(theta = theta, model = visitedModel, r = rA[icounter], post = postA[icounter])
              }
            if(sum(weightsA)==0)
              {
                weightsA <- rep(1/nbNeighborsA,nbNeighborsA)
                useless <- TRUE
              } else {
                weightsA <- weightsA/sum(weightsA)
              }
            proposedModelIndex <- sample(x = 1:icounter, size = 1, prob = weightsA)
            proposedModel <- neighborsA[,proposedModelIndex]
            iproposedModel <- c()
            for(t in S(proposedModel))
              {
                iproposedModel <- c(iproposedModel, ((t - 1)*K + 1):((t - 1)*K + proposedModel[t]))
              }
            proposedModel.mcmc[, iiter] <- proposedModel
            proposedTheta <- neighborThetaA[iproposedModel, proposedModelIndex]
            proposedTheta.mcmc[iproposedModel, iiter] <- proposedTheta
            
            candidatesB <- (1:p)[-S(proposedModel)]
            nbNeighborsB <- length(candidatesB)*K
            neighborsB <- matrix(nrow = p, ncol = nbNeighborsB, data = 0)
            neighborThetaB <- matrix(nrow = p*K, ncol = nbNeighborsB, data = 0)
            weightsB <- numeric(nbNeighborsB)
            postB <- numeric(nbNeighborsB)
            rRetour <- numeric(nbNeighborsB)
            icounter <- 0
            for(j in candidatesB)
              {
                for(k in 1:K)
                  {
                    icounter <- icounter + 1
                    iModel <- c()
                    visitedModel <- currentModel
                    visitedModel[j] <- k
                    neighborsB[,icounter] <- visitedModel
                    for(t in S(visitedModel))
                      {
                        iModel <- c(iModel, ((t - 1)*K + 1):((t - 1)*K + visitedModel[t]))
                      }
                    if(sum(visitedModel == currentModel) == p)
                      {
                        found <- icounter
                        iFoundModel <- iModel
                      }
                    designMatrixB <- dataFullMatrix[,iModel]
                    if(saveModels)
                      {
                        set <- findSeq(visitedModel,knownModels,modelsIndex)
                      }
                if(saveModels && !is.na(set))
                  {
                    mc <- knownModels[(set+p):(set+p+sum(visitedModel)-1)]
                    theta <- sqrt(sigma2)*rnorm(n = length(mc)) + mc
                    neighborThetaB[iModel,icounter] <- theta
                  } else {
                    thetaPseudoprior <- thetaFunc(designMatrixB)
                    theta <- thetaPseudoprior$theta
                    neighborThetaB[iModel,icounter] <- theta
                    mc <- thetaPseudoprior$mc
                    modelsIndex <- c(modelsIndex,length(knownModels)+1)
                    knownModels <- c(knownModels,visitedModel,mc)
                  }
                    postB[icounter] <- (t(theta - mc)%*%(theta - mc))/(2*sigma2)
                    rRetour[icounter] <- risk(theta,designMatrixB)
                    weightsB[icounter] <- modelWeight(theta = theta, model = visitedModel, r = rRetour[icounter], post = postB[icounter])
                  }
              }
            if(sum(weightsB)==0)
              {
                weightsB <- rep(1/nbNeighborsB,nbNeighborsB)
              } else {
                weightsB <- weightsB/sum(weightsB)
              }
            
            ratio <- ratioFunc(modelO = currentModel,
                             modelP = proposedModel,
                             rO = rRetour[found],
                             rP = rA[proposedModelIndex],
                             postO = postB[found],
                             postP = postA[proposedModelIndex])
            ratio <- min(1, ratio)
            ratio.mcmc[iiter] <- ratio
            choice <- rbinom(n = 1, size = 1, prob = ratio)
            if(choice == 1)
              {
                models.mcmc[, iiter] <- proposedModel
                theta.mcmc[iproposedModel, iiter] <- proposedTheta
                accept[iiter] <- TRUE
              } else {
                models.mcmc[, iiter] <- models.mcmc[, iiter - 1]
                theta.mcmc[iFoundModel, iiter] <- neighborThetaB[iFoundModel,found]
                accept[iiter] <- FALSE
              }
          }
      }
    cat('.\n',sep='')
    
    risques[niter] <- risk(theta.mcmc[,niter],dataFullMatrix)
    if(warn)
      {
        print("Many models visited by the Markov chain were discarded. The value of cst is probably too small.")
      }
    predictor <- apply(X = theta.mcmc[, (burnin + 1):niter], MARGIN = 1, FUN = mean)
    predictor.dictionary <- V%*%solve(Sing)%*%predictor
    res <- list()
    if(!missing(Xtest))
      {
        ntest <- nrow(Xtest)
        dataFullMatrixTest <- matrix(nrow = ntest, ncol = p*K, data = 0)
        for(i in 1:ntest)
          {
            for(j in 1:p)
              {
                for(k in 1:K)
                  {
                    dataFullMatrixTest[i, (j-1)*K + k] <- dictionary(k, Xtest[i,j])
                  }
              }
          }
        estimates <- dataFullMatrixTest%*%predictor.dictionary
        res <- c(res, list(predict = estimates))
      }
    res <- c(res, list(estimates = predictor.dictionary, ratio.mcmc = ratio.mcmc, accept = accept, models.mcmc = models.mcmc, risk.mcmc = risques))
    return(res)
    
  }
