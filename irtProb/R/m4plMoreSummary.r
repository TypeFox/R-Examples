`m4plMoreSummary` <-
function(x, out="result", thetaInitial=NULL){
 ## Future add of RMSE if thetaInitial != NULL
 if (!is.null(thetaInitial) && length(thetaInitial) == 1) thetaInitial <- rep(thetaInitial,length(x))
 extract <- function(x,i) x[[i]]
 if (length(x[[1]]$par) > 1) {
  rep                 <- length(x)
  model               <- names(x[[1]]$par)
  personParameters    <- data.frame(t(matrix(unlist(lapply(x,extract,1)),ncol=rep)))
  colnames(personParameters) <- model
  if (!is.null(thetaInitial)) {
   personParameters <- data.frame(personParameters, Error=personParameters[,1] - thetaInitial)
   }
  mPersonParameters   <- sapply(personParameters, mean, na.rm=TRUE)                   ##

  PersonParameters    <- rbind(MEAN   = mPersonParameters,
                               MEDIAN = sapply(personParameters, median, na.rm=TRUE), ##
                               SD     = sapply(personParameters, sd,     na.rm=TRUE), ##
                               N      = colSums(!is.na(personParameters)))
                               
  personSe            <- data.frame(t(matrix(unlist(lapply(x,extract,2)),ncol=rep)))
  colnames(personSe)  <- paste("Se", model, sep="")
  mPersonSe           <- sqrt(sapply(data.frame(personSe^2), mean, na.rm=TRUE))      ## 
  personSeTemp        <- personSe
  if (length(which(is.na(personSe))) > 0) personSeTemp <- personSeTemp[-which(is.na(personSe)),]
  
  PersonSe            <- rbind(MEAN   = mPersonSe,
                               MEDIAN = sapply(personSeTemp, median, na.rm=TRUE), ##
                               SD     = sapply(personSe,     sd, na.rm=TRUE),     ##
                               N      = colSums(!is.na(personSe)))
                               
  colnames(PersonSe)  <- model
  ## Action so that no parameters has a complete zero empirical variance
  wNEZERO <- which(sapply(personParameters,sd) <= 0); if (length(wNEZERO) > 0)    ##
             personParameters[1,wNEZERO] <- personParameters[1,wNEZERO] + 0.00000001
  eCorrelation        <- cor(personParameters, use = "pairwise.complete.obs")
  personCor           <- t(matrix(unlist(lapply(x,extract,3)),ncol=rep))
  nParameters         <- length(model)
  personCor           <- unlist(lapply(x,extract,3))
  personCor           <- matrix(sapply(data.frame(t(matrix(personCor,ncol=rep))), mean, na.rm=TRUE), ##
                         ncol=nParameters)
  ## Correction to be sure that NA values are returned
  personCor[is.nan(c(personCor))] <- NA
  personLL            <- t(matrix(unlist(lapply(x,extract,4)),ncol=rep))
  mPersonLL           <- sapply(data.frame(personLL), mean, na.rm=TRUE)           ##
  personLLTemp        <- data.frame(personLL)
  if (length(which(is.na(personLL))) > 0) personLLTemp <- personLLTemp[-which(is.na(personLL)),]
  
  PersonLL            <- rbind(MEAN   = mPersonLL,
                               MEDIAN = apply(personLLTemp,  2, median, na.rm=TRUE),
                               SD     = sapply(data.frame(personLL), sd, na.rm=TRUE), ##
                               N      = colSums(!is.na(personLL)))
                               
  #apply(personLLTemp,2,median, na.rm=TRUE)
  colnames(PersonLL)  <- colnames(personLL) <- c("LL","AIC","BIC")
  colnames(personCor) <- rownames(personCor) <- model
  if (out == "result") result <- data.frame(data.frame(personParameters),
                                            data.frame(personSe),
                                            data.frame(personLL))

  if (out == "report") report <- list(parameters    = PersonParameters,
                                      se            = PersonSe,
                                      logLikelihood = PersonLL,
                                      eCorrelation  = eCorrelation,
                                      tCorrelation  = personCor)
  }
 if (length(x[[1]]$par) == 1) {
  model             <- "T"
  nSubjects         <- length(x)
  personParameters  <- matrix(unlist(lapply(x,extract,1)), nrow=nSubjects)
  colnames(personParameters) <- model
  if (!is.null(thetaInitial)) {
   personParameters <- data.frame(cbind(personParameters, personParameters[,1] - thetaInitial))
   colnames(personParameters) <- c("T","Error")
   }
  mPersonParameters <- sapply(personParameters, mean, na.rm=TRUE)                       ##
  PersonParameters  <- rbind(MEAN   = mPersonParameters,
                             MEDIAN = sapply(personParameters, median, na.rm=TRUE),     ##
                             SD     = sapply(personParameters, sd,     na.rm=TRUE),     ##
                             N      = colSums(!is.na(personParameters)))
  personSe          <- matrix(unlist(lapply(x,extract,2)), nrow=nSubjects)
  mPersonSe         <- sqrt(sapply(data.frame(personSe^2), mean, na.rm=TRUE))           ##
  PersonSe          <- rbind(MEAN   = mPersonSe,
                             MEDIAN = sapply(data.frame(personSe), median, na.rm=TRUE), ##
                             SD     = sapply(data.frame(personSe), sd, na.rm=TRUE),     ##
                             N      = colSums(!is.na(personSe)))
  personLL          <- t(matrix(unlist(lapply(x,extract,4)), ncol=nSubjects))
  colnames(personLL)<- c("LL","AIC","BIC")
  mPersonLL         <- sapply(data.frame(personLL), mean, na.rm=TRUE)                   ##
  PersonLL          <- rbind(MEAN   = mPersonLL,
                             MEDIAN = sapply(data.frame(personLL), median, na.rm=TRUE), ##
                             SD     = sapply(data.frame(personLL), sd, na.rm=TRUE),     ##
                             N      = colSums(!is.na(personLL)))
  colnames(PersonSe) <- model
  #colnames(personParameters) <- model
  colnames(personSe)         <- paste("Se", model, sep="")

  if (out == "result") result <- data.frame(data.frame(personParameters),
                                            data.frame(personSe),
                                            data.frame(personLL))
  if (out == "report") report <- list(parameters    = PersonParameters,
                                      se            = PersonSe,
                                      logLikelihood = PersonLL,
                                      tCorrelation  = 1,
                                      tCorrelation  = 1)
  }
 if (out == "result") return(result)
 if (out == "report") return(report)
 }


