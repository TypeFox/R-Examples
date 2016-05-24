#' Risk Aversion Task (Sokol-Hessner, 2009)
#' 
#' @description 
#' Hierarchical Bayesian Modeling of the Risk Aversion Task with the following parameters: "rho" (risk aversion), "lambda" (loss aversion), and "tau" (inverse temp).
#' 
#' \strong{MODEL:}
#' Prospect Theory (Sokol-Hessner et al., 2009, PNAS)
#' 
#' @param data A .txt file containing the data to be modeled. Data columns should be labelled as follows: "subjID", "riskyGain", "riskyLoss", and "safeOption". See \bold{Details} below for more information. 
#' @param niter Number of iterations, including warm-up.
#' @param nwarmup Number of iterations used for warm-up only.
#' @param nchain Number of chains to be run.
#' @param ncore Integer value specifying how many CPUs to run the MCMC sampling on. Defaults to 1. 
#' @param nthin Every \code{i == nthin} sample will be used to generate the posterior distribution. Defaults to 1. A higher number can be used when auto-correlation within the MCMC sampling is high.
#' @param inits Character value specifying how the initial values should be generated. Options are "fixed" or "random" or your own initial values.
#' @param indPars Character value specifying how to summarize individual parameters. Current options are: "mean", "median", or "mode".
#' @param saveDir Path to directory where .RData file of model output (\code{modelData}) can be saved. Leave blank if not interested.
#' @param email Character value containing email address to send notification of completion. Leave blank if not interested. 
#' 
#' @return \code{modelData}  A class \code{'hBayesDM'} object with the following components:
#' \describe{
#'  \item{\code{model}}{Character string with the name of the model ("ra_prospect").}
#'  \item{\code{allIndPars}}{\code{'data.frame'} containing the summarized parameter 
#'    values (as specified by \code{'indPars'}) for each subject.}
#'  \item{\code{parVals}}{A \code{'list'} where each element contains posterior samples
#'    over different model parameters. }
#'  \item{\code{fit}}{A class \code{'stanfit'} object containing the fitted model.}
#' }
#' 
#' @importFrom rstan stan rstan_options extract
#' @importFrom modeest mlv
#' @importFrom mail sendmail
#' @importFrom stats median qnorm
#' @importFrom utils read.table
#'
#' @details 
#' This section describes some of the function arguments in greater detail.
#' 
#' \strong{data} should be assigned a character value specifying the full path and name of the file, including the file extension 
#' (e.g. ".txt"), that contains the behavioral data of all subjects of interest for the current analysis. 
#' The file should be a text (.txt) file whose rows represent trial-by-trial observations and columns 
#' represent variables. For the Risk Aversion Task, there should be four columns of data  with the labels 
#' "subjID", "riskyGain", "riskyLoss", and "safeOption". It is not necessary for the columns to be in this 
#' particular order, however it is necessary that they be labelled correctly and contain the information below:
#' \describe{
#'  \item{\code{"subjID"}}{A unique identifier for each subject within data-set to be analyzed.}
#'  \item{\code{"riskyGain"}}{Possible gain outcome of a risky option.}
#'  \item{\code{"riskyLoss"}}{Possible loss outcome of a risky option.}
#'  \item{\code{"safeOption"}}{Guaranteed amount of a safe option.}
#' } 
#' \strong{*}Note: The data.txt file may contain other columns of data (e.g. "Reaction_Time", "trial_number", etc.), but only the data with the column
#' names listed above will be used for analysis/modeling. As long as the columns above are present and labelled correctly,
#' there is no need to remove other miscellaneous data columns.    
#'  
#' \strong{nwarmup} is a numerical value that specifies how many MCMC samples should not be stored upon the 
#' beginning of each chain. For those familiar with Bayesian methods, this value is equivalent to a burn-in sample. 
#' Due to the nature of MCMC sampling, initial values (where the sampling chain begins) can have a heavy influence 
#' on the generated posterior distributions. The \strong{nwarmup} argument can be set to a high number in order to curb the 
#' effects that initial values have on the resulting posteriors.  
#' 
#' \strong{nchain} is a numerical value that specifies how many chains (i.e. independent sampling sequences) should be
#' used to draw samples from the posterior distribution. Since the posteriors are generated from a sampling 
#' process, it is good practice to run multiple chains to ensure that a representative posterior is attained. When
#' sampling is completed, the multiple chains may be checked for convergence with the \code{plot(myModel, type = "trace")}
#' command. The chains should resemble a "furry caterpillar". 
#' 
#' \strong{nthin} is a numerical value that specifies the "skipping" behavior of the MCMC samples being chosen 
#' to generate the posterior distributions. By default, \strong{nthin} is equal to 1, hence every sample is used to 
#' generate the posterior.
#' 
#' @export 
#' 
#' @examples 
#' \dontrun{
#' # Run the model and store results in "output"
#' output <- ra_prospect("example", 2000, 1000, 3, 3)
#' 
#' # Plot the posterior distributions of the hyper-parameters
#' plot(output)
#' 
#' # Show the WAIC and LOOIC model fit estimates 
#' printFit(output)
#' }

ra_prospect <- function(data     = NULL,
                    niter    = 4000, 
                    nwarmup  = 1000, 
                    nchain   = 1,
                    ncore    = 1, 
                    nthin    = 1,
                    inits    = "random",  
                    indPars  = "mean", 
                    saveDir  = NULL,
                    email    = NULL) {
  
  # Path to .stan model file
  modelPath <- system.file("stan", "ra_prospect.stan", package="hBayesDM")
  
  # To see how long computations take
  startTime <- Sys.time()    
  
  # For using example data
  if (data=="example") {
    data <- system.file("extdata", "ra_exampleData.txt", package = "hBayesDM")
  }
  
  # Load data
  if (file.exists(data)) {
    rawdata <- read.table( data, header = T )
  } else {
    stop("** The data file does not exist. Please check it again. **\n  e.g., data = '/MyFolder/SubFolder/dataFile.txt', ... **\n")
  }  
  
  # Individual Subjects
  subjList <- unique(rawdata[,"subjID"]) # list of subjects x blocks
  numSubjs <- length(subjList)           # number of subjects
  
  # Specify the number of parameters and parameters of interest 
  numPars <- 3
  POI     <- c("mu_rho", "mu_lambda", "mu_tau", 
               "rho" , "lambda", "tau",
               "log_lik")
  
  modelName <- "ra_prospect"
  
  # Information for user
  cat("\nModel name = ", modelName, "\n")
  cat("Data file  = ", data, "\n")
  cat("\nDetails:\n")
  cat(" # of chains                   = ", nchain, "\n")
  cat(" # of cores used               = ", ncore, "\n")
  cat(" # of MCMC samples (per chain) = ", niter, "\n")
  cat(" # of burn-in samples          = ", nwarmup, "\n")
  cat(" # of subjects                 = ", numSubjs, "\n")

  ################################################################################
  # THE DATA.  ###################################################################
  ################################################################################
  
  Tsubj <- as.vector( rep( 0, numSubjs ) ) # number of trials for each subject
  
  for ( sIdx in 1:numSubjs )  {
    curSubj     <- subjList[ sIdx ]
    Tsubj[sIdx] <- sum( rawdata$subjID == curSubj )  # Tsubj[N]
  }
  
  maxTrials <- max(Tsubj)
  
  # Information for user continued
  cat(" # of (max) trials per subject = ", maxTrials, "\n\n")
  
  # for multiple subjects
  gain    <- array(0, c(numSubjs, maxTrials) )
  loss    <- array(0, c(numSubjs, maxTrials) )
  cert    <- array(0, c(numSubjs, maxTrials))
  gamble  <- array(0, c(numSubjs, maxTrials))
  
  for (i in 1:numSubjs) {
    curSubj      <- subjList[i]
    useTrials    <- Tsubj[i]
    tmp          <- subset(rawdata, rawdata$subjID == curSubj)
    gain[i, 1:useTrials]    <- tmp[1:useTrials, "gain"]
    loss[i, 1:useTrials]    <- abs(tmp[1:useTrials, "loss"])  # absolute loss amount
    cert[i, 1:useTrials]    <- tmp[1:useTrials, "cert"]
    gamble[i, 1:useTrials]  <- tmp[1:useTrials, "gamble"]
  }
  
  dataList <- list(
    N       = numSubjs,
    T       = maxTrials,
    Tsubj   = Tsubj,
    numPars = numPars,
    gain    = gain,
    loss    = loss,
    cert    = cert,
    gamble  = gamble )
  
  # inits
  if (inits[1] != "random") {
    if (inits[1] == "fixed") {
      inits_fixed <- c(1.0, 1.0, 1.0)
    } else {
      if (length(inits)==numPars) {
        inits_fixed <- inits
      } else {
        stop("Check your inital values!")
      }
    }  
    genInitList <- function() {
      list(
        mu_p     = c( qnorm( inits_fixed[1]/2 ), qnorm( inits_fixed[2]/5 ), qnorm( inits_fixed[3]/5 ) ),
        sigma    = c(1, 1, 1),
        rho_p    = rep(qnorm( inits_fixed[1]/2 ), numSubjs),
        lambda_p = rep(qnorm( inits_fixed[2]/5 ), numSubjs),
        tau_p    = rep(qnorm( inits_fixed[3]/5 ), numSubjs) 
      )
    }
  } else {
    genInitList <- "random"
  }
  
  rstan::rstan_options(auto_write = TRUE)
  if (ncore > 1) {
    numCores <- parallel::detectCores()
    if (numCores < ncore) {
      options(mc.cores = numCores)
      warning("Number of cores specified for parallel computing greater than number of locally available cores. Using all locally available cores.")
    } else {
      options(mc.cores = ncore)
    }
  } else {
    options(mc.cores = 1)
  }
  
  cat("************************************\n")
  cat("** Building a model. Please wait. **\n")
  cat("************************************\n")
  
  # Fit the Stan model
  fit <- rstan::stan(file   = modelPath, 
                     data   = dataList, 
                     pars   = POI,
                     warmup = nwarmup,
                     init   = genInitList, 
                     iter   = niter, 
                     chains = nchain,
                     thin   = nthin)
  
  # Extract the Stan fit object
  parVals <- rstan::extract(fit, permuted=T)
  
  rho    <- parVals$rho
  lambda <- parVals$lambda
  tau    <- parVals$tau
  
  allIndPars <- array(NA, c(numSubjs, numPars))
  allIndPars <- as.data.frame(allIndPars)
  
  for (i in 1:numSubjs) {
    if (indPars=="mean") {
      allIndPars[i, ] <- c( mean(rho[, i]), 
                            mean(lambda[, i]), 
                            mean(tau[, i]) )    
    } else if (indPars=="median") {
      allIndPars[i, ] <- c( median(rho[, i]), 
                            median(lambda[, i]), 
                            median(tau[, i]) )
    } else if (indPars=="mode") {
      allIndPars[i, ] <- c( as.numeric(modeest::mlv(rho[, i], method="shorth")[1]), 
                            as.numeric(modeest::mlv(lambda[, i], method="shorth")[1]),
                            as.numeric(modeest::mlv(tau[, i], method="shorth")[1]) )
    }   
  }  
  
  allIndPars           <- cbind(allIndPars, subjList)
  colnames(allIndPars) <- c("rho", 
                            "lambda", 
                            "tau", 
                            "subjID")

  # Wrap up data into a list
  modelData        <- list(modelName, allIndPars, parVals, fit)
  names(modelData) <- c("model", "allIndPars","parVals","fit")    
  class(modelData) <- "hBayesDM"
  
  # Total time of computations
  endTime  <- Sys.time()
  timeTook <- endTime - startTime
  
  # If saveDir is specified, save modelData as a file. If not, don't save
  # Save each file with its model name and time stamp (date & time (hr & min))
  if (!is.null(saveDir)) {  
    currTime  <- Sys.time()
    currDate  <- Sys.Date()
    currHr    <- substr(currTime, 12, 13)
    currMin   <- substr(currTime, 15, 16)
    timeStamp <- paste0(currDate, "_", currHr, "_", currMin)
    save(modelData, file=file.path(saveDir, paste0(modelName, "_", timeStamp, ".RData"  ) ) )
  }
  
  # Send email to notify user of completion
  if (is.null(email)==F) {
    mail::sendmail(email, paste("model=", modelName, ", fileName = ", data),
                   paste("Check ", getwd(), ". It took ", as.character.Date(round(timeTook,2)), sep="") )
  }
  # Inform user of completion
  cat("\n************************************\n")
  cat("**** Model fitting is complete! ****\n")
  cat("************************************\n")
  
  return(modelData)
}
  
  
