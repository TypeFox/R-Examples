Estfun <-
    function(
        data,         # Matrix (n,1) or (n,2) with data
        Delta,        # numeric         : Time between observations
        ImputeN = 5,    # Positive integer: How many data points to impute
                                        # between consecutive observations? (5)
        seed,         # Integer: If positive sets seed, else random
                      # set as random
        GibbsN = 1000,       # Positive integer: Number of iterations of
                                        # the Gibbs sampler
        parKnown = list(),#list of named values for the known parameters
        Start = c(0,0,0,0,1,1), # Starting value for the drift
                                        # parameters in the Gibbs sampler
        diffPriorMean,
        diffPriorCovar,
        diffRW = diag(2), # Variance for MH random walk for diffusion
                                        # coefficients
        LatentPathStart, #
        Model = NULL, # Charater, specifying the model. Currently the
                                        # options are 'OU', 'CIR', 'FHN' and 'FHN5'.

        ## Details,      # logical         : Save detailed output from Gibbs
        ##               # sampler (FALSE
        ## SaveEach,     # positive integer: Save every x sample
        ##               # from Gibbs sampler (1)
        driftPriorMean,
        driftPriorCovar,
        driftRW,
        LatentMeanY0 = 0,
        LatentVarY0 = 1,
        RWrhoPaths = 1,
        RWrho2PathPoints = 1){
  anfang <- Sys.time()

  if(!(Model %in% c("OU","FHN","FHN5","CIR"))){
      stop("Model must be either 'OU', 'CIR', 'FHN' or 'FHN5'.")
    }
  ## if(!(pathMH %in% c(1,2))){
  ##     stop("Variable pathMH must be either 1 (independence sampler) or 2 (RW sampler) for latent path.")
  ## }
    data <- as.matrix(data)
    Start <- matrix(Start,nrow=1)
    LatentPathStart <- matrix(LatentPathStart,ncol=1)
    diffPriorMean <- matrix(diffPriorMean,nrow=1)

### Creating 3*g matrix parMat, where g is the number of parameters to
### estimate. First row is parameter number, second is value. Third
### row is either 1 for drift- or 2 for diffusion parameter.
### If parKnown==list(), parMat is 1*1 with value 0.
    Names1 <- tolower(names(parKnown))
    Names2 <- c(paste("drift",1:ShowModels(Model)$Ndrift,sep=""),
                paste("diff", 1:ShowModels(Model)$Ndiff,sep=""))

    if(length(Names1)==0){
        parMat <- matrix(0,ncol=1,nrow=1)
        NamOrder <- "None"
    }
    if(!is.null(Names1) & length(Names1)!=0){
        if(!all(Names1 %in% Names2)){
            stop('parKnown misspecified. Names must be driftX or diffX where X
         is an integer.\n See the function ShowModels.\n')
        }
        if(!all(!duplicated(Names1))){
            stop('parKnown misspecified. Names must be unique.\n')
        }

        tmp <- match(Names1,Names2)
        NamOrder <- Names1[order(tmp)]
        parMat <- matrix(NA,ncol=length(tmp),nrow=3)
        parMat[1,] <- sort(tmp)
        for(i in 1:length(tmp)){
            parMat[2,i] <- parKnown[[NamOrder[i]]]
            parMat[3,i] <- 1+(NamOrder[i] %in% paste("diff",1:ShowModels(Model)$Ndiff,sep=""))
        }
      }

### If neither driftPriorCovar or driftPriorMean are specified:
  if(is.null(driftPriorCovar) & is.null(driftPriorMean)){
### Use infinite variance for prior and sample directly
      MHdrift <- 0
      driftPriorMean <- matrix(rep(1,4),ncol=4)
      driftPriorCovar <- diag(4)
      driftRW <- diag(4)
  } else
### If driftPriorCovar and driftPriorMean are specified and driftRW is
### not:
      if(!is.null(driftPriorCovar) &
         !is.null(driftPriorMean) &
         is.null(driftRW)){
### Use finite variance for prior and sample directly
          MHdrift <- 1
          driftRW <- diag(4)
      } else
### If driftPriorCovar, driftPriorMean and driftRW are all specified
          if(!is.null(driftPriorCovar) &
             !is.null(driftPriorMean) &
             !is.null(driftRW)){
### Use MH to sample drift
              MHdrift <- 2
          } else{
              stop("Please specify how to update drift parameters!")
          }

  A <- .Call("Estfun",
             Delta,
             data,
             ImputeN,
             GibbsN,
             seed,
             Start,
             parMat,
             diffPriorMean,
             diffPriorCovar,
             diffRW,
             LatentPathStart,
             Model,
             driftPriorMean,
             driftPriorCovar,
             driftRW,
             MHdrift,
             LatentMeanY0,
             LatentVarY0,
             RWrhoPaths,
             RWrho2PathPoints,
#             DUP=FALSE,
             PACKAGE = "BIPOD")

  colnames(A$Drift) <- paste("driftpar",1:dim(A$Drift)[2])
  colnames(A$Diff)<- paste("diffpar",1:dim(A$Diff)[2])

    ende <- Sys.time()
    A$Info <- list(Runtime.secs =
                   as.numeric(difftime(ende,anfang,units="secs")),
                   Model = Model,
                   N     = dim(data)[1],
                   Seed  = seed,
                   Coords= paste(c("One coordinate","Two coordinates"),"observed")[A$flag+1],
                   M     = ImputeN,
                   J     = GibbsN,
                   Delta = Delta,
                   Start = Start,
                   ParFixed = NamOrder)

    colnames(A$Info$Start) <- Names2
    A$flag <- NULL
    class(A) <- "BIPOD"
    return(A)
}
