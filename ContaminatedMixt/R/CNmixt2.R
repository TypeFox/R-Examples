#########################################################
## Model Selection:                                    ##
## Number of mixture components and parsimonious model ##
#########################################################

CNmixt <- function(
  X,                            # matrix of data
  G,                            # vector with number of groups to be evaluated
  model=NULL,                   # models to be considered in model selection
  initialization="mixt",        # initialization procedure: "random.soft", "random.hard", "manual", or "mixt"
  alphafix=NULL,                   # vector of dimension G with proportion of good observations in each group
  alphamin=0.5,                 # vector of minimum proportions of good data 
  etafix=NULL,                     # vector of dimension G with degree of contamination in each group
  etamax=1000,                   # maximum value of eta
  seed=NULL,
  start.z=NULL,                 # (n x k)-matrix of soft or hard classification: it is used only if initialization="manual"    
  start.v=NULL,                 # (n x 2 x k)-array of soft or hard classification in each group: it is used only if initialization="manual"  	
  #veo=FALSE,
  start=0,                      # initialization for the package mixture
  ind.label=NULL,               # indexes of the labelled observations
  label=NULL,                   # groups of the labelled observations
  iter.max=1000,                # maximum number of iterations in the EM-algorithm
  threshold=1.0e-03,            # stopping rule in the Aitken rule
  parallel = FALSE,
  eps=1.0e-100
){
  call=match.call()
  initialization <- match.arg(initialization,c("mixt","kmeans","random.soft","random.hard","manual"))
  if(is.data.frame(X)) 
    X <- as.matrix(X) 
  
  if(nrow(X)<ncol(X)){
    warning("Dimensionality of data exceeds the sample size: it may result in model-fitting failure") 
  }
  
  n <- length(X)
  
  gridG     <- G
  numG      <- length(gridG)
  nummodel  <- length(model)
  
  IC <- array(NA,c(numG,nummodel,3),dimnames=list(paste(gridG,"groups",sep=" "),model,c("loglik","BIC","ICL")))
  
  cont <- 0
  par  <- list()
  
  ###########
  modelXnormNames <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","EEV","VVE","VEV","EVV","VVV")
  eqmod <- data.frame(name=modelXnormNames, number=c(1,1,2,2,2,2,3,3,3,3,3,3,3,3))
  if(is.null(model)) model <- modelXnormNames 
  else model <- match.arg(model, modelXnormNames, several.ok = TRUE)
  
  lm <- list(k=G)
  lm$model <- model
  mm <- expand.grid(lm)
  if (!is.null(model) & 1 %in% G){
    mm <- merge(mm,eqmod,by.x="model",by.y="name")
    mm1 <- mm[mm$k==1,]
    mm2 <- mm1[!duplicated(subset(mm1,select= -model, drop=FALSE)),]
    
    if (nrow(mm1) > nrow(mm2)){
      cat("With G = 1, some models are equivalent, so only one model from each set of equivalent models will be run.\n")
    }
    mm <- rbind(mm2,mm[mm$k!=1,])
  }
  
  mm <- mm[order(mm$k),,drop=FALSE]
  
  job <- function(i){
    cat("\nEstimating model")
    if (!is.null(mm$model[i])) cat(paste0(" ",mm$model[i]))
    cat(paste0(" with G = ",mm$k[i],":"))
     .CNmixtG(
      X=X,  		                      
      G=mm$k[i],                            
      initialization=initialization,      
      modelname=mm$model[i],                     
      alphafix=alphafix,
      alphamin=alphamin, #rep(alphamin,mm$k[i]),
      etafix=etafix,
      etamax=etamax,
      seed=seed,
      start.z=start.z,                 		
      start.v=start.v,                   	
      #veo=veo,
      start=start,                      
      ind.label=ind.label,               
      label=label,                   
      iter.max=iter.max,                
      threshold=threshold,
      eps=eps
    )  
  }
  
  if(parallel){
    cores <- getOption("cl.cores", detectCores())
    cat(paste("\n Using",cores,"cores\n"))
    cl <- makeCluster(cores)
    #clusterExport(cl,envir=environment())
    par <- parLapply(cl=cl,1:nrow(mm),function(i) job(i))
    stopCluster(cl)
  }
  else {
    par <- lapply(1:nrow(mm),function(i) job(i))
  }
  i<- 1
  cat("\n")
  while (!i > length(par)){
    if (! is.null(par[[i]]$error)){
      cat(paste(par[[i]]$error,"\n"))
      par[[i]] <- NULL
    }
    i<- i + 1
  }
  if (!is.null(par)){
    res <-
      structure(
        list(
          call      = call,
          models = par
        ),              
        class = "ContaminatedMixt"
      )
      print(res) 
      invisible(res) }
  else {
    cat("No model was estimated.\n")
    return(NULL)
  }
}

