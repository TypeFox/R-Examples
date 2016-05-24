svmmajcrossval <- function(X,y,search.grid=list(lambda=2^seq(5,-5,length.out=19)),...,convergence=1e-4,weights.obs=1,check.positive=TRUE,
                      print.progress=FALSE,ngroup=5,groups=NULL)
  { 
    return.model <- TRUE 
#PERFORMING CROSS VALIDATION USING A LIST OF PARAMETER SETTINGS
#INITIALISES ALL COMBINATION OF PARAMETER SETTINGS AND
#CALLS THE FUNCTION CROSSVAL.SVMMAJ

    if(!is.factor(y)) y <-factor(y)
    classes      <- sort(levels(y))
    y            <- sign( (y==classes[2]) - .5 )
    n            <-  length(y)

    if(!is.data.frame(X)){
    #OBTAIN THE DIMENSION OF X
        if(!is.matrix(X)) X  <- matrix(X,nrow=n)
        X    <- data.frame(X=X)
    }
    #==================================================
    #INITIALISE THE LEVELS OF ALL PARAMETERS
    #AND TRANSFORM INTO GRIDPOINTS
    #--------------------------------------------------
    params       <- names(search.grid)
    nparam       <- length(search.grid)
    svmmajinput  <- formals(svmmaj.default)
    G            <- expand.grid(search.grid)
    
    #===================================================
    #CHECKING FOR NONNEGATIVE VALUES
    #---------------------------------------------------
    #CHECK FOR NONPOSITIVE VALUES (IF NEEDED)
    if(check.positive){
        if(any(weights.obs<0)) stop('weights should be nonnegative')
        if(any(search.grid[['lambda']]<0))
                               stop('lambda should be nonnegative')
    }
   #====================================================
   #INITIALISING WEIGHTS OF OBJECTS
   #----------------------------------------------------
    if(!is.null(names(weights.obs)) && length(weights.obs)==2)
        weights.obs <- weights.obs[sort(names(weights.obs))]

    w <- rep(1,length(y))
    if(length(weights.obs)==2){
       w[y==-1] <- weights.obs[[1]]
       w[y==1]  <- weights.obs[[2]]
    } else if(length(weights.obs)==length(y))
       w[] <- weights.obs

    #===================================================
    #INITIALIZE GROUP SIZE
    #---------------------------------------------------
    ngroup <- trunc(ngroup)
    if(is.null(groups)){
      if (ngroup < 2) stop("ngroup should be greater than or equal to 2")
      if (ngroup > n) stop("ngroup should be less than or equal to the number of observations")

      if (ngroup == n) {
      #Jackknife cross validation
          groups <- 1:n
          leave.out <- 1
      } else if (ngroup < n) {
      #Fixed number of groups
           leave.out <- trunc(n/ngroup)
           groups    <- sample(rep(1:ngroup, floor(n/ngroup) )[1:(n-leave.out)])
      }
    } else {
        #Prespecified groups
        ngroup    <- max(groups)
        leave.out <- trunc(n/ngroup)
    }
    #===================================================
    #PERFORM CROSS-VALIDATION
    #---------------------------------------------------

     losses  <- rowSums(sapply(1:ngroup,
     #PER GROUP
     function(i) {
        outputs = NULL
        if(print.progress) cat('Group',i,'of',ngroup,':',NROW(G),' gridpoints\n')
        
        formals(svmmaj.default)[params]     <-  as.list(G[1,])
        outputs <- svmmaj.default(X[groups!=i, ], y[groups!=i],weight.obs = w[groups!=i],
                    initial.point = NULL, check.positive=FALSE,convergence=convergence,...)
        lMat <- vector('numeric',NROW(G))  
        #Predict out-of-sample
        qhat    <-  predict.svmmaj(outputs, X[groups==i,],y[groups==i])
        #Calculating number of correctly predicted objects
        lMat[1] <-  sum(w[groups==i]* outputs$hinge(qhat,y[groups==i])$loss)#sum(w[groups==i]*(qhat*y[groups==i]<0 ))
        #PER GRIDPOINT
        if(NROW(G)>1)
          for(j in 2:NROW(G)) {
            #UPDATE INPUT PARAMETER SETTINGS
            outputs$call[params]          <- G[j,]
            outputs$call['initial.point'] <- expression(outputs$theta)
            #UPDATE OUTPUTS
            outputs <- update(outputs)
            if(print.progress){ 
                cat(j,'')
                if(j %% 10 == 0) cat('\n')
            }
            #Predict out-of-sample
            qhat    <-  predict.svmmaj(outputs, X[groups==i,],y[groups==i])
#            gc()
            #Calculating number of correctly predicted objects
            lMat[j] <- sum(w[groups==i]* outputs$hinge(qhat,y[groups==i])$loss)  #sum(w[groups==i]*(qhat*y[groups==i]<0 ))
          }
      if(print.progress) cat('\n')
      return(lMat)
    }))
    #=================================================
    #SAVE THE OPTIMAL VALUES
    #-------------------------------------------------
    index            <- which.min(losses)
    param.opt        <- G[index,,drop=FALSE]
    names(param.opt) <- params
    G                <- cbind(G,loss=losses)
    gc()

    #=================================================
    #RETURN THE RESULTS
    #-------------------------------------------------
    results <- list( param.grid  = G,
                     loss.opt = losses[index],
                     param.opt   = param.opt
                    )
                    
    if(return.model){
      formals(svmmaj.default)[params]     <-  as.list(param.opt) 
      results$model <- svmmaj.default(X, y, weights.obs = w,
                    initial.point = NULL, check.positive=FALSE,convergence=1e-8,...) 
    } 
    return(results)     
}
