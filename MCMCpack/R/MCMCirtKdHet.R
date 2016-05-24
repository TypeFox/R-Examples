MCMCirtKdHet <-
function(datamatrix, dimensions, item.constraints=list(),
   burnin = 1000, mcmc = 1000, thin=1, verbose = 0, seed = NA,
   alphabeta.start = NA, b0 = 0, B0=0.04, c0=0, d0=0, store.item = FALSE,
   store.ability=TRUE,store.sigma=TRUE, drop.constant.items=TRUE) {
    
    echo.name <- "MCMCirtKdHet"
    
    # translate variable names to MCMCordfactanal naming convention
    x <- as.matrix(datamatrix) 
    factors <- dimensions
    lambda.constraints <- item.constraints
    lambda.start <- alphabeta.start
    l0 <- b0
    L0 <- B0
    sigma.c0 <- c0
    sigma.d0 <- d0
    store.lambda <- store.item
    store.scores <- store.ability
    drop.constantvars <- drop.constant.items

    # extract X and variable names from the model formula and frame 
          
    if (is.matrix(x)){
    	if (drop.constantvars==TRUE){
    	   x.col.var <- apply(x, 2, var, na.rm=TRUE)
        keep.inds <- x.col.var>0
        keep.inds[is.na(keep.inds)] <- FALSE      
        x <- x[,keep.inds]
      }
      X <- as.data.frame(x)
      xvars <- dimnames(X)[[2]]
      xobs <- dimnames(X)[[1]]
      N <- nrow(X)    # number of observations
      K <- ncol(X)    # number of manifest variables
      for (i in 1:K){
        X[,i] <- as.integer(X[,i])
        if (sum(X[,i] != 0 & X[,i] != 1,na.rm=TRUE) != 0){
      	stop("Data must be 0, 1, and NA only.\n")
      	} 
        X[is.na(X[,i]), i] <- -999
      }
      X <- as.matrix(X)           
    }
    else { 
      stop("Please provide data as a matrix.\n")
      }

    ## take care of the case where X has no row names
    if (is.null(xobs)){
      xobs <- 1:N
    }
    
    #check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    
    ## setup constraints on Lambda
    holder <- build.factor.constraints(lambda.constraints, X, K, factors+1)
    Lambda.eq.constraints <- holder[[1]]
    Lambda.ineq.constraints <- holder[[2]]
    X.names <- holder[[3]]
  
    ## setup prior on Lambda
    holder <- form.factload.norm.prior(l0, L0, K, factors+1, X.names)
    Lambda.prior.mean <- holder[[1]]
    Lambda.prior.prec <- holder[[2]]
    Lambda.prior.mean[,1] <- Lambda.prior.mean[,1] * -1 
    
    
   
    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]


    ## Starting values for Lambda
    Lambda <- matrix(0, K, factors+1)
    if (is.na(lambda.start)){# sets Lambda to equality constraints & 0s
      for (i in 1:K){
        for (j in 1:(factors+1)){
          if (Lambda.eq.constraints[i,j]==-999){
            if(Lambda.ineq.constraints[i,j]==0){
              if (j==1){
                  probit.out <- glm(as.factor(X[X[,i]!=-999,i])~1,
                                    family=binomial(link=probit))
                  probit.beta <- coef(probit.out)
                  Lambda[i,j] <- probit.beta[1]
              }
            }
            if(Lambda.ineq.constraints[i,j]>0){
              Lambda[i,j] <- 1.0
            }
            if(Lambda.ineq.constraints[i,j]<0){
              Lambda[i,j] <- -1.0
            }          
          }
          else Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else if (is.matrix(lambda.start)){
      if (nrow(lambda.start)==K && ncol(lambda.start)==(factors+1))
        Lambda  <- lambda.start
      else {
        cat("Starting values not of correct size for model specification.\n")
        stop("Please respecify and call ", echo.name, "() again\n")
      }
    }
    else if (length(lambda.start)==1 && is.numeric(lambda.start)){
      Lambda <- matrix(lambda.start, K, factors+1)
      for (i in 1:K){
        for (j in 1:(factors+1)){
          if (Lambda.eq.constraints[i,j] != -999)
            Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else {
      cat("Starting values neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call ", echo.name, "() again\n")
    }


    ## define holder for posterior sample
    if (store.scores == FALSE && store.lambda == FALSE && store.sigma == FALSE){
      stop("Please specify parameters to be stored.\n")
    }
    else if (store.scores == TRUE && store.lambda == FALSE && store.sigma == FALSE){
      sample <- matrix(data=0, mcmc/thin, (factors+1)*N)
    }
    else if(store.scores == FALSE && store.lambda == TRUE && store.sigma == FALSE) {
      sample <- matrix(data=0, mcmc/thin, K*(factors+1))
    }
    else if(store.scores==TRUE && store.lambda==TRUE && store.sigma == FALSE) {
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+(factors+1)*N)
    } 
    else if (store.scores == FALSE && store.lambda == FALSE && store.sigma == TRUE){
      sample <- matrix(data=0, mcmc/thin, N)
    }
    else if (store.scores == TRUE && store.lambda == FALSE && store.sigma == TRUE){
      sample <- matrix(data=0, mcmc/thin, (factors+1)*N + N)
    }
    else if(store.scores == FALSE && store.lambda == TRUE && store.sigma == TRUE) {
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+N)
    }
    else if(store.scores==TRUE && store.lambda==TRUE && store.sigma == TRUE) {
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+(factors+1)*N + N)
    }


    ## create templates
    posterior <- .C("irtKdHetpost",
                    samdata = as.double(sample),
                    samrow = as.integer(nrow(sample)),
                    samcol = as.integer(ncol(sample)),
                    Xdata = as.integer(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    Lamstartdata = as.double(Lambda),
                    Lamstartrow = as.integer(nrow(Lambda)),
                    Lamstartcol = as.integer(ncol(Lambda)),
                    Lameqdata = as.double(Lambda.eq.constraints),
                    Lameqrow = as.integer(nrow(Lambda.eq.constraints)),
                    Lameqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lamineqdata = as.double(Lambda.ineq.constraints),
                    Lamineqrow = as.integer(nrow(Lambda.ineq.constraints)),
                    Lamineqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lampmeandata = as.double(Lambda.prior.mean),
                    Lampmeanrow = as.integer(nrow(Lambda.prior.mean)),
                    Lampmeancol = as.integer(ncol(Lambda.prior.prec)),
                    Lampprecdata = as.double(Lambda.prior.prec),
                    Lampprecrow = as.integer(nrow(Lambda.prior.prec)),
                    Lamppreccol = as.integer(ncol(Lambda.prior.prec)),
                    storelambda = as.integer(store.lambda),
                    storescores = as.integer(store.scores),
                    storesigma = as.integer(store.sigma),
                    sigmapriorc = as.double(sigma.c0),
                    sigmapriord = as.double(sigma.d0),
                    package="MCMCpack"
                    )

	output <- mcmc(data=matrix(posterior$samdata, posterior$samrow, posterior$samcol,byrow=FALSE),start=1, end=mcmc, thin=thin)
	
	    par.names <- NULL
    if (store.lambda==TRUE){
   
        alpha.hold <- paste("alpha", X.names, sep=".")
        beta.hold <- paste("beta", X.names, sep = ".")
        beta.hold <- rep(beta.hold, factors, each=factors)
        beta.hold <- paste(beta.hold, 1:factors, sep=".")
                
        Lambda.names <- t(cbind(matrix(alpha.hold, K, 1), 
                                matrix(beta.hold,K,factors,byrow=TRUE)))  
        dim(Lambda.names) <- NULL
      
      par.names <- c(par.names, Lambda.names)
    }
    
    if (store.scores==TRUE){
      phi.names <- paste(paste("theta",
                               rep(xobs, each=(factors+1)), sep="."),
                         rep(0:factors,(factors+1)), sep=".")
      par.names <- c(par.names, phi.names)      
      
    }
    
        if (store.sigma==TRUE){
      sigma.names <- paste("sigma", rep(xobs), sep=".")
      par.names <- c(par.names, sigma.names)
    }
    
    varnames(output) <- par.names 

    # get rid of columns for constrained parameters
    output.df <- as.data.frame(as.matrix(output))

    output.var <- sd(output.df)
    output.df <- output.df[,output.var != 0]
    output <- mcmc(as.matrix(output.df), start=burnin+1, end=burnin+mcmc,
                   thin=thin)
    
    # add constraint info so this isn't lost
    attr(output, "constraints") <- lambda.constraints
    attr(output, "n.manifest") <- K
    attr(output, "n.factors") <- factors
    attr(output, "title") <-
        "MCMCpack Heteroskedastic K-Dimensional Item Response Theory Model Posterior Sample"

    return(output)

    
  }

