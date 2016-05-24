synth <-
function(           data.prep.obj = NULL,
                      X1 = NULL,
                      X0 = NULL,
                      Z0 = NULL,
                      Z1 = NULL,
                      custom.v = NULL,
                      optimxmethod = c("Nelder-Mead","BFGS"),
                      genoud = FALSE,
                      quadopt = "ipop",
                      Margin.ipop = 0.0005,
                      Sigf.ipop = 5,
                      Bound.ipop = 10,
                      verbose = FALSE,
                       ...
                      )
  { 
    # Retrieve dataprep objects
    if(is.null(data.prep.obj) == FALSE)
      {
        cat("\nX1, X0, Z1, Z0 all come directly from dataprep object.\n\n")
        X1 <- data.prep.obj$X1
        Z1 <- data.prep.obj$Z1
        X0 <- data.prep.obj$X0
        Z0 <- data.prep.obj$Z0
      } else {
        cat("X1,X0,Z1,Z0 were individually input (not dataprep object.)\n\n")
      }

     # routine checks
     store <- list(X1=X1,X0=X0,Z1=Z1,Z0=Z0)
     for(i in 1:4){
     if(is.null(store[[i]]))
      {stop(paste("\n",names(store)[i],"is missing \n"))}
     if(sum(is.na(store[[i]]))>0)
      {stop(paste("\n NAs in",names(store)[i],"\n"))}
     if(is.matrix(store[[i]]) == FALSE)
      {stop(paste("\n",names(store)[i],"is not a matrix object\n"))}
     }

    # geometry checks
    if(ncol(X1)!=1){stop("\n Please specify only one treated unit: X1 has to have ncol= 1")}
    if(ncol(Z1)!=1){stop("\n Please specify only one treated unit: Z1 has to have ncol= 1")}

    if(ncol(X0)<2){stop("\n Please specify at least two control units: X0 has to have ncol >= 2 ")}
    if(ncol(Z0)<2){stop("\n Please specify only one treated unit: Z0 has to have ncol >= 2")}

    if(nrow(Z0)!=nrow(Z1)){stop("\n Different number of periods for treated and controls: nrow(Z0) unequal nrow(Z1)")}
    if(nrow(X0)!=nrow(X1)){stop("\n Different number of predictors for treated and controls: nrow(X0) unequal nrow(X1)")}

    if(nrow(X0)==0){stop("No predictors specified. Please specify at least on predictor")}
    if(nrow(Z0)==0){stop("No periods specified for Z1 and Z0. Please specify at least on period")}

    if(0 %in% apply(X0,1,sd))
     {stop("\n At least one predictor in X0 has no variation across control units. Please remove this predictor.")}

    # collinearity check
   # check <- try(solve(t(X0)%*%X0),silent=TRUE)
   # if(class(check)=="try-error")
   #  {stop("\n Some of the predictors in X0 are collinear (t(X0)%*%X0) not invertible")}

    # Normalize X
    nvarsV <- dim(X0)[1]
    big.dataframe <- cbind(X0, X1)
    divisor <- sqrt(apply(big.dataframe, 1, var))
    scaled.matrix <-
      t(t(big.dataframe) %*% ( 1/(divisor) *
                              diag(rep(dim(big.dataframe)[1], 1)) ))

    X0.scaled <- scaled.matrix[,c(1:(dim(X0)[2]))]
    if(is.vector(X0.scaled)==TRUE)
     {X0.scaled <- t(as.matrix(X0.scaled))}
    X1.scaled <- scaled.matrix[,dim(scaled.matrix)[2]]


    # check if custom v weights are supplied or
    # if only on predictor is specified,
    # we jump to quadratic optimization over W weights
    # if not start optimization over V and W
    if(is.null(custom.v) & nrow(X0) != 1)
      {
      
      # two attemps for best V are made: 
      # equal weights and regression based starting values 
     cat("\n****************",
         "\n searching for synthetic control unit  \n","\n"
        )
      
      if(genoud == TRUE) # if user wants genoud as well
      {
      require(rgenoud)                # we run genoud first
      cat("\n****************",
         "\n genoud() requested for optimization\n","\n"
        )

      rgV.genoud <- genoud(
                             fn.V, 
                             nvarsV, 
                             X0.scaled = X0.scaled,
                             X1.scaled = X1.scaled,
                             Z0 = Z0,
                             Z1 = Z1,
                             quadopt = quadopt,
                             margin.ipop = Margin.ipop,
                             sigf.ipop = Sigf.ipop,
                             bound.ipop = Bound.ipop
                             )
      SV1 <- rgV.genoud$par  # and use these as starting values

      cat("\n****************",
         "\n genoud() finished, now running local optimization using optim()\n","\n"
        )

      } else {
      # if we don't use genoud first: set of starting values: equal weights
      SV1 <- rep(1/nvarsV,nvarsV)
      }
      
      # now we run optimization
      all.methods <- FALSE
      if(sum(optimxmethod %in% c("All"))==1){ all.methods <- TRUE }
     rgV.optim.1 <- optimx(par=SV1, fn=fn.V,
                             gr=NULL, hess=NULL, 
                             method=optimxmethod, itnmax=NULL, hessian=FALSE,
                             control=list(kkt=FALSE,
                                          starttests=FALSE,
                                          dowarn=FALSE,
                                          all.methods=all.methods),
                             X0.scaled = X0.scaled,
                             X1.scaled = X1.scaled,
                             Z0 = Z0,
                             Z1 = Z1,
                             quadopt = quadopt,
                             margin.ipop = Margin.ipop,
                             sigf.ipop = Sigf.ipop,
                             bound.ipop = Bound.ipop
                            )
      # get minimum
      if(verbose==TRUE){print(rgV.optim.1)}
      rgV.optim.1 <- collect.optimx(rgV.optim.1,"min")
      
      # second set of starting values: regression method 
      # will sometimes not work because of collinear Xs
      # so it's wrapped in a try command
      Xall <- cbind(X1.scaled,X0.scaled)
      Xall <- cbind(rep(1,ncol(Xall)),t(Xall))
      Zall <- cbind(Z1,Z0)
      Beta <- try(solve(t(Xall)%*%Xall)%*%t(Xall)%*%t(Zall),silent=TRUE)
      
      # if inverses did not work, we
      # stick with first results    
      if(inherits(Beta,"try-error")) 
       {
        rgV.optim <- rgV.optim.1
       } else {
      # otherwise we run a second optimization with regression starting values
        Beta <- Beta[-1,]
        V    <- Beta%*%t(Beta)
        SV2  <- diag(V)
        SV2 <- SV2 / sum(SV2)
  
      rgV.optim.2 <- optimx(par=SV2, fn=fn.V,
                             gr=NULL, hess=NULL, 
                             method=optimxmethod, itnmax=NULL, hessian=FALSE,
                             control=list(kkt=FALSE,
                                          starttests=FALSE,
                                          dowarn=FALSE,
                                          all.methods=all.methods),
                             X0.scaled = X0.scaled,
                             X1.scaled = X1.scaled,
                             Z0 = Z0,
                             Z1 = Z1,
                             quadopt = quadopt,
                             margin.ipop = Margin.ipop,
                             sigf.ipop = Sigf.ipop,
                             bound.ipop = Bound.ipop
                            )
      # get minimum
      if(verbose==TRUE){print(rgV.optim.2)}
      rgV.optim.2 <- collect.optimx(rgV.optim.2,"min")
  
      # ouput
      if(verbose == TRUE){
      cat("\n Equal weight loss is:",rgV.optim.1$value,"\n")
      cat("\n Regression Loss is:",rgV.optim.2$value,"\n")
      }       
      # and keep the better optim results    
      if(rgV.optim.1$value < rgV.optim.2$value) 
       {
        rgV.optim <- rgV.optim.1
       } else {
        rgV.optim <- rgV.optim.2
       }
      } # close if statement for second regression based optimization attempt
     
      # final V weights from optimization
      solution.v   <- abs(rgV.optim$par)/sum(abs(rgV.optim$par))
     } else { # jump here if only optimize over W


     cat("\n****************",
         "\n optimization over w weights: computing synthtic control unit \n","\n\n"
         )

     if(nrow(X0)==1)
      {
       custom.v <- 1 # this is the case where only one predictor is specified V is the identity matrix
      } else {
     # if user choose to supply v manually:
     cat("\n****************",
         "\n v weights supplied manually: computing synthtic control unit \n","\n\n"
         )
     if(length(custom.v) != nvarsV){stop("custom.V misspecified: length(custom.V) != nrow(X1)")}
     if(mode(custom.v)!="numeric"){stop("custom.V must be numeric")}
      }

    # enter solution.V
    rgV.optim  <- NULL
    solution.v <- abs(custom.v)/sum(custom.v)
    
  } # close else statment for by-passing V optimization


    # last step: now recover solution.w
    V <- diag(x=as.numeric(solution.v),nrow=nvarsV,ncol=nvarsV)
    H <- t(X0.scaled) %*% V %*% (X0.scaled)
    a <- X1.scaled
    c <- -1*c(t(a) %*% V %*% (X0.scaled) )
    A <- t(rep(1, length(c)))
    b <- 1
    l <- rep(0, length(c))
    u <- rep(1, length(c))
    r <- 0

     if(quadopt=="ipop"){
    res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r,
                margin = Margin.ipop, maxiter = 1000, sigf = Sigf.ipop, bound = Bound.ipop)
    solution.w <- as.matrix(primal(res))
    } else {
    # LowRankQP
     if(quadopt=="LowRankQP"){
      res <- LowRankQP(Vmat=H,dvec=c,Amat=A,bvec=1,uvec=rep(1,length(c)),method="LU")
      solution.w <- as.matrix(res$alpha)
     } 
    }

    rownames(solution.w) <- colnames(X0)
    colnames(solution.w) <- "w.weight"
    names(solution.v) <- rownames(X0)

    loss.w <- t(X1.scaled - X0.scaled %*% solution.w) %*%
      V %*% (X1.scaled - X0.scaled %*% solution.w)

    loss.v <-
      t(Z1 - Z0 %*% as.matrix(solution.w)) %*%
        (Z1 - Z0 %*% as.matrix(solution.w)) 
    loss.v <- loss.v/nrow(Z0)      
 
    # produce viewable output
    cat("\n****************",
        "\n****************",
        "\n****************",
        "\n\nMSPE (LOSS V):", loss.v,
  #      "\n\nLOSS (W):", loss.w,
        "\n\nsolution.v:\n", round(as.numeric(solution.v), 10),
        "\n\nsolution.w:\n", round(as.numeric(solution.w), 10),
        "\n\n"
        )
        
    optimize.out <- list(
                         solution.v = solution.v,
                         solution.w = solution.w,
                         loss.v = loss.v,
                         loss.w = loss.w,
                         custom.v = custom.v,
                         rgV.optim = rgV.optim
                         )

    return(invisible(optimize.out))

  }

