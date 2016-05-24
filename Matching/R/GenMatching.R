FastMatchC <- function(N, xvars, All, M, cdd, ww, Tr, Xmod, weights)
  {
    ret <- .Call("FastMatchC", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                 as.double(cdd), as.double(ww), as.double(Tr),
                 as.double(Xmod), as.double(weights),
                 PACKAGE="Matching")
    return(ret)
  }

MatchGenoudStage1 <- function(Tr=Tr, X=X, All=All, M=M, weights=weights,
                              tolerance)
  {
    N  <- nrow(X)
    xvars <- ncol(X)
    
# if SATC is to be estimated the treatment indicator is reversed    
    if (All==2)
      Tr <- 1-Tr

# check on the number of matches, to make sure the number is within the limits
# feasible given the number of observations in both groups.
    if (All==1)
      {
        M <- min(M,min(sum(Tr),sum(1-Tr)));        
      } else {
        M <- min(M,sum(1-Tr));
      }

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
    Mu.X  <- matrix(0, xvars, 1)
    Sig.X <- matrix(0, xvars, 1)

    weights.sum <- sum(weights)

    for (k in 1:xvars)
      {
        Mu.X[k,1] <- sum(X[,k]*weights)/weights.sum;
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weights)/weights.sum-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        if(Sig.X[k,1] < tolerance)
          Sig.X[k,1] <- tolerance
        X[,k]=eps/Sig.X[k,1]
      } #end of k loop

    ret <- list(Tr=Tr, X=X, All=All, M=M, N=N)
    return(ret)
  } #end of MatchGenoudStage1


###############################################################################
## For Caliper!
##
###############################################################################

MatchGenoudStage1caliper <- function(Tr=Tr, X=X, All=All, M=M, weights=weights,
                                     exact=exact, caliper=caliper,
                                     distance.tolerance, tolerance)
  {
    N  <- nrow(X)
    xvars <- ncol(X)
    weights.orig  <- as.matrix(weights)

    if (!is.null(exact))
      {
        exact = as.vector(exact)
        nexacts = length(exact)
        if ( (nexacts > 1) & (nexacts != xvars) )
          {
            warning("length of exact != ncol(X). Ignoring exact option")
            exact <- NULL
          } else if (nexacts==1 & (xvars > 1) ){
            exact <- rep(exact, xvars)
          }
      }

    if (!is.null(caliper))
      {
        caliper = as.vector(caliper)
        ncalipers = length(caliper)
        if ( (ncalipers > 1) & (ncalipers != xvars) )
          {
            warning("length of caliper != ncol(X). Ignoring caliper option")
            caliper <- NULL
          } else if (ncalipers==1 & (xvars > 1) ){
            caliper <- rep(caliper, xvars)
          }
      }

    if (!is.null(caliper))
      {
        ecaliper <- vector(mode="numeric", length=xvars)
        sweights  <- sum(weights.orig)
        for (i in 1:xvars)
          {
            meanX  <- sum( X[,i]*weights.orig )/sweights
            sdX  <- sqrt(sum( (X[,i]-meanX)^2 )/sweights)
            ecaliper[i]  <- caliper[i]*sdX
          }
      } else {
        ecaliper <- NULL
      }

    if (!is.null(exact))
      {
        if(is.null(caliper))
          {
            max.diff <- abs(max(X)-min(X) + distance.tolerance * 100)
            ecaliper <- matrix(max.diff, nrow=xvars, ncol=1)
          }
        
        for (i in 1:xvars)
          {
            if (exact[i])
              ecaliper[i] <- distance.tolerance;
          }
      }        
    
# if SATC is to be estimated the treatment indicator is reversed    
    if (All==2)
      Tr <- 1-Tr

# check on the number of matches, to make sure the number is within the limits
# feasible given the number of observations in both groups.
    if (All==1)
      {
        M <- min(M,min(sum(Tr),sum(1-Tr)));        
      } else {
        M <- min(M,sum(1-Tr));
      }

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
    Mu.X  <- matrix(0, xvars, 1)
    Sig.X <- matrix(0, xvars, 1)

    weights.sum <- sum(weights)

    for (k in 1:xvars)
      {
        Mu.X[k,1] <- sum(X[,k]*weights)/weights.sum;
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weights)/weights.sum-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        if(Sig.X[k,1] < tolerance)
          Sig.X[k,1] <- tolerance        
        X[,k]=eps/Sig.X[k,1]
      } #end of k loop

    ret <- list(Tr=Tr, X=X, All=All, M=M, N=N, ecaliper=ecaliper)
    return(ret)
  } #end of MatchGenoudStage1caliper


###############################################################################
## GenMatch
##
###############################################################################

GenMatch <- function(Tr, X, BalanceMatrix=X, estimand="ATT", M=1,
                     weights=NULL,
                     pop.size = 100, max.generations=100,
                     wait.generations=4, hard.generation.limit=FALSE,
                     starting.values=rep(1,ncol(X)),
                     fit.func="pvals",
                     MemoryMatrix=TRUE,
                     exact=NULL, caliper=NULL, replace=TRUE, ties=TRUE,
                     CommonSupport=FALSE,nboots=0, ks=TRUE, verbose=FALSE,
                     distance.tolerance=0.00001,
                     tolerance=sqrt(.Machine$double.eps),
                     min.weight=0,
                     max.weight=1000,
                     Domains=NULL,
                     print.level=2,
                     project.path=NULL,
                     paired=TRUE,
                     loss=1,
                     data.type.integer=FALSE,
                     restrict=NULL,
                     cluster=FALSE,
                     balance=TRUE, ...)
  {

    requireNamespace("rgenoud")

    Tr <- as.double(Tr)
    X  <- as.matrix(X)
    BalanceMatrix  <- as.matrix(BalanceMatrix)

    if(length(Tr) != nrow(X))
      {
        stop("length(Tr) != nrow(X)")
      }
    if(!is.function(fit.func))
      {
        if(nrow(BalanceMatrix) != length(Tr))
          {
            stop("nrow(BalanceMatrix) != length(Tr)")
          }
      }

    if (is.null(weights))
      {
        weights <- rep(1,length(Tr))
        weights.flag <- FALSE
      } else {
        weights.flag <- TRUE
        weights <- as.double(weights)
        if( length(Tr) != length(weights))
          {
            stop("length(Tr) != length(weights)")
          }            
      }

    isna  <- sum(is.na(Tr)) + sum(is.na(X)) + sum(is.na(weights)) + sum(is.na(BalanceMatrix))
    if (isna!=0)
      {
        stop("GenMatch(): input includes NAs")
        return(invisible(NULL))
      }    

    #check inputs
    if (sum(Tr !=1 & Tr !=0) > 0) {
      stop("Treatment indicator must be a logical variable---i.e., TRUE (1) or FALSE (0)")
    }
    if (var(Tr)==0) {
      stop("Treatment indicator ('Tr') must contain both treatment and control observations")
    }    
    if (distance.tolerance < 0)
      {
        warning("User set 'distance.tolerance' to less than 0.  Resetting to the default which is 0.00001.")
        distance.tolerance <- 0.00001
      }
    #CommonSupport
    if (CommonSupport !=1 & CommonSupport !=0) {
      stop("'CommonSupport' must be a logical variable---i.e., TRUE (1) or FALSE (0)")
    }    
    if(CommonSupport==TRUE)
      {
        tr.min <- min(X[Tr==1,1])
        tr.max <- max(X[Tr==1,1])

        co.min <- min(X[Tr==0,1])
        co.max <- max(X[Tr==0,1])

        if(tr.min >= co.min)
          {
            indx1 <- X[,1] < (tr.min-distance.tolerance)
          } else {
            indx1 <- X[,1] < (co.min-distance.tolerance)            
          }

        if(co.max <= tr.max)
          {        
            indx2 <- X[,1] > (co.max+distance.tolerance)
          } else {
            indx2 <- X[,1] > (tr.max+distance.tolerance)            
          }

        indx3 <- indx1==0 & indx2==0
        Tr <- as.double(Tr[indx3])
        X  <- as.matrix(X[indx3,])
        BalanceMatrix <- as.matrix(BalanceMatrix[indx3,])
        weights <- as.double(weights[indx3])
      }#end of CommonSupport
    
    if (pop.size < 0 | pop.size!=round(pop.size) )
      {
        warning("User set 'pop.size' to an illegal value.  Resetting to the default which is 100.")
        pop.size <- 100
      }
    if (max.generations < 0 | max.generations!=round(max.generations) )
      {
        warning("User set 'max.generations' to an illegal value.  Resetting to the default which is 100.")
        max.generations <-100
      }
    if (wait.generations < 0 | wait.generations!=round(wait.generations) )
      {
        warning("User set 'wait.generations' to an illegal value.  Resetting to the default which is 4.")
        wait.generations <- 4
      }
    if (hard.generation.limit != 0 & hard.generation.limit !=1 )
      {
        warning("User set 'hard.generation.limit' to an illegal value.  Resetting to the default which is FALSE.")
        hard.generation.limit <- FALSE
      }
    if (data.type.integer != 0 & data.type.integer !=1 )
      {
        warning("User set 'data.type.integer' to an illegal value.  Resetting to the default which is TRUE.")
        data.type.integer <- TRUE
      }
    if (MemoryMatrix != 0 & MemoryMatrix !=1 )
      {
        warning("User set 'MemoryMatrix' to an illegal value.  Resetting to the default which is TRUE.")
        MemoryMatrix <- TRUE
      }                
    if (nboots < 0 | nboots!=round(nboots) )
      {
        warning("User set 'nboots' to an illegal value.  Resetting to the default which is 0.")
        nboots <- 0
      }
    if (ks != 0 & ks !=1 )
      {
        warning("User set 'ks' to an illegal value.  Resetting to the default which is TRUE.")
        ks <- TRUE
      }
    if (verbose != 0 & verbose !=1 )
      {
        warning("User set 'verbose' to an illegal value.  Resetting to the default which is FALSE.")
        verbose <- FALSE
      }
    if (min.weight < 0)
      {
        warning("User set 'min.weight' to an illegal value.  Resetting to the default which is 0.")
        min.weight <- 0
      }
    if (max.weight < 0)
      {
        warning("User set 'max.weight' to an illegal value.  Resetting to the default which is 1000.")
        max.weight <- 1000
      }
    if (print.level != 0 & print.level !=1 & print.level !=2 & print.level !=3)
      {
        warning("User set 'print.level' to an illegal value.  Resetting to the default which is 2.")
        print.level <- 2
      }    
    if (paired != 0 & paired !=1 )
      {
        warning("User set 'paired' to an illegal value.  Resetting to the default which is TRUE.")
        paired <- FALSE
      }    
    ##from Match()
    if (tolerance < 0)
      {
        warning("User set 'tolerance' to less than 0.  Resetting to the default which is 0.00001.")
        tolerance <- 0.00001
      }
    if (M < 1)
      {
        warning("User set 'M' to less than 1.  Resetting to the default which is 1.")
        M <- 1
      }
    if ( M!=round(M) )
      {
        warning("User set 'M' to an illegal value.  Resetting to the default which is 1.")
        M <- 1        
      }
    if (replace!=FALSE & replace!=TRUE)
      {
        warning("'replace' must be TRUE or FALSE.  Setting to TRUE")
        replace <- TRUE
      }
    if(replace==FALSE)
      ties <- FALSE
    if (ties!=FALSE & ties!=TRUE)
      {
        warning("'ties' must be TRUE or FALSE.  Setting to TRUE")
        ties <- TRUE
      }

    #print warning if pop.size, max.generations and wait.generations are all set to their original values
    if(pop.size==100 & max.generations==100 & wait.generations==4)
      {
        warning("The key tuning parameters for optimization were are all left at their default values.  The 'pop.size' option in particular should probably be increased for optimal results.  For details please see the help page and http://sekhon.berkeley.edu/papers/MatchingJSS.pdf")
      }

    #loss function
    if (is.double(loss))
      {
        if (loss==1)  {
          loss.func=sort
          lexical=ncol(BalanceMatrix)
          if(ks)
            lexical=lexical+lexical
        } else if(loss==2) {
          loss.func=min
          lexical=0
        } else{
          stop("unknown loss function")
        }
      } else if (is.function(loss)) {
        loss.func=loss
        lexical=1
      } else {
        stop("unknown loss function")
      }

    #set lexical for fit.func
    if (is.function(fit.func))
      {
        lexical = 1
      } else if (fit.func=="qqmean.max" | fit.func=="qqmedian.max" | fit.func=="qqmax.max")   {
        lexical=ncol(BalanceMatrix)
      } else if (fit.func!="qqmean.mean" & fit.func!="qqmean.max" &
                 fit.func!="qqmedian.median" & fit.func!="qqmedian.max"
                 & fit.func!="pvals") {
        stop("invalid 'fit.func' argument")
      } else  if (!fit.func=="pvals") {
        lexical = 0
      }

    if(replace==FALSE)
      {
        #replace==FALE, needs enough observation
        #ATT
        orig.weighted.control.nobs <- sum(weights[Tr!=1])
        orig.weighted.treated.nobs <- sum(weights[Tr==1])
        if(estimand=="ATC")
          {
            if (orig.weighted.treated.nobs < orig.weighted.control.nobs)
              {
                warning("replace==FALSE, but there are more (weighted) control obs than treated obs.  Some obs will be dropped.  You may want to estimate ATC instead")
              }
          } else if(estimand=="ATE") 
            {
              #ATE
              if (orig.weighted.treated.nobs > orig.weighted.control.nobs)
                {
                  warning("replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.")
                }              
              if (orig.weighted.treated.nobs < orig.weighted.control.nobs)
                {
                  warning("replace==FALSE, but there are more (weighted) control obs than treated obs.  Some control obs will not be matched.  You may want to estimate ATT instead.")
                }
            } else {
              #ATT
              if (orig.weighted.treated.nobs > orig.weighted.control.nobs)
                {
                  warning("replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.")
                }
            }
        
        #we need a restrict matrix if we are going to not do replacement
        if(is.null(restrict))
          {
            restrict <- t(as.matrix(c(0,0,0)))
          }
      }#end of replace==FALSE    

    

    #check the restrict matrix input
    if(!is.null(restrict))
      {
        if(!is.matrix(restrict))
          stop("'restrict' must be a matrix of restricted observations rows and three columns: c(i,j restriction)")

        if(ncol(restrict)!=3 )
          stop("'restrict' must be a matrix of restricted observations rows and three columns: c(i,j restriction)")

        restrict.trigger <- TRUE
      }  else {
        restrict.trigger <- FALSE
      }

    if(!is.null(caliper) | !is.null(exact) | restrict.trigger | !ties)
      {
        GenMatchCaliper.trigger <- TRUE
      } else {
        GenMatchCaliper.trigger <- FALSE
      }

    isunix  <- .Platform$OS.type=="unix"
    if (is.null(project.path))
      {
        if (print.level < 3 & isunix)
          {
            project.path="/dev/null"
          } else {
            project.path=paste(tempdir(),"/genoud.pro",sep="")
            
            #work around for rgenoud bug
            #if (print.level==3)
            #print.level <- 2
          }
      } 
    
    nvars <- ncol(X)
    balancevars <- ncol(BalanceMatrix)

    if (is.null(Domains))
      {
        Domains <- matrix(min.weight, nrow=nvars, ncol=2)
        Domains[,2] <- max.weight
      } else {
        indx <- (starting.values < Domains[,1]) | (starting.values > Domains[,2])
        starting.values[indx] <- round( (Domains[indx,1]+Domains[indx,2])/2 )
      }

    # create All
    if (estimand=="ATT")
      {
        All  <- 0
      } else if(estimand=="ATE") {
        All  <- 1
      } else if(estimand=="ATC") {
        All  <- 2
      } else {
        All  <- 0
        warning("User set 'estimand' to an illegal value.  Resetting to the default which is 'ATT'")
      }

    #stage 1 Match, only needs to be called once    
    if(!GenMatchCaliper.trigger)
      {

        s1 <- MatchGenoudStage1(Tr=Tr, X=X, All=All, M=M, weights=weights,
                                tolerance=tolerance);
        s1.Tr <- s1$Tr
        s1.X <- s1$X
        s1.All <- s1$All
        s1.M <- s1$M
        s1.N <- s1$N
        rm(s1)
      } else {
        s1 <- MatchGenoudStage1caliper(Tr=Tr, X=X, All=All, M=M, weights=weights,
                                       exact=exact, caliper=caliper,
                                       distance.tolerance=distance.tolerance,
                                       tolerance=tolerance)
        s1.Tr <- s1$Tr
        s1.X <- s1$X
        s1.All <- s1$All
        s1.M <- s1$M
        s1.N <- s1$N
        s1.ecaliper  <- s1$ecaliper
        
        if (is.null(s1.ecaliper))
          {
            caliperFlag  <- 0
            Xorig  <- 0
            CaliperVec  <- 0
          } else {
            caliperFlag  <- 1
            Xorig  <- X
            CaliperVec  <- s1$ecaliper
          }
        rm(s1)        
      } #GenMatchCaliper.trigger

    genoudfunc  <- function(x)
      {
        wmatrix <- diag(x, nrow=nvars)
        if ( min(eigen(wmatrix, symmetric=TRUE, only.values=TRUE)$values) < tolerance )
            wmatrix <- wmatrix + diag(nvars)*tolerance
        
        ww <- chol(wmatrix)

        if(!GenMatchCaliper.trigger)
          {
            
            if (weights.flag==TRUE)
              {
                FastMatchC.internal <- function(N, xvars, All, M, cdd, ww, Tr, Xmod, weights)
                  {
                    ret <- .Call("FastMatchC", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                                 as.double(cdd), as.double(ww), as.double(Tr),
                                 as.double(Xmod), as.double(weights),
                                 PACKAGE="Matching")
                    return(ret)
                  }
                
                rr <- FastMatchC.internal(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                                          cdd=distance.tolerance, ww=ww, Tr=s1.Tr, Xmod=s1.X,
                                          weights=weights)
              } else {
                FasterMatchC.internal <- function(N, xvars, All, M, cdd, ww, Tr, Xmod, weights)
                  {
                    ret <- .Call("FasterMatchC", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                                 as.double(cdd), as.double(ww), as.double(Tr),
                                 as.double(Xmod), 
                                 PACKAGE="Matching")
                    return(ret)
                  }
                
                rr <- FasterMatchC.internal(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                                             cdd=distance.tolerance, ww=ww, Tr=s1.Tr, Xmod=s1.X)
              } #end of weights.flag
          } else {
            if (weights.flag==TRUE)
              {
                MatchLoopC.internal <- function(N, xvars, All, M, cdd, caliperflag, replace, ties, ww, Tr, Xmod, weights, CaliperVec,
                                                Xorig, restrict.trigger, restrict)
                  {
                    
                    if(restrict.trigger)
                      {
                        restrict.nrow <- nrow(restrict)
                      } else {
                        restrict.nrow <- 0
                      }    
                    
                    ret <- .Call("MatchLoopC", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                                 as.double(cdd), as.integer(caliperflag), as.integer(replace), as.integer(ties), as.double(ww), as.double(Tr),
                                 as.double(Xmod), as.double(weights), as.double(CaliperVec), as.double(Xorig),
                                 as.integer(restrict.trigger), as.integer(restrict.nrow), as.double(restrict),
                                 #next line is sets the DiagWeightMatrixFlag
                                 as.double(1),                             
                                 PACKAGE="Matching")
                    return(ret)
                  } #end of MatchLoopC.internal
                
                rr <- MatchLoopC.internal(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                                          cdd=distance.tolerance,
                                          caliperflag=caliperFlag,
                                          replace=replace, ties=ties,
                                          ww=ww, Tr=s1.Tr, Xmod=s1.X, weights=weights,
                                          CaliperVec=CaliperVec, Xorig=Xorig,
                                          restrict.trigger=restrict.trigger, restrict=restrict)
              } else {
                MatchLoopCfast.internal <- function(N, xvars, All, M, cdd, caliperflag, replace, ties, ww, Tr, Xmod, CaliperVec, Xorig,
                                                    restrict.trigger, restrict)
                  {
                    
                    if(restrict.trigger)
                      {
                        restrict.nrow <- nrow(restrict)
                      } else {
                        restrict.nrow <- 0
                      }    
                    
                    ret <- .Call("MatchLoopCfast", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                                 as.double(cdd), as.integer(caliperflag), as.integer(replace), as.integer(ties), as.double(ww), as.double(Tr),
                                 as.double(Xmod), as.double(CaliperVec), as.double(Xorig),
                                 as.integer(restrict.trigger), as.integer(restrict.nrow), as.double(restrict),
                                 #next line is the DiagWeightMatrixFlag
                                 as.double(1),                              
                                 PACKAGE="Matching")
                    return(ret)
                  } #end of MatchLoopCfast.internal
                
                
                rr <- MatchLoopCfast.internal(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                                              cdd=distance.tolerance,
                                              caliperflag=caliperFlag,
                                              replace=replace, ties=ties,
                                              ww=ww, Tr=s1.Tr, Xmod=s1.X, 
                                              CaliperVec=CaliperVec, Xorig=Xorig,
                                              restrict.trigger=restrict.trigger, restrict=restrict)
              } #end of weights.flag
            
            #no matches
            if(rr[1,1]==0) {
              warning("no valid matches found in GenMatch evaluation") 
              return(rep(-9999, balancevars*2))
            }
            
            rr <- rr[,c(4,5,3)]            
          } #Caliper.Trigger

        #should be the same as GenBalance() in GenBalance.R but we need to include it here because of
        #cluster scoping issues.
        GenBalance.internal <-
          function(rr, X, nvars=ncol(X), nboots = 0, ks=TRUE, verbose = FALSE, paired=TRUE)
          {

            #CUT-AND-PASTE from GenBalance.R, the functions before GenBalance. but get rid of warn *switch*
            MATCHpt <- function(q, df, ...)
              {
                #don't know how general it is so let's try to work around it.
                ret=pt(q,df, ...)
                
                if (is.na(ret)) {
                  ret   <- pt(q, df, ...)
                  if(is.na(ret))
                    warning("pt() generated NaN. q:",q," df:",df,"\n",date())
                }
                
                return(ret)
              } #end of MATCHpt

            Mt.test.pvalue  <- function(Tr, Co, weights)
              {
                v1  <- Tr-Co
                estimate  <- sum(v1*weights)/sum(weights)
                var1  <- sum( ((v1-estimate)^2)*weights )/( sum(weights)*sum(weights) )
                
                if (estimate==0 & var1==0)
                  {
                    return(1)
                  }
                
                statistic  <- estimate/sqrt(var1)
                #    p.value    <- (1-pnorm(abs(statistic)))*2
                p.value    <- (1-MATCHpt(abs(statistic), df=sum(weights)-1))*2
                
                return(p.value)
              } #end of Mt.test.pvalue

            Mt.test.unpaired.pvalue  <- function(Tr, Co, weights)
              {
                obs <- sum(weights)
                
                mean.Tr <- sum(Tr*weights)/obs
                mean.Co <- sum(Co*weights)/obs
                estimate <- mean.Tr-mean.Co
                var.Tr  <- sum( ( (Tr - mean.Tr)^2 )*weights)/(obs-1)
                var.Co  <- sum( ( (Co - mean.Co)^2 )*weights)/(obs-1)
                dim <- sqrt(var.Tr/obs + var.Co/obs)
                
                if (estimate==0 & dim==0)
                  {
                    return(1)
                  }
                
                statistic  <- estimate/dim
                
                a1 <- var.Tr/obs
                a2 <- var.Co/obs
                dof <- ((a1 + a2)^2)/( (a1^2)/(obs - 1) + (a2^2)/(obs - 1) )    
                p.value    <- (1-MATCHpt(abs(statistic), df=dof))*2    
                
                return(p.value)
              } #end of Mt.test.unpaired.pvalue

            ks.fast <- function(x, y, n.x, n.y, n)
              {
                w <- c(x, y)
                z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
                z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]        
                
                return( max(abs(z)) )
              } #ks.fast

            index.treated <- rr[,1]
            index.control <- rr[,2]
            weights <- rr[,3]
            
            tol  <- .Machine$double.eps*100  
            storage.t <- c(rep(9,nvars))
            storage.k <- c(rep(9,nvars))
            fs.ks     <- matrix(nrow=nvars, ncol=1)
            s.ks      <- matrix(nrow=nvars, ncol=1)  
            bbcount   <- matrix(0, nrow=nvars, ncol=1)
            dummy.indx  <- matrix(0, nrow=nvars, ncol=1)
            
            w  <- c(X[,1][index.treated], X[,1][index.control])
            obs <- length(w)
            n.x  <- length(X[,1][index.treated])
            n.y  <- length(X[,1][index.control])
            cutp <- round(obs/2)  
            w  <- matrix(nrow=obs, ncol=nvars)
            
            for (i in 1:nvars)
              {
                w[,i] <- c(X[,i][index.treated], X[,i][index.control])
                
                if(paired)
                  {
                    t.out <- Mt.test.pvalue(X[,i][index.treated],
                                            X[,i][index.control],
                                            weights = weights)
                  } else {
                    t.out <- Mt.test.unpaired.pvalue(X[,i][index.treated],
                                                     X[,i][index.control],
                                                     weights = weights)
                  }
                
                storage.t[i] <- t.out            
        
                dummy.indx[i]  <- length(unique(X[,i])) < 3
                
                if (!dummy.indx[i] & ks & nboots > 9)
                  {
                    fs.ks[i]  <- ks.fast(X[,i][index.treated], X[,i][index.control],
                                         n.x=n.x, n.y=n.y, n=obs)
                  } else if(!dummy.indx[i] & ks)
                    {
                      
                      storage.k[i] <- Mks.test(X[,i][index.treated], X[,i][index.control])$p.value

                    }
              }#end of i loop

            
            if (ks & nboots > 9)
              {
                n.x  <- cutp
                n.y  <- obs-cutp
                for (b in 1:nboots)
                  {
                    sindx  <- sample(1:obs, obs, replace = TRUE)
                    
                    for (i in 1:nvars)
                      {
                        
                        if (dummy.indx[i])
                          next;
                        
                        X1tmp <- w[sindx[1:cutp],i ]
                        X2tmp <- w[sindx[(cutp + 1):obs], i]
                        s.ks[i] <- ks.fast(X1tmp, X2tmp, n.x=n.x, n.y=n.y, n=obs)
                        if (s.ks[i] >= (fs.ks[i] - tol) )
                          bbcount[i]  <-  bbcount[i] + 1
                      }#end of i loop
                  } #end of b loop
                
                for (i in 1:nvars)
                  {
                    
                    if (dummy.indx[i])
                      {
                        storage.k[i]  <- 9
                        next;
                      }
                    
                    storage.k[i]  <- bbcount[i]/nboots
                    
                  }
                storage.k[storage.k==9]=storage.t[storage.k==9]
                output <- c(storage.t, storage.k)
              } else if(ks){
                storage.k[storage.k==9]=storage.t[storage.k==9]                
                output <- c(storage.t, storage.k)
              } else {
                output <- storage.t
              }
            
            if(sum(is.na(output)) > 0) {
              output[is.na(output)] = 2
              warning("output has NaNs")
            }
            
            if (verbose == TRUE)
              {
                cat("\n")
                for (i in 1:nvars)
                  {
                    cat("\n", i, " t-test p-val  =", storage.t[i], "\n" )
                    if(ks)
                      cat(" ", i, "  ks-test p-val = ", storage.k[i], " \n",sep="")
                  }
                cat("\nsorted return vector:\n", sort(output), "\n")
                cat("number of return values:", length(output), "\n")
              }
            
            return(output)
          } #end of GenBalance.internal

        GenBalanceQQ.internal <- function(rr, X, summarystat="mean", summaryfunc="mean")
          {
            index.treated <- rr[,1]
            index.control <- rr[,2]
    
            nvars <- ncol(X)
            qqsummary   <- c(rep(NA,nvars))
            
            for (i in 1:nvars)
              {    
                
                qqfoo <- qqstats(X[,i][index.treated], X[,i][index.control], standardize=TRUE)
                
                if (summarystat=="median")
                  {
                    qqsummary[i] <- qqfoo$mediandiff
                  } else if (summarystat=="max")  {
                    qqsummary[i] <- qqfoo$maxdiff    
                  } else {
                    qqsummary[i] <- qqfoo$meandiff
                  }
                
              } #end of for loop
            
            
            if (summaryfunc=="median")
              {
                return(median(qqsummary))
              } else if (summaryfunc=="max")  {
                return(sort(qqsummary, decreasing=TRUE))
              } else if (summaryfunc=="sort")  {
                return(sort(qqsummary, decreasing=TRUE))
              } else {
                return(mean(qqsummary))
              }    
          } #end of GenBalanceQQ.internal


        if (is.function(fit.func)) {
          a <- fit.func(rr, BalanceMatrix)
          return(a)
        } else if (fit.func=="pvals")
          {
            a <- GenBalance.internal(rr=rr, X=BalanceMatrix, nvars=balancevars, nboots=nboots,
                                     ks=ks, verbose=verbose, paired=paired)
            a <- loss.func(a)
            return(a)
          } else if (fit.func=="qqmean.mean") {
            a <- GenBalanceQQ.internal(rr=rr, X=BalanceMatrix, summarystat="mean", summaryfunc="mean")
            return(a)
          } else if (fit.func=="qqmean.max") {
            a <- GenBalanceQQ.internal(rr=rr, X=BalanceMatrix, summarystat="mean", summaryfunc="max")
            return(a)
          } else if (fit.func=="qqmax.mean") {
            a <- GenBalanceQQ.internal(rr=rr, X=BalanceMatrix, summarystat="max", summaryfunc="mean")
            return(a)            
          } else if (fit.func=="qqmax.max") {
            a <- GenBalanceQQ.internal(rr=rr, X=BalanceMatrix, summarystat="max", summaryfunc="max")
            return(a)                        
          } else if (fit.func=="qqmedian.median") {
            a <- GenBalanceQQ.internal(rr=rr, X=BalanceMatrix, summarystat="median", summaryfunc="median")
            return(a)
          } else if (fit.func=="qqmedian.max") {
            a <- GenBalanceQQ.internal(rr=rr, X=BalanceMatrix, summarystat="median", summaryfunc="max")
            return(a)
          } 
      } #end genoudfunc

    #cluster info
    clustertrigger=1
    if (is.logical(cluster))
      {
        if (cluster==FALSE)  {
          clustertrigger=0
        } else {
          stop("cluster option must be either FALSE, an object of the 'cluster' class (from the 'parallel' package) or a list of machines so 'genoud' can create such an object")
        }
      }
    
    if(clustertrigger) {
      parallel.exists = requireNamespace("parallel")
      if (!parallel.exists) {
        stop("The 'cluster' feature cannot be used unless the package 'parallel' can be loaded.")
      }
    } 

    if(clustertrigger)
      {

        GENclusterExport <- function (cl, list, envir = .GlobalEnv) 
          {
            gets <- function(n, v) {
              assign(n, v, envir = envir)
              NULL
            }
            for (name in list) {
              parallel::clusterCall(cl, gets, name, get(name))
            }
          }        

        if (class(cluster)[1]=="SOCKcluster" | class(cluster)[1]=="PVMcluster" | class(cluster)[1]=="spawnedMPIcluster" | class(cluster)[1]=="MPIcluster") {
          clustertrigger=1
          cl <- cluster
          cl.genoud <- cl
        } else {
          clustertrigger=2
          cluster <- as.vector(cluster)
          cat("Initializing Cluster\n")
          cl <- parallel::makePSOCKcluster(cluster)
          cl.genoud <- cl
        }      
      } else {
        cl.genoud <- FALSE
      }#end of clustertrigger

    if (clustertrigger > 0)
      {
        #create restrict.summary, because passing the entire restrict matrix is too much
        
        parallel::clusterEvalQ(cl, library("Matching"))
        GENclusterExport(cl, c("s1.N", "s1.All", "s1.M", "s1.Tr", "s1.X", "nvars",
                               "tolerance", "distance.tolerance", "weights",
                               "BalanceMatrix", "balancevars", "nboots", "ks", "verbose", "paired", "loss.func",
                               "fit.func"))

        if(GenMatchCaliper.trigger) {
          GENclusterExport(cl, c("caliperFlag", "CaliperVec", "Xorig", "restrict.trigger", "restrict","replace"))
        }
        
        GENclusterExport(cl, "genoudfunc")
      }

    do.max <- FALSE    
    if(!is.function(fit.func))
      {
        if (fit.func=="pvals") 
          do.max <- TRUE
      }

    rr <- rgenoud::genoud(genoudfunc, nvars=nvars, starting.values=starting.values,
                 pop.size=pop.size, max.generations=max.generations,
                 wait.generations=wait.generations, hard.generation.limit=hard.generation.limit,
                 Domains=Domains,
                 MemoryMatrix=MemoryMatrix,
                 max=do.max, gradient.check=FALSE, data.type.int=data.type.integer,
                 hessian=FALSE,
                 BFGS=FALSE, project.path=project.path, print.level=print.level,
                 lexical=lexical,
                 cluster=cl.genoud,
                 balance=balance,
                 ...)
    
    wmatrix <- diag(rr$par, nrow=nvars)
    
    if ( min(eigen(wmatrix, symmetric=TRUE, only.values=TRUE)$values) < tolerance )
      wmatrix <- wmatrix + diag(nvars)*tolerance
        
    ww <- chol(wmatrix)

    if(!GenMatchCaliper.trigger)
      {
        mout <- FastMatchC(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                           cdd=distance.tolerance, ww=ww, Tr=s1.Tr, Xmod=s1.X,
                           weights=weights)
        rr2 <- list(value=rr$value, par=rr$par, Weight.matrix=wmatrix, matches=mout, ecaliper=NULL)
      } else {
        mout <- MatchLoopC(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                           cdd=distance.tolerance,
                           caliperflag=caliperFlag,
                           replace=replace, ties=ties,
                           ww=ww, Tr=s1.Tr, Xmod=s1.X, weights=weights,
                           CaliperVec=CaliperVec, Xorig=Xorig,
                           restrict.trigger=restrict.trigger, restrict=restrict,
                           DiagWeightMatrixFlag=1)

        #no matches
        if(mout[1,1]==0) {
          warning("no valid matches found by GenMatch") 
        }        
        
        rr2 <- list(value=rr$value, par=rr$par, Weight.matrix=wmatrix, matches=mout, ecaliper=CaliperVec)
      }

    if (clustertrigger==2)
      parallel::stopCluster(cl)    
    
    class(rr2) <- "GenMatch"
    return(rr2)
  } #end of GenMatch
