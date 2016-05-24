# Jasjeet S. Sekhon <sekhon@berkeley.edu>
# HTTP://sekhon.berkeley.edu/
# UC Berkeley

# Match(): function to estimate treatments using a matching estimator.
# Currently only the ability to estimate average treatment effects
# using the approach of Abadie and Imbens is implemented.  In the
# future, quantile treatment effects will be implemented along with
# the ability to use robust estimation when estimating the propensity
# score. MatchBalance(), and balanceUV() test for balance.

Match  <- function(Y=NULL,Tr,X,Z=X,V=rep(1,length(Y)), estimand="ATT", M=1,
                   BiasAdjust=FALSE,exact=NULL,caliper=NULL, replace=TRUE, ties=TRUE,
                   CommonSupport=FALSE,Weight=1,Weight.matrix=NULL, weights=NULL,
                   Var.calc=0, sample=FALSE, restrict=NULL, match.out=NULL,
                   distance.tolerance=0.00001, tolerance=sqrt(.Machine$double.eps),
                   version="standard")
  {

    BiasAdj  <- as.double(BiasAdjust)
    sample  <- as.double(sample)

    if ( (BiasAdj != 0) & (BiasAdj != 1) )
      {
        warning("User set 'BiasAdjust' to a non-logical value.  Resetting to the default which is FALSE.")        
        BiasAdj <- 0
      }
    
    #we don't need to use a Y
    if (is.null(Y))
      {
        Y = rep(0, length(Tr))
        version <- "fast"

        if(BiasAdj)
          {
            warning("'BiasAdjust' set to FALSE because Y is NULL")
            BiasAdj <- FALSE
          }
      }
    
    Y  <- as.double(Y)
    Tr <- as.double(Tr)
    X  <- as.matrix(X)
    Z  <- as.matrix(Z)
    V  <- as.matrix(V)

    orig.nobs  <- length(Y)
    nobs  <- orig.nobs
    xvars <- ncol(X)    
    orig.tr.nobs <- length(Tr)
    if (orig.tr.nobs != orig.nobs)
      {
        stop("length(Y) != length(Tr)")
      }
    if( orig.tr.nobs != nrow(X))
      {
        stop("length(Tr) != nrow(X)")
      }    
    
    if( orig.nobs != nrow(X))
      {
        stop("length(Y) != nrow(X)")
      }
    if( orig.nobs != nrow(V))
      {
        stop("length(Y) != nrow(V)")
      }
    if( orig.nobs != nrow(Z))
      {
        stop("length(Y) != nrow(Z)")
      }

    if (is.null(weights))
      {
        weights <- rep(1,length(Y))
        weights.flag <- FALSE
      } else {
        weights.flag <- TRUE
        weights <- as.double(weights)
        if( orig.tr.nobs != length(weights))
          {
            stop("length(Tr) != length(weights)")
          }        
      }

    isna  <- sum(is.na(Y)) + sum(is.na(Tr)) + sum(is.na(X)) + sum(is.na(weights)) + sum(is.na(Z)) + sum(is.na(V))
    if (isna!=0)
      {
        stop("Match(): input includes NAs")
        return(invisible(NULL))
      }

    if (sum(Tr !=1 & Tr !=0) > 0) {
      stop("Treatment indicator ('Tr') must be a logical variable---i.e., TRUE (1) or FALSE (0)")
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
        Y  <- as.double(Y[indx3])
        Tr <- as.double(Tr[indx3])
        X  <- as.matrix(X[indx3,])
        Z  <- as.matrix(Z[indx3,])
        V  <- as.matrix(V[indx3,])
        weights <- as.double(weights[indx3])

        #let's recalculate these for common support
        orig.nobs  <- length(Y)
        nobs  <- orig.nobs
      }#end of CommonSupport

    #check additional inputs
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
    if ( M != round(M) )
      {
        warning("User set 'M' to an illegal value.  Resetting to the default which is 1.")
        M <- 1        
      }
    if (Var.calc < 0)
      {
        warning("User set 'Var.calc' to less than 0.  Resetting to the default which is 0.")
        Var.calc <- 0
      }
    if ( (sample != 0) & (sample != 1) )
      {
        warning("User set 'sample' to a non-logical value.  Resetting to the default which is FALSE.")        
        sample <- 0
      }
    if (Weight != 1 & Weight != 2 & Weight != 3)
      {
        warning("User set 'Weight' to an illegal value.  Resetting to the default which is 1.")        
        Weight <- 1
      }
    if (version!="fast" & version != "standard" & version != "legacy" & version != "Matchby" & version != "MatchbyAI")
      {
        warning("User set 'version' to an illegal value.  Resetting to the default which is 'standard'.")        
        version <- "standard"
      }

    if(version=="Matchby")
      {
        version <- "fast"
        Matchby.call <- TRUE
        MatchbyAI <- FALSE
      } else if(version=="MatchbyAI")
      {
        version <- "standard"
        Matchby.call <- TRUE
        MatchbyAI <- TRUE
      } else {
        Matchby.call <- FALSE
        MatchbyAI <- FALSE
      }

    if (Var.calc !=0 & version=="fast")
      {
        warning("Var.calc cannot be estimate when version=='fast'")
        Var.calc=0
      }
    if (BiasAdj!=FALSE & version=="fast")
      {
        warning("Bias Adjustment cannot be estimated when version=='fast'")
        BiasAdj=0
      }
    if (replace!=FALSE & replace!=TRUE)
      {
        warning("'replace' must be TRUE or FALSE.  Setting to TRUE")
        replace <- TRUE
      }
    if(replace==FALSE)
      {
        ties <- FALSE        
        version="fast"
        if (version=="legacy")
          warning("'version' is set to 'fast' because replace==FALSE")
      }
    if (ties!=FALSE & ties!=TRUE)
      {
        warning("'ties' must be TRUE or FALSE.  Setting to TRUE")
        ties <- TRUE
      }
    if(ties==FALSE)
      {
        version="fast"
        if (version=="legacy")
          warning("'version' is set to 'fast' because ties==FALSE")

        if(BiasAdjust==TRUE) {
          warning("Bias Adjustment can only be estimated when ties==TRUE and replace=TRUE.  Setting BiasAdjust=FALSE")
          BiasAdjust <- FALSE
          BiasAdj  <- 0
        }
      }

    if (!is.null(match.out) & class(match.out) != "Match") {
      warning("match.out object not of class 'Match'")
      return(invisible(NULL))
    }

    ccc  <- tolerance
    cdd  <- distance.tolerance

    orig.treated.nobs  <- sum(Tr==1)
    orig.control.nobs  <- sum(Tr==0)
    orig.wnobs  <- sum(weights)
    orig.weighted.treated.nobs <- sum( weights[Tr==1] )
    orig.weighted.control.nobs <- sum( weights[Tr==0] )    
    weights.orig  <- as.matrix(weights)
    zvars <- ncol(Z);

    estimand.orig <- estimand
    if (estimand=="ATT")
      {
        estimand  <- 0

        if(BiasAdj==1 & orig.treated.nobs<zvars)
          {
            warning("Fewer treated obs than variables in 'Z': BiasAdjust set to FALSE")
            BiasAdj=0            
          }

      } else if(estimand=="ATE") {
        estimand  <- 1

        if(BiasAdj==1 & orig.nobs<zvars)
          {
            warning("Fewer obs than variables in 'Z': BiasAdjust set to FALSE")
            BiasAdj=0            
          }
      } else if(estimand=="ATC") {
        estimand  <- 2

        if(BiasAdj==1 & orig.control.nobs<zvars)
          {
            warning("Fewer control obs than variables in 'Z': BiasAdjust set to FALSE")
            BiasAdj=0            
          }

      } else {
        estimand  <- 0
        warning("User set 'estimand' to an illegal value.  Resetting to the default which is 'ATT'")
      }

    if (!is.null(Weight.matrix))
      {

        if(class(Weight.matrix)=="GenMatch")
          {
            Weight.matrix = Weight.matrix$Weight.matrix
          }
        
        if (Weight==2)
          {
            warning("User supplied 'Weight.matrix' is being used even though 'Weight' is not set equal to 3")
          }
        Weight  <- 3
      } else {
        Weight.matrix <- dim(X)[2]
      }

    if(Var.calc > orig.weighted.treated.nobs)
      {
        warning("'Var.calc' > the number of treated obs: 'Var.calc' reset to ",
                orig.weighted.treated.nobs, immediate.=Matchby.call)
        Var.calc <- orig.weighted.treated.nobs
      }
    if(Var.calc > orig.weighted.control.nobs)
      {
        warning("'Var.calc' > the number of control obs: 'Var.calc' reset to ",
                orig.weighted.control.nobs,  immediate.=Matchby.call)
        Var.calc <- orig.weighted.control.nobs
      }
    
    if(orig.nobs > 20000 & version!="fast" & !Matchby.call)
      {
        warning("The version='fast' option is recommended for large datasets if speed is desired.  For additional speed, you may also consider using the ties=FALSE option.", immediate.=TRUE)
      }

    #check the restrict matrix input
    if(!is.null(restrict))
      {
        if(!is.matrix(restrict))
          stop("'restrict' must be a matrix of restricted observations rows and three columns: c(i,j restriction)")

        if(ncol(restrict)!=3 )
          stop("'restrict' must be a matrix of restricted observations rows and three columns: c(i,j restriction)")
      }

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
            max.diff <- abs(max(X)-min(X) + tolerance * 100)
            ecaliper <- matrix(max.diff, nrow=xvars, ncol=1)
          }
        
        for (i in 1:xvars)
          {
            if (exact[i])
              ecaliper[i] <- tolerance;
          }
      }

    if(replace==FALSE)
      {
        #replace==FALE, needs enough observation
        #ATT
        orig.weighted.control.nobs <- sum(weights[Tr!=1])
        if(estimand==0)
          {
            if (orig.weighted.treated.nobs > orig.weighted.control.nobs)
              {
                warning("replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.")
              }
          } else if(estimand==1) 
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
            } else 
        {
              #ATC
          if (orig.weighted.treated.nobs < orig.weighted.control.nobs)
            {
              warning("replace==FALSE, but there are more (weighted) control obs than treated obs.  Some obs will be dropped.  You may want to estimate ATC instead")
            }
        }
        
        #we need a restrict matrix if we are going to not do replacement
        if(is.null(restrict))
          {
            restrict <- t(as.matrix(c(0,0,0)))
          }
        if(version!="fast" &  version!="standard")
          {
            warning("reverting to 'standard' version because replace=FALSE")
            version="standard"
          }
      }#end of replace==FALSE

#    if(version=="fast" & is.null(ecaliper) & sum(weights==1)==orig.nobs)
    if(version=="fast" | version=="standard")
      {
        if(!is.null(match.out))
          {
            ret <- RmatchLoop(Y=Y, Tr=Tr, X=X, Z=Z, V=V, All=estimand, M=M, BiasAdj=BiasAdj,
                              Weight=Weight, Weight.matrix=Weight.matrix, Var.calc=Var.calc,
                              weight=weights, SAMPLE=sample, ccc=ccc, cdd=cdd,
                              ecaliper=ecaliper, exact=exact, caliper=caliper,
                              restrict=restrict, MatchLoopC.indx=match.out$MatchLoopC,
                              weights.flag=weights.flag,  replace=replace, ties=ties,
                              version=version, MatchbyAI=MatchbyAI)
          } else {
            ret <- RmatchLoop(Y=Y, Tr=Tr, X=X, Z=Z, V=V, All=estimand, M=M, BiasAdj=BiasAdj,
                              Weight=Weight, Weight.matrix=Weight.matrix, Var.calc=Var.calc,
                              weight=weights, SAMPLE=sample, ccc=ccc, cdd=cdd,
                              ecaliper=ecaliper, exact=exact, caliper=caliper,
                              restrict=restrict, weights.flag=weights.flag,  replace=replace, ties=ties,
                              version=version, MatchbyAI=MatchbyAI)   
          }
      } else {
        ret <- Rmatch(Y=Y, Tr=Tr, X=X, Z=Z, V=V, All=estimand, M=M, BiasAdj=BiasAdj,
                      Weight=Weight, Weight.matrix=Weight.matrix, Var.calc=Var.calc,
                      weight=weights, SAMPLE=sample, ccc=ccc, cdd=cdd,
                      ecaliper=ecaliper, restrict=restrict)
      }

    if(is.null(ret$est))
      {
        if(!Matchby.call)
          {
            if(ret$valid < 1)
              {
                if (ret$sum.caliper.drops > 0) {
                  warning("'Match' object contains no valid matches (probably because of the caliper or the exact option).") 
                } else {
                  warning("'Match' object contains no valid matches")
                }
              } else {
                if (ret$sum.caliper.drops > 0) {
                  warning("'Match' object contains only 1 valid match (probably because of the caliper or the exact option).") 
                } else {
                  warning("'Match' object contains only one valid match")
                }            
              }
          } #endof if(!Matchby.call)

        z <- NA
        class(z)  <- "Match"    
        return(z)
      }
    
    indx <-  cbind(ret$art.data[,1],  ret$art.data[,2],  ret$W)

    index.treated  <- indx[,1]
    index.control  <- indx[,2]
    weights        <- indx[,3]
    sum.caliper.drops <- ret$sum.caliper.drops

   #RESET INDEX.TREATED        
    indx  <- as.matrix(cbind(index.treated,index.control))
    if (estimand==0) {
      #"ATT"
      index.treated  <- indx[,1]
      index.control  <- indx[,2]
    } else if(estimand==1) {
      #"ATE"
      tmp.index.treated  <- indx[,1]
      tmp.index.control  <- indx[,2]
      
      tl  <- length(tmp.index.treated)
      index.treated <- vector(length=tl, mode="numeric")
      index.control <- vector(length=tl, mode="numeric")
      trt  <- Tr[tmp.index.treated]==1
      for (i in 1:tl)
        {
          if (trt[i]) {
            index.treated[i]  <- tmp.index.treated[i]
            index.control[i]  <- tmp.index.control[i]
          } else {
            index.treated[i]  <- tmp.index.control[i]
            index.control[i]  <- tmp.index.treated[i]
          }
        }
    } else if(estimand==2) {
      #"ATC"
      index.treated  <- indx[,2]
      index.control  <- indx[,1]
    }        


    mdata  <- list()
    mdata$Y  <- c(Y[index.treated],Y[index.control])
    mdata$Tr <- c(Tr[index.treated],Tr[index.control])
    mdata$X  <- rbind(X[index.treated,],X[index.control,])
    mdata$orig.weighted.treated.nobs <- orig.weighted.treated.nobs

    #naive standard errors
    mest  <- sum((Y[index.treated]-Y[index.control])*weights)/sum(weights)
    v1  <- Y[index.treated] - Y[index.control]
    varest  <- sum( ((v1-mest)^2)*weights)/(sum(weights)*sum(weights))
    se.standard  <- sqrt(varest)

    wnobs <- sum(weights)

    if(estimand==0)
      {
        #ATT
        actual.drops <- orig.weighted.treated.nobs-wnobs
      } else if (estimand==1)
        {
          #ATE
          actual.drops <- orig.wnobs-wnobs
        } else {
          #ATC
          actual.drops <- (orig.wnobs-orig.weighted.treated.nobs)-wnobs 
        }

    #What obs were dropped?
    index.dropped <-  NULL #nothing was dropped
    if (sum.caliper.drops > 0 )
      {
        if(estimand.orig=="ATT")
          {
            matched.index <- which(Tr==1)
            matched <- !(matched.index %in% index.treated)        
            
          } else if(estimand.orig=="ATC")
            {
              matched.index <- which(Tr==0)
              matched <- !(matched.index %in% index.control)
            } else if(estimand.orig=="ATE")
              {
                matched.index <- 1:length(Tr)
                matched <- !(matched.index %in% c(index.treated,index.control))
              } 
        index.dropped <- matched.index[matched]  #obs not matched
      } #end of sum.caliper.drops > 0
            
    z  <- list(est=ret$est, se=ret$se, est.noadj=mest, se.standard=se.standard,
               se.cond=ret$se.cond, 
               mdata=mdata, 
               index.treated=index.treated, index.control=index.control, index.dropped=index.dropped,
               weights=weights, orig.nobs=orig.nobs, orig.wnobs=orig.wnobs,
               orig.treated.nobs=orig.treated.nobs,
               nobs=nobs, wnobs=wnobs,
               caliper=caliper, ecaliper=ecaliper, exact=exact,
               ndrops=actual.drops, ndrops.matches=sum.caliper.drops, 
               MatchLoopC=ret$MatchLoopC, version=version,
               estimand=estimand.orig)

    if(MatchbyAI)
      {
        z$YCAUS   <- ret$YCAUS
        z$ZCAUS   <- ret$ZCAUS
        z$Kcount  <- ret$Kcount
        z$KKcount <- ret$KKcount
        z$Sigs    <- ret$Sigs
      }

    class(z)  <- "Match"    
    return(z)
  } #end of Match

summary.Match  <- function(object, ..., full=FALSE, digits=5)
  {
    if(!is.list(object)) {
      warning("'Match' object contains less than two valid matches.  Cannot proceed.")
      return(invisible(NULL))
    }
    
    if (class(object) != "Match") {
      warning("Object not of class 'Match'")
      return(invisible(NULL))
    }

    if(object$version!="fast")
      {
        cat("\n")
        cat("Estimate... ",format(object$est,digits=digits),"\n")
        cat("AI SE...... ",format(object$se,digits=digits),"\n")
        cat("T-stat..... ",format(object$est/object$se,digits=digits),"\n")
        cat("p.val...... ",format.pval((1-pnorm(abs(object$est/object$se)))*2,digits=digits),"\n")
        cat("\n")
      } else {
        cat("\n")
        cat("Estimate... ",format(object$est,digits=digits),"\n")
        cat("SE......... ",format(object$se.standard,digits=digits),"\n")
        cat("T-stat..... ",format(object$est/object$se.standard,digits=digits),"\n")
        cat("p.val...... ",format.pval((1-pnorm(abs(object$est/object$se.standard)))*2,digits=digits),"\n")
        cat("\n")        
      }

    if(full)
      {
        cat("Est noAdj.. ",format(object$est.noadj,digits=digits),"\n")        
        cat("SE......... ",format(object$se.standard,digits=digits),"\n")
        cat("T-stat..... ",format(object$est.noadj/object$se.standard,digits=digits),"\n")
        cat("p.val...... ",format.pval((1-pnorm(abs(object$est.noadj/object$se.standard)))*2,digits=digits),"\n")
        cat("\n")
      }


    if(object$orig.wnobs!=object$orig.nobs)
      cat("Original number of observations (weighted)... ", round(object$orig.wnobs, 3),"\n")
    cat("Original number of observations.............. ", object$orig.nobs,"\n")
    if(object$mdata$orig.weighted.treated.nobs!=object$orig.treated.nobs)
      cat("Original number of treated obs (weighted).... ", round(object$mdata$orig.weighted.treated.nobs, 3),"\n")
    if(object$estimand!="ATC")
      {
        cat("Original number of treated obs............... ", object$orig.treated.nobs,"\n")
      } else {
        cat("Original number of control obs............... ",
            object$orig.nobs-object$orig.treated.nobs,"\n")        
      }
    cat("Matched number of observations............... ", round(object$wnobs, 3),"\n")
    cat("Matched number of observations  (unweighted). ", length(object$index.treated),"\n")

    cat("\n")
    if(!is.null(object$exact))
      {
        cat("Number of obs dropped by 'exact' or 'caliper' ", object$ndrops.matches, "\n")
        if (object$ndrops.matches!=round(object$ndrops))
          cat("Weighted #obs dropped by 'exact' or 'caliper' ", round(object$ndrops, 3),"\n")        
        cat("\n")        
      }else if(!is.null(object$caliper))
      {
        cat("Caliper (SDs)........................................  ",object$caliper,"\n")            
        cat("Number of obs dropped by 'exact' or 'caliper' ", object$ndrops.matches, "\n")
        if (object$ndrops.matches!=round(object$ndrops))
          cat("Weighted #obs dropped by 'exact' or 'caliper' ", round(object$ndrops, 3),"\n")        
        cat("\n")
      }
    
    z <- list()
    class(z) <- "summary.Match"
    return(invisible(z))
  } #end of summary.Match

print.summary.Match <- function(x, ...)
  {
    invisible(NULL)
  }


Rmatch <- function(Y, Tr, X, Z, V, All, M, BiasAdj, Weight, Weight.matrix, Var.calc, weight,
                   SAMPLE, ccc, cdd, ecaliper=NULL, restrict=NULL)
  {
    sum.caliper.drops <- 0
    X.orig <- X

#are we using the restriction matrix?
    if(is.matrix(restrict)) {
      restrict.trigger <- TRUE
    } else {
      restrict.trigger <- FALSE
    }
    
# if SATC is to be estimated the treatment indicator is reversed  
    if (All==2) {
      Tr <- 1-Tr
    }

    

# check on the number of matches, to make sure the number is within the limits
# feasible given the number of observations in both groups.
    if (All==1)
      {
        M <- min(M,min(sum(Tr),sum(1-Tr)));        
      } else {
        M <- min(M,sum(1-Tr));
      }

# two slippage parameters that are used to determine whether distances are equal
# distances less than ccc or cdd are interpeted as zero.
# these are passed in.  ccc, cdd


# I. set up
# I.a. vector for which observations we want the average effect
# iot_t is the vector with weights in the average treatment effects
# iot_c is the vector of indicators for potential controls

    if (All==1)
      {
        iot.t <- weight;
        iot.c <- as.matrix(rep(1,length(Tr)))
      } else {
        iot.t <- Tr*weight;
        iot.c <- 1-Tr    
      }

# I.b. determine sample and covariate vector sizes
    N  <- nrow(X)
    Kx <- ncol(X)
    Kz <- ncol(Z)

# K covariates
# N observations
#    Nt <- sum(Tr)
#    Nc <- sum(1-Tr)
#    on <- as.matrix(rep(1,N))

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
# The matrix AA enables one to transform the user supplied weight matrix 
# to take account of this transformation.  BUT THIS IS NOT USED!!
# Mu_X and Sig_X keep track of the original mean and variances
#    AA    <- diag(Kx)
    Mu.X  <- matrix(0, Kx, 1)
    Sig.X <- matrix(0, Kx, 1)

    for (k in 1:Kx)
      {
        Mu.X[k,1] <- sum(X[,k]*weight)/sum(weight)
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weight)/sum(weight)-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        if(Sig.X[k,1] < ccc)
          Sig.X[k,1] <- ccc
        X[,k]=eps/Sig.X[k,1]
#        AA[k,k]=Sig.X[k,1]
      } #end of k loop

#    Nv <- nrow(V)
    Mv <- ncol(V)
    Mu.V  <- matrix(0, Mv, 1)
    Sig.V <- matrix(0, Mv, 1)

    for (j in 1:Mv)
      {
        Mu.V[j,1]= ( t(V[,j])%*%weight ) /sum(weight)
#        dv <- V[,j]-Mu.V[j,1]
        sv <- sum(V[,j]*V[,j]*weight)/sum(weight)-Mu.V[j,1]^2
        if (sv > 0)
          {
            sv <- sqrt(sv)
          } else {
            sv <- 0
          }
        sv <- sv * sqrt(N/(N-1))
        Sig.V[j,1] <- sv
      } #end of j loop

# I.d. define weight matrix for metric, taking into account normalization of
# regressors.
# If the eigenvalues of the regressors are too close to zero the Mahalanobis metric
# is not used and we revert back to the default inverse of variances.
    if (Weight==1)
      {
        Weight.matrix=diag(Kx)
      } else if (Weight==2) {
        if (min (eigen( t(X)%*%X/N, only.values=TRUE)$values) > ccc)
          {
            Weight.matrix= solve(t(X)%*%X/N) 
          } else {
            Weight.matrix <- diag(Kx)
          }
      }
      # DO NOT RESCALE THE Weight.matrix!!
      #else if (Weight==3)
      #  {
      #    Weight.matrix <- AA %*% Weight.matrix %*% AA
      #  }

#    if (exact==1)
#      {
#        Weight.matrix <- cbind(Weight.matrix, matrix(0,nrow=Kx,ncol=Mv))
#        Weight.matrix <- rbind(Weight.matrix, cbind(matrix(0,nrow=Mv,ncol=Kx),
#                               1000*solve(diag(as.vector(Sig.V*Sig.V), nrow=length(Sig.V)))))
#        Weight.matrix <- as.matrix(Weight.matrix)
#        X <- cbind(X,V)
#        Mu.X  <- rbind(Mu.X, matrix(0, nrow(Mu.V), 1))
#        Sig.X <- rbind(Sig.X, matrix(1, nrow(Sig.V), 1))
#      } #end of exact

    Nx <- nrow(X)
    Kx <- ncol(X)

    if ( min(eigen(Weight.matrix, only.values=TRUE)$values) < ccc )
      Weight.matrix <- Weight.matrix + diag(Kx)*ccc

    # I.fg. initialize matrices before looping through sample
    YCAUS  <- as.matrix(rep(0, N))
    SCAUS  <- as.matrix(rep(0, N))
    XCAUS  <- matrix(0, nrow=N, ncol=Kx)
    ZCAUS  <- matrix(0, nrow=N, ncol=Kz)
    Kcount <- as.matrix(rep(0, N))
    KKcount <- as.matrix(rep(0, N))
    MMi     <- as.matrix(rep(0, N))

# II. Loop through all observations that need to be matched.
    INN <- as.matrix(1:N)
    ww <- chol(Weight.matrix) # so that ww*ww=w.m
#    TT <- as.matrix(1:N)

    # initialize some data objects
    DD <- NULL
    I  <- NULL
    IM <- NULL
    IT <- NULL
    IX <- NULL
    IZ <- NULL
    IY <- NULL
    W  <- NULL
    ADist <- NULL
    WWi <- NULL
    Xt <- NULL
    Zt <- NULL
    Yt <- NULL
    Xc <- NULL
    Zc <- NULL
    Yc <- NULL

    for (i in 1:N)
      {
        #treatment indicator for observation to be matched        
        TREATi <- Tr[i]

        # proceed with all observations if All==1
        # but only with treated observations if All=0        
        if ( (TREATi==1 & All!=1) | All==1 )
          {
            # covariate value for observation to be matched                        
            xx <- t(as.matrix(X[i,]))

            # covariate value for observation to be matched            
            zz <- t(as.matrix(Z[i,]))

            # outcome value for observation to be matched                        
            yy <- Y[i]

            #JSS: check *
            foo <- as.matrix(rep(1, Nx))
            DX <- (X - foo %*% xx) %*% t(ww)

            if (Kx>1)
              {
                #JSS
                foo <- t(DX*DX)
                Dist <- as.matrix(apply(foo, 2, sum))
              } else {
                Dist <- as.matrix(DX*DX)
              } #end of Kx

            # Dist distance to observation to be matched
            # is N by 1 vector

            #use the restriction matrix
            if (restrict.trigger)
              {
                for (j in 1:nrow(restrict))
                  {
                    if (restrict[j,1]==i)
                      {
                        if (restrict[j,3] < 0) {
                          Dist[restrict[j,2]] = .Machine$double.xmax
                        } else {
                          Dist[restrict[j,2]] = restrict[j,3]
                        }
                      } else if (restrict[j,2]==i) {
                        if (restrict[j,3] < 0) {
                          Dist[restrict[j,1]] = .Machine$double.xmax
                        } else {
                          Dist[restrict[j,1]] = restrict[j,3]
                        }
                      }
                  }
              } #end if restrict.trigger

            # set of potential matches (all observations with other treatment)
            # JSS, note:logical vector
            POTMAT <- Tr == (1-TREATi) 

            # X's for potential matches
#            XPOT <- X[POTMAT,]
            DistPot <- Dist[POTMAT,1]

#            TTPotMat <- TT[POTMAT,1]
            weightPot <- as.matrix(weight[POTMAT,1])

            # sorted distance of potential matches            
            S <- sort(DistPot)
            L <- order(DistPot)
            weightPot.sort <- weightPot[L,1]
            weightPot.sum <- cumsum(weightPot.sort)
            tt <- 1:(length(weightPot.sum))
            MMM <- min(tt[weightPot.sum>=M])

            # distance at last match
            Distmax <- S[MMM]

            # selection of actual matches 
            #logical index
            if (restrict.trigger)
              {
                ACTMAT <- POTMAT & ( (Dist <= (Distmax+cdd))  & (Dist < .Machine$double.xmax) )

                if (sum(ACTMAT) < 1)
                  next;                
              } else {
                ACTMAT <- POTMAT & ( Dist <= (Distmax+cdd) )
              }
            
            Ii <- i * matrix(1, nrow=sum(ACTMAT), ncol=1)
            IMi <- as.matrix(INN[ACTMAT,1])
            
            if(!is.null(ecaliper))
              {
                for (j in 1:length(Ii))
                  {
                    for( x in 1:Kx)
                      {
                        diff <- abs(X.orig[i,x] - X.orig[IMi[j], x])
                        if (diff > ecaliper[x])
                          {
#                            print(diff)

                            ACTMAT[IMi[j]] <- FALSE
                            sum.caliper.drops <- sum.caliper.drops+1
                            
                            break
                          }
                      } #x loop
                  } #j loop

                if (sum(ACTMAT) < 1)
                  next;
              } #ecaliper check

            # distance to actual matches            
            ACTDIST <- as.matrix(Dist[ACTMAT,1])

            # counts how many times each observation is matched.
            Kcount <- Kcount + weight[i] * weight*ACTMAT/sum(ACTMAT*weight)
            
            KKcount <- KKcount+weight[i,1] * weight*weight*ACTMAT /
              (sum(ACTMAT*weight)*sum(ACTMAT*weight))

            # counts how many times each observation is matched.
            # counts number of matches for observation i
            # Unless there are ties this should equal M

            Mi <- sum(weight*ACTMAT)            
            MMi[i,1] <- Mi
            Wi <- as.matrix(weight[ACTMAT,1])

            # mean of Y's for actual matches
            ymat <- t(Y[ACTMAT,1]) %*% Wi/Mi

            # mean of X's for actual matches
            # mean of Z's for actual matches
            if (length(Wi)>1)
              {
                xmat <- t(t(X[ACTMAT,]) %*% Wi/Mi)
                zmat <- t(t(Z[ACTMAT,]) %*% Wi/Mi)
              } else {
                xmat <- t(X[ACTMAT,]) * as.double(Wi)/Mi
                zmat <- t(Z[ACTMAT,]) * as.double(Wi)/Mi
              }

            # estimate causal effect on y for observation i
            YCAUS[i,1] <- TREATi %*% (yy-ymat)+(1-TREATi) %*% (ymat-yy)

            # difference between x and actual matches for observation i
            XCAUS[i,] <- TREATi %*% (xx-xmat)+(1-TREATi) %*% (xmat-xx)

            ZCAUS[i,] <- TREATi %*% (zz-zmat)+(1-TREATi) %*% (zmat-zz)

            # collect results
            I  <- rbind(I, i * matrix(1, nrow=sum(ACTMAT), ncol=1))
            DD <- rbind(DD, ACTDIST)
            IM <- rbind(IM, as.matrix(INN[ACTMAT,1]))
            IT <- rbind(IT, TREATi * as.matrix(rep(1, sum(ACTMAT))))
            IX <- rbind(IX, as.matrix(rep(1, sum(ACTMAT))) %*% xx)
            IZ <- rbind(IZ, as.matrix(rep(1, sum(ACTMAT))) %*% zz)
            IY <- rbind(IY, as.matrix(rep(1, sum(ACTMAT))) * yy)

            # weight for matches
            W <- as.matrix(c(W, weight[i,1] * Wi/Mi))
            ADist <- as.matrix(c(ADist, ACTDIST))
            WWi   <- as.matrix(c(WWi, Wi))

            if (TREATi==1)
              {

                if (ncol(X) > 1)
                  {
                    # covariates for treated
                    Xt <- rbind(Xt, as.matrix(rep(1, sum(ACTMAT))) %*% xx)                    
                    # covariate for matches
                    Xc <- rbind(Xc, X[ACTMAT,])
                    
                  } else {
                    # covariates for treated
                    Xt <- as.matrix(c(Xt, as.matrix(rep(1, sum(ACTMAT))) %*% xx))
                    # covariate for matches
                    Xc <- as.matrix(c(Xc, X[ACTMAT,]))
                  }

                if (ncol(Z) > 1)
                  {
                    # covariates for treated                
                    Zt <- rbind(Zt, as.matrix(rep(1, sum(ACTMAT))) %*% zz)
                    # covariate for matches
                    Zc <- rbind(Zc, Z[ACTMAT,])
                  } else {
                    # covariates for treated                
                    Zt <- as.matrix(c(Zt, as.matrix(rep(1, sum(ACTMAT))) %*% zz))
                    # covariate for matches                
                    Zc <- as.matrix(c(Zc, Z[ACTMAT,]))
                  }

                # outcome for treated                              
                Yt <- as.matrix(c(Yt, yy * as.matrix(rep(1, sum(ACTMAT))) ))
                # outcome for matches
                Yc <- as.matrix(c(Yc, Y[ACTMAT,1]))
              } else {

                if (ncol(X) > 1)
                  {
                    # covariates for controls
                    Xc <- rbind(Xc, as.matrix(rep(1, sum(ACTMAT))) %*% xx)
                    # covariate for matches
                    Xt <- rbind(Xt, X[ACTMAT,])
                  } else {
                    # covariates for controls                
                    Xc <- as.matrix(c(Xc, as.matrix(rep(1, sum(ACTMAT))) %*% xx))
                    # covariate for matches
                    Xt <- as.matrix(c(Xt, X[ACTMAT,]))
                  }

                if (ncol(Z) > 1)
                  {
                    # covariates for controls                
                    Zc <- rbind(Zc, as.matrix(rep(1, sum(ACTMAT))) %*% zz)
                    # covariate for matches                
                    Zt <- rbind(Zt, Z[ACTMAT,])
                  } else {
                    # covariates for controls                
                    Zc <- as.matrix(c(Zc, as.matrix(rep(1, sum(ACTMAT))) %*% zz))
                    # covariate for matches                
                    Zt <- as.matrix(c(Zt, Z[ACTMAT,]))
                  }
                # outcome for controls                               
                Yc <- as.matrix(c(Yc, as.matrix(rep(1, sum(ACTMAT))) * yy))
                # outcome for matches
                Yt <- as.matrix(c(Yt, Y[ACTMAT,1]))

              } #end of TREATi
          } #end of if
      } #i loop

    # transform matched covariates back for artificial data set
    Xt.u <- Xt
    Xc.u <- Xc
    IX.u <- IX
    for (k in 1:Kx)
      {
        Xt.u[,k] <- Mu.X[k,1]+Sig.X[k,1] * Xt.u[,k]
        Xc.u[,k] <- Mu.X[k,1]+Sig.X[k,1] * Xc.u[,k]
        IX.u[,k] <- Mu.X[k,1]+Sig.X[k,1] * IX.u[,k]
      }
    if (All!=1)
      {
        I  <- as.matrix(I[IT==1,])
        IM <- as.matrix(IM[IT==1,])
        IT <- as.matrix(IT[IT==1,])
        IY <- as.matrix(IY[IT==1,])
        Yc <- as.matrix(Yc[IT==1,])
        Yt <- as.matrix(Yt[IT==1,])
        W  <- as.matrix(W[IT==1,])
        ADist <- as.matrix(ADist[IT==1,])
        WWi   <- as.matrix(WWi[IT==1,])
        IX.u  <- as.matrix(IX.u[IT==1,])
        Xc.u  <- as.matrix(Xc.u[IT==1,])
        Xt.u  <- as.matrix(Xt.u[IT==1,])
        Xc    <- as.matrix(Xc[IT==1,])
        Xt <- as.matrix(Xt[IT==1,])
        Zc <- as.matrix(Zc[IT==1,])
        Zt <- as.matrix(Zt[IT==1,])
        IZ <- as.matrix(IZ[IT==1,])
      } #end of if

    if (length(I) < 1)
      {
        return(list(sum.caliper.drops=sum.caliper.drops,valid=0))
      } else if(length(I) < 2)
        {
          return(list(sum.caliper.drops=sum.caliper.drops,valid=1))          
        }

    if (BiasAdj==1)
      {
        # III. Regression of outcome on covariates for matches
        if (All==1)
          {
            # number of observations
            NNt <- nrow(Z)
            # add intercept        
            ZZt <- cbind(rep(1, NNt), Z)
            # number of covariates
            Nx <- nrow(ZZt)
            Kx <- ncol(ZZt)
            xw <- ZZt*(sqrt(Tr*Kcount) %*% t(as.matrix(rep(1,Kx))))
            
            foo <- min(eigen(t(xw)%*%xw, only.values=TRUE)$values)
            foo <- as.double(foo<=ccc)
            foo2 <- apply(xw, 2, sd)

            options(show.error.messages = FALSE)
            wout <- NULL
            try(wout <- solve( t(xw) %*% xw + diag(Kx) * ccc * (foo) * max(foo2)) %*%
                (t(xw) %*% (Y*sqrt(Tr*Kcount))))
            if(is.null(wout))
              {
                wout2 <- NULL
                try(wout2 <- ginv( t(xw) %*% xw + diag(Kx) * ccc * (foo) * max(foo2)) %*%
                    (t(xw) %*% (Y*sqrt(Tr*Kcount))))
                if(!is.null(wout2))
                  {
                    wout <-wout2
                    warning("using generalized inverse to calculate Bias Adjustment probably because of singular 'Z'")
                  }
              }
            options(show.error.messages = TRUE)
            if(is.null(wout))
              {
                warning("unable to calculate Bias Adjustment probably because of singular 'Z'")
                BiasAdj <- 0
              } else {
                NW <- nrow(wout)
#                KW <- ncol(wout)
                Alphat <- wout[2:NW,1]
              }
          } else {
            Alphat <- matrix(0, nrow=Kz, ncol=1)
          } #end if ALL
      }

    if(BiasAdj==1)
      {
        # III.b.  Controls
        NNc <- nrow(Z)
        ZZc <- cbind(matrix(1, nrow=NNc, ncol=1),Z)
        Nx <- nrow(ZZc)
        Kx <- ncol(ZZc)
        
        xw <- ZZc*(sqrt((1-Tr)*Kcount) %*% matrix(1, nrow=1, ncol=Kx))
        
        foo <- min(eigen(t(xw)%*%xw, only.values=TRUE)$values)
        foo <- as.double(foo<=ccc)
        foo2 <- apply(xw, 2, sd)

        options(show.error.messages = FALSE)
        wout <- NULL        
        try(wout <- solve( t(xw) %*% xw + diag(Kx) * ccc * (foo) * max(foo2)) %*%
            (t(xw) %*% (Y*sqrt((1-Tr)*Kcount))))
        if(is.null(wout))
          {
            wout2 <- NULL
            try(wout2 <- ginv( t(xw) %*% xw + diag(Kx) * ccc * (foo) * max(foo2)) %*%
                (t(xw) %*% (Y*sqrt((1-Tr)*Kcount))))
            if(!is.null(wout2))
              {
                wout <-wout2
                warning("using generalized inverse to calculate Bias Adjustment probably because of singular 'Z'")
              }
          }        
        options(show.error.messages = TRUE)
        if(is.null(wout))
          {
            warning("unable to calculate Bias Adjustment probably because of singular 'Z'")
            BiasAdj <- 0
          } else {
            NW <- nrow(wout)
#            KW <- ncol(wout)
            Alphac <- as.matrix(wout[2:NW,1])
          }
      }
    

    if(BiasAdj==1)
      {
        # III.c. adjust matched outcomes using regression adjustment for bias adjusted matching estimator

        SCAUS <- YCAUS-Tr*(ZCAUS %*% Alphac)-(1-Tr)*(ZCAUS %*% Alphat)
        # adjusted control outcome
        Yc.adj <- Yc+BiasAdj * (IZ-Zc) %*% Alphac
        # adjusted treated outcome
        Yt.adj <- Yt+BiasAdj*(IZ-Zt) %*% Alphat
        Tau.i <- Yt.adj - Yc.adj
      } else {
        Yc.adj <- Yc
        Yt.adj <- Yt
        Yt.adj <- Yt
        Tau.i <- Yt.adj - Yc.adj        
      }

    art.data <- cbind(I,IM,IT,DD,IY,Yc,Yt,W,WWi,ADist,IX.u,Xc.u,Xt.u,
                      Yc.adj,Yt.adj,Tau.i)


    # III. If conditional variance is needed, initialize variance vector
    # and loop through all observations

    Nx <- nrow(X)
    Kx <- ncol(X)

#   ww <- chol(Weight.matrix)
#   NN <- as.matrix(1:N)
    if (Var.calc>0)
      {
        Sig <- matrix(0, nrow=N, ncol=1)
        # overall standard deviation of outcome
        # std <- sd(Y)
        for (i in 1:N)
          {
            # treatment indicator observation to be matched
            TREATi <- Tr[i,1]
            # covariate value for observation to be matched
            xx <- X[i,]
            # potential matches are all observations with the same treatment value
            POTMAT <- (Tr==TREATi)
            POTMAT[i,1] <- 0
            weightPOT <- as.matrix(weight[POTMAT==1,1])
            DX <- (X - matrix(1, Nx,1) %*% xx) %*% t(ww)
            if (Kx>1)
              {
                foo <- apply(t(DX*DX), 2, sum)
                Dist <- as.matrix(foo)
              } else {
                Dist <- DX*DX
              }

            # distance to observation to be matched

            # Distance vector only for potential matches
            DistPot <- Dist[POTMAT==1,1]
            # sorted distance of potential matches
            S <- sort(DistPot)
            L <- order(DistPot)
            weightPOT.sort <- weightPOT[L,1]
            weightPOT.sum <- cumsum(weightPOT.sort)
            tt <-  1:(length(weightPOT.sum))
            MMM <- min(tt[weightPOT.sum >= Var.calc])
            MMM <- min(MMM,length(S))
            Distmax=S[MMM]

            # distance of Var_calc-th closest match
            ACTMAT <- (POTMAT==1) & (Dist<= (Distmax+ccc))

            # indicator for actual matches, that is all potential
            # matches closer than, or as close as the Var_calc-th
            # closest

            Yactmat <- as.matrix(c(Y[i,1], Y[ACTMAT,1]))
            weightactmat <- as.matrix(c(weight[i,1], weight[ACTMAT,1]))
            fm <- t(Yactmat) %*% weightactmat/sum(weightactmat)
            sm <- sum(Yactmat*Yactmat*weightactmat)/sum(weightactmat)
            sigsig <- (sm-fm %*% fm)*sum(weightactmat)/(sum(weightactmat)-1)
            
            # standard deviation of actual matches
            Sig[i,1] <- sqrt(sigsig)
          }# end of i loop
        #variance estimate 
        Sigs <- Sig*Sig
      } #end of var.calc > 0

    # matching estimator
    est <- t(W) %*% Tau.i/sum(W)
#    est.t <- sum((iot.t*Tr+iot.c*Kcount*Tr)*Y)/sum(iot.t*Tr+iot.c*Kcount*Tr)
#    est.c <- sum((iot.t*(1-Tr)+iot.c*Kcount*(1-Tr))*Y)/sum(iot.t*(1-Tr)+iot.c*Kcount*(1-Tr))

    if (Var.calc==0)
      {
        eps <- Tau.i - as.double(est)
        eps.sq <- eps*eps
        Sigs <- 0.5 * matrix(1, N, 1) %*% (t(eps.sq) %*% W)/sum(W)
#        sss <- sqrt(Sigs[1,1])
      } #end of Var.calc==0

    SN <- sum(iot.t)
    var.sample <- sum((Sigs*(iot.t+iot.c*Kcount)*(iot.t+iot.c*Kcount))/(SN*SN))

    if (All==1)
      {
        var.pop <- sum((Sigs*(iot.c*Kcount*Kcount+2*iot.c*Kcount-iot.c*KKcount))/(SN*SN))
      } else {
        var.pop=sum((Sigs*(iot.c*Kcount*Kcount-iot.c*KKcount))/(SN*SN))
      }

    if (BiasAdj==1)
      {
        dvar.pop <- sum(iot.t*(SCAUS-as.double(est))*(SCAUS-as.double(est)))/(SN*SN)
      } else {
        dvar.pop <- sum(iot.t*(YCAUS-as.double(est))*(YCAUS-as.double(est)))/(SN*SN)
      }

    var.pop <- var.pop + dvar.pop

    if (SAMPLE==1)
      {
        var <- var.sample
      } else {
        var <- max(var.sample, var.pop)
        var <- var.pop
      }

    var.cond <- max(var.sample,var.pop)-var.sample
    se <- sqrt(var)
    se.cond <- sqrt(var.cond)

#    Sig <- sqrt(Sigs)
#    aug.data <- cbind(Y,Tr,X,Z,Kcount,Sig,weight)

    if (All==2)
      est <- -1*est

#    if (exact==1)
#      {
#        Vt.u <- Xt.u[,(Kx-Mv+1):Kx]
#        Vc.u <- Xc.u[,(Kx-Mv+1):Kx]
#        Vdif <- abs(Vt.u-Vc.u)
#
#        if (Mv>1)
#          Vdif <- as.matrix(apply(t(Vdif), 2, sum))
#
#        em[1,1] <- length(Vdif)
#        em[2,1] <- sum(Vdif>0.000001)
#      }#end of exact==1

    return(list(est=est, se=se, se.cond=se.cond, W=W,
                sum.caliper.drops=sum.caliper.drops,
                art.data=art.data))
  }# end of Rmatch


#
# See Rosenbaum and Rubin (1985) and Smith and Todd Rejoinder (JofEconometrics) p.9
#
sdiff.pooled  <- function(Tr, Co, weights=rep(1,length(Co)),
                          weights.Tr=rep(1,length(Tr)),
                          weights.Co=rep(1,length(Co)),
                          match=FALSE)
  {
    if (!match)
      {
        obs.Tr <- sum(weights.Tr)
        obs.Co <- sum(weights.Co)
#        obs.total <- obs.Tr+obs.Co
        
        mean.Tr <- sum(Tr*weights.Tr)/obs.Tr
        mean.Co <- sum(Co*weights.Co)/obs.Co
        diff  <- mean.Tr - mean.Co

#match Rubin
#        mean.total <- sum(Tr*weights.Tr)/obs.total + sum(Co*weights.Co)/obs.total
#        var.total  <- sum( ( (Tr - mean.total)^2 )*weights.Tr)/(obs.total-1) +
#          sum( ( (Co - mean.total)^2 )*weights.Co)/(obs.total-1)
        var.pooled <- ( sum( ( (Tr - mean.Tr)^2)*weights.Tr)/(obs.Tr-1) +
          sum( ( (Co - mean.Co)^2 )*weights.Co)/(obs.Co-1) )/2 
        

        if(var.pooled==0 & diff==0)
          {
            sdiff <- 0
          } else {
            sdiff <- 100*diff/sqrt( var.pooled )
          }
        
      } else{
        diff  <- sum( (Tr-Co)*weights ) /sum(weights)
        mean.Tr <- sum(Tr*weights)/sum(weights)
        mean.Co <- sum(Co*weights)/sum(weights)

#match Rubin
        obs <- sum(weights)

#       obs for total        
#       obs = sum(weights)*2
#        mean.total <- (mean.Tr + mean.Co)/2
#        var.total  <- sum( ( (Tr - mean.total)^2 )*weights)/(obs-1) +
#          sum( ( (Co - mean.total)^2 )*weights)/(obs-1)
        
        var.pooled  <- ( sum( ( (Tr - mean.Tr)^2 )*weights)/(obs-1) +
          sum( ( (Co - mean.Co)^2 )*weights)/(obs-1) )/2

        if(var.pooled==0 & diff==0)
          {
            sdiff <- 0
          } else {
            sdiff <- 100*diff/sqrt(var.pooled)
          }        
      }

    return(sdiff)
  }


#
# STANDARDIZED BY THE VARIANCE OF THE TREATMENT GROUP
# See sdiff.rubin for Rosenbaum and Rubin (1985) and Smith and Todd Rejoinder (JofEconometrics) p.9
#
sdiff  <- function(Tr, Co, weights=rep(1,length(Co)),
                   weights.Tr=rep(1,length(Tr)),
                   weights.Co=rep(1,length(Co)),
                   match=FALSE,
                   estimand="ATT")

  {
    if (!match)
      {
        obs.Tr <- sum(weights.Tr)
        obs.Co <- sum(weights.Co)
        
        mean.Tr <- sum(Tr*weights.Tr)/obs.Tr
        mean.Co <- sum(Co*weights.Co)/obs.Co
        diff  <- mean.Tr - mean.Co

        if(estimand=="ATC")
          {
            var  <- sum( ( (Co - mean.Co)^2 )*weights.Co)/(obs.Co-1)
          } else {
            var  <- sum( ( (Tr - mean.Tr)^2 )*weights.Tr)/(obs.Tr-1)
          }

        if(var==0 & diff==0)
          {
            sdiff=0
          } else {
            sdiff <- 100*diff/sqrt( (var) )
          }
      } else{
        mean.Tr <- sum(Tr*weights)/sum(weights)
        mean.Co <- sum(Co*weights)/sum(weights)
        diff  <- mean.Tr - mean.Co

        if(estimand=="ATC")
          {
            var  <- sum( ( (Co - mean.Co)^2 )*weights)/(sum(weights)-1)
          } else {
            var  <- sum( ( (Tr - mean.Tr)^2 )*weights)/(sum(weights)-1)
          }

        if(var==0 & diff==0)
          {
            sdiff <- 0
          } else {
            sdiff <- 100*diff/sqrt( (var) )
          }        
      }

    return(sdiff)
  }

# function runs sdiff and t.test
#
balanceUV  <- function(Tr, Co, weights=rep(1,length(Co)), exact=FALSE, ks=FALSE, nboots=1000,
                       paired=TRUE, match=FALSE,
                       weights.Tr=rep(1,length(Tr)), weights.Co=rep(1,length(Co)),
                       estimand="ATT")
  {
    ks.out  <- NULL

    if (!match)
      {
        sbalance.pooled  <- sdiff.pooled(Tr=Tr, Co=Co,
                                       weights.Tr=weights.Tr,
                                       weights.Co=weights.Co,
                                       match=FALSE)

        sbalance.constvar  <- sdiff(Tr=Tr, Co=Co,
                                    weights.Tr=weights.Tr,
                                    weights.Co=weights.Co,
                                    match=FALSE,
                                    estimand=estimand)        

        obs.Tr <- sum(weights.Tr)
        obs.Co <- sum(weights.Co)        
        
        mean.Tr <- sum(Tr*weights.Tr)/obs.Tr
        mean.Co <- sum(Co*weights.Co)/obs.Co
        var.Tr  <- sum( ( (Tr - mean.Tr)^2 )*weights.Tr)/(obs.Tr-1)
        var.Co  <- sum( ( (Co - mean.Co)^2 )*weights.Co)/(obs.Co-1)
        var.ratio  <- var.Tr/var.Co

        qqsummary <- qqstats(x=Tr, y=Co, standardize=TRUE)
        qqsummary.raw <- qqstats(x=Tr, y=Co, standardize=FALSE)

        tt  <- Mt.test.unpaired(Tr, Co,
                                weights.Tr=weights.Tr,
                                weights.Co=weights.Co)

        if (ks)
          ks.out <- ks.boot(Tr, Co,nboots=nboots)
      } else {
        sbalance.pooled  <- sdiff(Tr=Tr, Co=Co, weights=weights, match=TRUE)
        sbalance.constvar  <- sdiff(Tr=Tr, Co=Co, weights=weights, match=TRUE, estimand=estimand)
        
        mean.Tr  <- sum(Tr*weights)/sum(weights);
        mean.Co  <- sum(Co*weights)/sum(weights);
        var.Tr  <- sum( ( (Tr - mean.Tr)^2 )*weights)/sum(weights);
        var.Co  <- sum( ( (Co - mean.Co)^2 )*weights)/sum(weights);
        var.ratio  <- var.Tr/var.Co

        qqsummary <- qqstats(x=Tr, y=Co, standardize=TRUE)
        qqsummary.raw <- qqstats(x=Tr, y=Co, standardize=FALSE)

        if(paired)
          {
            tt  <- Mt.test(Tr, Co, weights)
          } else {
            tt  <- Mt.test.unpaired(Tr, Co, weights.Tr=weights, weights.Co=weights)
          }
        
        if (ks)
          ks.out  <- ks.boot(Tr, Co, nboots=nboots)
      }        

    ret  <- list(sdiff=sbalance.constvar,
                 sdiff.pooled=sbalance.pooled,
                 mean.Tr=mean.Tr,mean.Co=mean.Co,
                 var.Tr=var.Tr,var.Co=var.Co, p.value=tt$p.value,
                 var.ratio=var.ratio,
                 ks=ks.out, tt=tt,
                 qqsummary=qqsummary,
                 qqsummary.raw=qqsummary.raw)

    class(ret) <-  "balanceUV"
    return(ret)
  } #balanceUV

summary.balanceUV  <- function(object, ..., digits=5)
  {

    if (class(object) != "balanceUV") {
      warning("Object not of class 'balanceUV'")
      return(NULL)
    }

    cat("mean treatment........", format(object$mean.Tr, digits=digits),"\n")
    cat("mean control..........", format(object$mean.Co, digits=digits),"\n")
    cat("std mean diff.........", format(object$sdiff, digits=digits),"\n\n")    

    cat("mean raw eQQ diff.....", format(object$qqsummary.raw$meandiff, digits=digits),"\n")
    cat("med  raw eQQ diff.....", format(object$qqsummary.raw$mediandiff, digits=digits),"\n")
    cat("max  raw eQQ diff.....", format(object$qqsummary.raw$maxdiff, digits=digits),"\n\n")

    cat("mean eCDF diff........", format(object$qqsummary$meandiff, digits=digits),"\n")
    cat("med  eCDF diff........", format(object$qqsummary$mediandiff, digits=digits),"\n")
    cat("max  eCDF diff........", format(object$qqsummary$maxdiff, digits=digits),"\n\n")

    cat("var ratio (Tr/Co).....", format(object$var.ratio, digits=digits),"\n")
    cat("T-test p-value........", format.pval(object$tt$p.value,digits=digits), "\n")            
    if (!is.null(object$ks))
      {
        if(!is.na(object$ks$ks.boot.pvalue))
          {
            cat("KS Bootstrap p-value..", format.pval(object$ks$ks.boot.pvalue, digits=digits), "\n")
          }
        cat("KS Naive p-value......", format(object$ks$ks$p.value, digits=digits), "\n")                        
        cat("KS Statistic..........", format(object$ks$ks$statistic, digits=digits), "\n")
      }
    cat("\n")
  } #end of summary.balanceUV


#print Before and After balanceUV together
PrintBalanceUV  <- function(BeforeBalance, AfterBalance, ..., digits=5)
  {

    if (class(BeforeBalance) != "balanceUV") {
      warning("BeforeBalance not of class 'balanceUV'")
      return(NULL)
    }

    if (class(AfterBalance) != "balanceUV") {
      warning("AfterBalance not of class 'balanceUV'")
      return(NULL)
    }

    space.size <- digits*2
#    space <- rep(" ",space.size)

    cat("                      ", "Before Matching", "\t \t After Matching\n")
    cat("mean treatment........", format(BeforeBalance$mean.Tr, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$mean.Tr, digits=digits, width=space.size), 
        "\n")
    cat("mean control..........", format(BeforeBalance$mean.Co, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$mean.Co, digits=digits, width=space.size),
        "\n")
    cat("std mean diff.........", format(BeforeBalance$sdiff, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$sdiff, digits=digits, width=space.size),
        "\n\n")
    
    cat("mean raw eQQ diff.....", format(BeforeBalance$qqsummary.raw$meandiff, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$qqsummary.raw$meandiff, digits=digits, width=space.size),
        "\n")
    cat("med  raw eQQ diff.....", format(BeforeBalance$qqsummary.raw$mediandiff, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$qqsummary.raw$mediandiff, digits=digits, width=space.size),
        "\n")
    cat("max  raw eQQ diff.....", format(BeforeBalance$qqsummary.raw$maxdiff, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$qqsummary.raw$maxdiff, digits=digits, width=space.size),
        "\n\n")

    cat("mean eCDF diff........", format(BeforeBalance$qqsummary$meandiff, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$qqsummary$meandiff, digits=digits, width=space.size),
        "\n")
    cat("med  eCDF diff........", format(BeforeBalance$qqsummary$mediandiff, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$qqsummary$mediandiff, digits=digits, width=space.size),
        "\n")
    cat("max  eCDF diff........", format(BeforeBalance$qqsummary$maxdiff, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$qqsummary$maxdiff, digits=digits, width=space.size),
        "\n\n")

    cat("var ratio (Tr/Co).....", format(BeforeBalance$var.ratio, digits=digits, width=space.size), "\t \t",
        format(AfterBalance$var.ratio, digits=digits, width=space.size),
        "\n")
    cat("T-test p-value........", format(format.pval(BeforeBalance$tt$p.value,digits=digits), justify="right", width=space.size), "\t \t",
        format(format.pval(AfterBalance$tt$p.value,digits=digits), justify="right", width=space.size),
        "\n")            
    if (!is.null(BeforeBalance$ks))
      {
        if(!is.na(BeforeBalance$ks$ks.boot.pvalue))
          {
            cat("KS Bootstrap p-value..", format(format.pval(BeforeBalance$ks$ks.boot.pvalue, digits=digits),  justify="right",width=space.size), "\t \t",
                format(format.pval(AfterBalance$ks$ks.boot.pvalue, digits=digits),  justify="right", width=space.size),
                "\n")
          }
        cat("KS Naive p-value......", format(format.pval(BeforeBalance$ks$ks$p.value, digits=digits),  justify="right",width=space.size), "\t \t",
            format(format.pval(AfterBalance$ks$ks$p.value, digits=digits),  justify="right",width=space.size),
            "\n")                        
        cat("KS Statistic..........", format(BeforeBalance$ks$ks$statistic, digits=digits, width=space.size), "\t \t",
            format(AfterBalance$ks$ks$statistic, digits=digits, width=space.size),
            "\n")
      }
    cat("\n")        
  } #end of PrintBalanceUV


#removed as of 0.99-7 (codetools)
#McNemar  <- function(Tr, Co, weights=rep(1,length(Tr)))

McNemar2 <- function (Tr, Co, correct = TRUE, weights=rep(1,length(Tr)))
{
  x  <- Tr
  y  <- Co
  if (is.matrix(x)) {
    stop("this version of McNemar cannot handle x being a matrix")
  } else {
    if (is.null(y)) 
      stop("if x is not a matrix, y must be given")
    if (length(x) != length(y)) 
      stop("x and y must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- factor(x[OK])
    y <- factor(y[OK])
    r <- nlevels(x)
    if ((r < 2) || (nlevels(y) != r))
      {
        stop("x and y must have the same number of levels (minimum 2)")
      }
  }
  tx <- table(x, y)
  facs  <- levels(x)
  txw  <- tx
  for(i in 1:r)
    {
      for(j in 1:r)
        {
          indx  <- x==facs[i] & y==facs[j]
          txw[i,j]  <- sum(weights[indx]);
        }
    }
  pdiscordant  <- sum( ( (x!=y)*weights )/sum(weights) )
  x  <- txw
  
  PARAMETER <- r * (r - 1)/2
  METHOD <- "McNemar's Chi-squared test"
  if (correct && (r == 2) && any(x - t(x))) {
    y <- (abs(x - t(x)) - 1)
    METHOD <- paste(METHOD, "with continuity correction")
  } else y <- x - t(x)
  x <- x + t(x)
  STATISTIC <- sum(y[upper.tri(x)]^2/x[upper.tri(x)])
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  names(STATISTIC) <- "McNemar's chi-squared"
  names(PARAMETER) <- "df"

  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
               p.value = PVAL, method = METHOD, data.name = DNAME,
               pdiscordant=pdiscordant)
  class(RVAL) <- "htest"
  return(RVAL)
}

ks<-function(x,y,w=F,sig=T){
#  Compute the Kolmogorov-Smirnov test statistic
#
# Code for the Kolmogorov-Smirnov test is adopted from the Splus code
# written by Rand R. Wilcox for his book titled "Introduction to
# Robust Estimation and Hypothesis Testing." Academic Press, 1997.
#
#  Also see 
#@book( knuth1998,
#  author=       {Knuth, Donald E.},
#  title=        {The Art of Computer Programming, Vol. 2: Seminumerical Algorithms},
#  edition=      "3rd",
#  publisher=    "Addison-Wesley", address= "Reading: MA", 
#  year=         1998
#)
# and Wilcox 1997 
#
#  w=T computes the weighted version instead. 
#  sig=T indicates that the exact significance level is to be computed.
#  If there are ties, the reported significance level is exact when
#  using the unweighted test, but for the weighted test the reported
#  level is too high.
#
#  This function uses the functions ecdf, kstiesig, kssig and kswsig
#
#  This function returns the value of the test statistic, the approximate .05
#  critical value, and the exact significance level if sig=T.
#
#  Missing values are automatically removed
#

  ecdf<-function(x,val){
  #  compute empirical cdf for data in x evaluated at val
  #  That is, estimate P(X <= val) 
  #
    ecdf<-length(x[x<=val])/length(x)
    ecdf
  }#end ecdf
  
  x<-x[!is.na(x)]
  y<-y[!is.na(y)]
  w<-as.logical(w)
  sig<-as.logical(sig)
  tie<-logical(1)
  tie<-F
  siglevel<-NA
  z<-sort(c(x,y))  # Pool and sort the observations
  for (i in 2:length(z))if(z[i-1]==z[i])tie<-T #check for ties
  v<-1   # Initializes v
  for (i in 1:length(z))v[i]<-abs(ecdf(x,z[i])-ecdf(y,z[i]))
  ks<-max(v)
  crit<-1.36*sqrt((length(x)+length(y))/(length(x)*length(y))) # Approximate
#                                                       .05 critical value
  if(!w && sig && !tie)siglevel<-kssig(length(x),length(y),ks)
  if(!w && sig && tie)siglevel<-kstiesig(x,y,ks)
  if(w){
    crit<-(max(length(x),length(y))-5)*.48/95+2.58+abs(length(x)-length(y))*.44/95
    if(length(x)>100 || length(y)>100){
      print("When either sample size is greater than 100,")
      print("the approximate critical value can be inaccurate.")
      print("It is recommended that the exact significance level be computed.")
    }
    for (i in 1:length(z)){
      temp<-(length(x)*ecdf(x,z[i])+length(y)*ecdf(y,z[i]))/length(z)
      temp<-temp*(1.-temp)
      v[i]<-v[i]/sqrt(temp)
    }
    v<-v[!is.na(v)]
    ks<-max(v)*sqrt(length(x)*length(y)/length(z))
    if(sig)siglevel<-kswsig(length(x),length(y),ks)
    if(tie && sig){
      print("Ties were detected. The reported significance level")
      print("of the weighted Kolmogorov-Smirnov test statistic is not exact.")
    }}
                                        #round off siglevel in a nicer way
  if(is.double(ks) & is.double(crit) & !is.na(ks) & !is.na(crit))
    {
      if (is.na(siglevel) & ks < crit)
        {
          siglevel  <- 0.99999837212332
        }
      if (is.double(siglevel) & !is.na(siglevel))
        {
          if (siglevel < 0)
            siglevel  <- 0
        }
    }
  list(test=ks,critval=crit,pval=siglevel)
}


kssig<-function(m,n,val){
#
#    Compute significance level of the  Kolmogorov-Smirnov test statistic
#    m=sample size of first group
#    n=sample size of second group
#    val=observed value of test statistic
#
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 0:m){
    for (j in 0:n)cmat[i+1,j+1]<-abs(i/m-j/n)
  }
  cmat<-ifelse(cmat<=val,1e0,0e0)
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kssig<-1.-umat[m+1,n+1]/exp(term)
  return(kssig)
}

kstiesig<-function(x,y,val){
#
#    Compute significance level of the  Kolmogorov-Smirnov test statistic
#    for the data in x and y.
#    This function allows ties among the  values.
#    val=observed value of test statistic
#
m<-length(x)
n<-length(y)
z<-c(x,y)
z<-sort(z)
cmat<-matrix(0,m+1,n+1)
umat<-matrix(0,m+1,n+1)
for (i in 0:m){
for (j in 0:n){
if(abs(i/m-j/n)<=val)cmat[i+1,j+1]<-1e0
k<-i+j
if(k > 0 && k<length(z) && z[k]==z[k+1])cmat[i+1,j+1]<-1
}
}
for (i in 0:m){
for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
}
term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
kstiesig<-1.-umat[m+1,n+1]/exp(term)
kstiesig
}



kswsig<-function(m,n,val){
#
#    Compute significance level of the weighted
#    Kolmogorov-Smirnov test statistic
#
#    m=sample size of first group
#    n=sample size of second group
#    val=observed value of test statistic
#
  mpn<-m+n
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 1:m-1){
    for (j in 1:n-1)cmat[i+1,j+1]<-abs(i/m-j/n)*sqrt(m*n/((i+j)*(1-(i+j)/mpn)))
  }
  cmat<-ifelse(cmat<=val,1,0)
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kswsig<-1.-umat[m+1,n+1]/exp(term)
  kswsig
}

Mks.test.handler <- function(w)
  {
    #suppress the following warning:
    #In ks.test() :
    # p-values will be approximate in the presence of ties
    if( any( grepl( "ties", w) ) )
      invokeRestart( "muffleWarning" )

    #invoke as:
    #withCallingHandlers( ks.test(round(x1), round(x2)), warning = Mks.test.handler )
  }

Mks.test <- function(...)
  {
    RVAL <- withCallingHandlers( ks.test(...), warning = Mks.test.handler )
    return(RVAL)
  }

Mt.test  <- function(Tr, Co, weights)
  {
    v1  <- Tr-Co
    estimate  <- sum(v1*weights)/sum(weights)
    var1  <- sum( ((v1-estimate)^2)*weights )/( sum(weights)*sum(weights) )

    parameter  <- Inf

    #get rid of NA for t.test!
    if (estimate==0 & var1==0)
      {
        statistic = 0        
        p.value = 1
      }  else {
        statistic  <- estimate/sqrt(var1)
        
        p.value    <- (1-MATCHpt(abs(statistic), df=sum(weights)-1))*2
      }

    z  <- list(statistic=statistic, parameter=parameter, p.value=p.value,
               estimate=estimate)
    return(z)
  } #end of Mt.test

Mt.test.unpaired  <- function(Tr, Co,
                              weights.Tr=rep(1,length(Tr)),
                              weights.Co=rep(1,length(Co)))
  {
    obs.Tr <- sum(weights.Tr)
    obs.Co <- sum(weights.Co)
    
    mean.Tr <- sum(Tr*weights.Tr)/obs.Tr
    mean.Co <- sum(Co*weights.Co)/obs.Co
    estimate <- mean.Tr-mean.Co
    var.Tr  <- sum( ( (Tr - mean.Tr)^2 )*weights.Tr)/(obs.Tr-1)
    var.Co  <- sum( ( (Co - mean.Co)^2 )*weights.Co)/(obs.Co-1)
    dim <- sqrt(var.Tr/obs.Tr + var.Co/obs.Co)

    parameter  <- Inf

    #get rid of NA for t.test!
    if (estimate==0 & dim==0)
      {
        statistic = 0
        p.value = 1
      }  else {
        statistic  <- estimate/dim

        a1 <- var.Tr/obs.Tr
        a2 <- var.Co/obs.Co
        dof <- ((a1 + a2)^2)/( (a1^2)/(obs.Tr - 1) + (a2^2)/(obs.Co - 1) )    
        
        p.value    <- (1-MATCHpt(abs(statistic), df=dof))*2
      }

    z  <- list(statistic=statistic, parameter=parameter, p.value=p.value,
               estimate=estimate)
    return(z)
  } #end of Mt.test.unpaired

MatchBalance <- function(formul, data=NULL, match.out=NULL, ks=TRUE, 
                         nboots=500, weights=NULL,
                         digits=5, paired=TRUE, print.level=1)
  {
    if(!is.list(match.out) & !is.null(match.out)) {
      warning("'Match' object contains no valid matches")
      match.out  <- NULL
    }

    if ( (class(match.out) != "Match") & (class(match.out) != "Matchby") & (!is.null(match.out)) ) {
      warning("Object not of class 'Match'")
      match.out  <- NULL
    }

    orig.na.action <- as.character(options("na.action"))
    options("na.action"=na.pass)            
    if (is.null(data))
      {
        xdata <- as.data.frame(get.xdata(formul,datafr=environment(formul)))
        Tr <- as.double(get.ydata(formul,datafr=environment(formul)))

      } else {
        data  <- as.data.frame(data)

        xdata  <- as.data.frame(get.xdata(formul, data))
        Tr  <- as.double(get.ydata(formul, data))
      }
    options("na.action"=orig.na.action)

    if (is.null(weights))
      weights <- rep(1,length(Tr))

    if(!is.numeric(weights))
      stop("'weights' must be a numeric vector")

    if( sum(is.na(xdata))!=0 | sum(is.na(Tr))!=0)
      {

        if(orig.na.action!="na.omit" & orig.na.action!="na.exclude" & orig.na.action!="na.fail")
          warning("'na.action' should be set to 'na.omit', 'na.exclude' or 'na.fail' see 'help(na.action)'")

        if (orig.na.action=="na.fail")
          {
            stop("NA's found in data input.")            
          } else {
            warning("NA's found in data input.  IT IS HIGHLY RECOMMENDED THAT YOU TEST IF THE MISSING VALUES ARE BALANCED INSTEAD OF JUST DELETING THEM.")
          }

        indx1 <- apply(is.na(xdata),1,sum)==0 & is.na(Tr)==0
        Tr <- Tr[indx1]
        xdata = xdata[indx1,]
        weights <- weights[indx1]
      } #end of na

    if (sum(Tr !=1 & Tr !=0) > 0) {
      stop("Treatment indicator must be a logical variable---i.e., TRUE (1) or FALSE (0)")
    }

    nvars  <- ncol(xdata)
    names.xdata  <- names(xdata)

    findx  <- 1
    if (sum(xdata[,1]==rep(1,nrow(xdata)))==nrow(xdata))
      {
        findx  <- 2
      }

    if(nboots < 10 & nboots > 0)
      nboots <- 10
    
    if (ks)
      {
        ks.bm <- KSbootBalanceSummary(index.treated=(Tr==0),
                                      index.control=(Tr==1),
                                      X=xdata[,findx:nvars],
                                      nboots=nboots)

        if (!is.null(match.out))
          {
            ks.am <- KSbootBalanceSummary(index.treated=match.out$index.treated,
                                          index.control=match.out$index.control,
                                          X=xdata[,findx:nvars],
                                          nboots=nboots)
          }
      } 

    BeforeMatchingBalance <- list()
    AfterMatchingBalance <- list()

    BMsmallest.p.value <- 1
    BMsmallest.number <- 1
    BMsmallest.name <- names.xdata[findx]

    AMsmallest.p.value <- NULL
    AMsmallest.number <- NULL
    AMsmallest.name <- NULL
    
    if (!is.null(match.out))
      {
        AMsmallest.p.value <- 1
        AMsmallest.number <- 1
        AMsmallest.name <- names.xdata[findx]        
      }

    for( i in findx:nvars)
      {
        count <- i-findx+1
        if(print.level > 0)
          cat("\n***** (V",count,") ", names.xdata[i]," *****\n",sep="")
        
        ks.do  <- FALSE
        is.dummy  <- length(unique( xdata[,i] )) < 3
        if (ks & !is.dummy)
          ks.do  <- TRUE

        BeforeMatchingBalance[[count]]  <-  balanceUV(xdata[,i][Tr==1], xdata[,i][Tr==0], nboots=0,
                                                      weights.Tr=weights[Tr==1], weights.Co=weights[Tr==0],
                                                      match=FALSE)
        
        if (BeforeMatchingBalance[[count]]$tt$p.value < BMsmallest.p.value)
          {
            BMsmallest.p.value <- BeforeMatchingBalance[[count]]$tt$p.value
            BMsmallest.number <- count
            BMsmallest.name <- names.xdata[i]            
          } else if (BeforeMatchingBalance[[count]]$tt$p.value == BMsmallest.p.value)
            {
              BMsmallest.number <- c(BMsmallest.number,count)
              BMsmallest.name <- c(BMsmallest.name,names.xdata[i])
            }
        
        if (ks.do)
          {
            BeforeMatchingBalance[[count]]$ks <- list()
            BeforeMatchingBalance[[count]]$ks$ks <- list()
            BeforeMatchingBalance[[count]]$ks$ks$p.value <- ks.bm$ks.naive.pval[count]
            BeforeMatchingBalance[[count]]$ks$ks$statistic <- ks.bm$ks.stat[count]              
            if (nboots > 0)
              {
                BeforeMatchingBalance[[count]]$ks$ks.boot.pvalue <- ks.bm$ks.boot.pval[count]

                if (ks.bm$ks.boot.pval[count] < BMsmallest.p.value)
                  {
                    BMsmallest.p.value <- ks.bm$ks.boot.pval[count]
                    BMsmallest.number <- count
                    BMsmallest.name <- names.xdata[i]            
                  } else if ( (ks.bm$ks.boot.pval[count] == BMsmallest.p.value) & (sum(BMsmallest.number==count)==0) )
                    {
                      BMsmallest.number <- c(BMsmallest.number,count)
                      BMsmallest.name <- c(BMsmallest.name,names.xdata[i])
                    }
              } else {
                BeforeMatchingBalance[[count]]$ks$ks.boot.pvalue <- NA

                if (ks.bm$ks.naive.pval[count] < BMsmallest.p.value)
                  {
                    BMsmallest.p.value <- ks.bm$ks.naive.pval[count]
                    BMsmallest.number <- count
                    BMsmallest.name <- names.xdata[i]            
                  } else if ( (ks.bm$ks.naive.pval[count] == BMsmallest.p.value) & (sum(BMsmallest.number==count)==0) )
                    {
                      BMsmallest.number <- c(BMsmallest.number,count)
                      BMsmallest.name <- c(BMsmallest.name,names.xdata[i])
                    }              
              }
            
          } else {
            BeforeMatchingBalance[[count]]$ks <- NULL
          }
        
        if (!is.null(match.out))
          {
            AfterMatchingBalance[[count]]  <- balanceUV(xdata[,i][match.out$index.treated],
                                                        xdata[,i][match.out$index.control],
                                                        weights=match.out$weights, nboots=0,
                                                        paired=paired, match=TRUE)
            
            if (AfterMatchingBalance[[count]]$tt$p.value < AMsmallest.p.value)
              {
                AMsmallest.p.value <- AfterMatchingBalance[[count]]$tt$p.value
                AMsmallest.number <- count
                AMsmallest.name <- names.xdata[i]            
              } else if ( (AfterMatchingBalance[[count]]$tt$p.value == AMsmallest.p.value) & (sum(AMsmallest.number==count)==0) )
                    {
                      AMsmallest.number <- c(AMsmallest.number,count)
                      AMsmallest.name <- c(AMsmallest.name,names.xdata[i])
                    }
            
            if (ks.do)
              {                
                AfterMatchingBalance[[count]]$ks <- list()
                AfterMatchingBalance[[count]]$ks$ks <- list()
                AfterMatchingBalance[[count]]$ks$ks$p.value <- ks.am$ks.naive.pval[count]
                AfterMatchingBalance[[count]]$ks$ks$statistic <- ks.am$ks.stat[count]
                if (nboots > 0)
                  {
                    AfterMatchingBalance[[count]]$ks$ks.boot.pvalue <- ks.am$ks.boot.pval[count]

                    if (ks.am$ks.boot.pval[count] < AMsmallest.p.value)
                      {
                        AMsmallest.p.value <- ks.am$ks.boot.pval[count]
                        AMsmallest.number <- count
                        AMsmallest.name <- names.xdata[i]            
                      } else if ( (ks.am$ks.boot.pval[count] == AMsmallest.p.value) & (sum(AMsmallest.number==count)==0) )
                        {
                          AMsmallest.number <- c(AMsmallest.number,count)
                          AMsmallest.name <- c(AMsmallest.name,names.xdata[i])
                        }
                  } else {
                    AfterMatchingBalance[[count]]$ks$ks.boot.pvalue <- NA

                    if (ks.am$ks.naive.pval[count] < AMsmallest.p.value)
                      {
                        AMsmallest.p.value <- ks.am$ks.naive.pval[count]
                        AMsmallest.number <- count
                        AMsmallest.name <- names.xdata[i]            
                      } else if ( (ks.am$ks.naive.pval[count] == AMsmallest.p.value) & (sum(AMsmallest.number==count)==0) )
                        {
                          AMsmallest.number <- c(AMsmallest.number,count)
                          AMsmallest.name <- c(AMsmallest.name,names.xdata[i])
                        }              
                  }                    
              } else {
                AfterMatchingBalance[[count]]$ks <- NULL
              }
            if(print.level > 0)
              PrintBalanceUV(BeforeMatchingBalance[[count]], AfterMatchingBalance[[count]], digits=digits)
          } else {
            if(print.level > 0)
              {
                cat("before matching:\n")
                summary(BeforeMatchingBalance[[count]], digits=digits)
              }
          } #end of if match.out
      } #end of for loop
    
    if(print.level & ( (nvars-findx+1) > 1))
      {
        cat("\n")

        if (BMsmallest.p.value < 1)
          {
            cat("Before Matching Minimum p.value:", format.pval(BMsmallest.p.value, digits=digits),"\n")
            cat("Variable Name(s):",BMsmallest.name, " Number(s):",BMsmallest.number,"\n\n")
          } else {
            cat("Before Matching Minimum p.value: 1\n\n")
          }

        if (!is.null(match.out))
          {
            if(AMsmallest.p.value < 1)
              {
                cat("After Matching Minimum p.value:", format.pval(AMsmallest.p.value, digits=digits),"\n")
                cat("Variable Name(s):",AMsmallest.name, " Number(s):",AMsmallest.number,"\n\n")
              } else {
                cat("After Matching Minimum p.value: 1\n\n")
              }
          } #end of !is.null(match.out)
      }#end of print.level & (nvars > 1)
  
    return(invisible(list(BeforeMatching=BeforeMatchingBalance,
                          AfterMatching=AfterMatchingBalance,
                          BMsmallest.p.value=BMsmallest.p.value,
                          BMsmallestVarName=BMsmallest.name,
                          BMsmallestVarNumber=BMsmallest.number,
                          AMsmallest.p.value=AMsmallest.p.value,
                          AMsmallestVarName=AMsmallest.name,
                          AMsmallestVarNumber=AMsmallest.number)))
  } #end of MatchBalance


get.xdata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "term.labels"))==0 & attr(t1, "intercept")==0) {
    m <- NULL;  # no regressors specified for the model matrix
  } else {
    m <- model.matrix(formul, data=datafr, drop.unused.levels = TRUE)
  }
  return(m);
}


# get.ydata:
# Return response vector corresponding to the formula in formul
# 
get.ydata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "response"))==0) {
    m <- NULL;  # no response variable specified
  }  else {
    m <- model.response(model.frame(formul, data=datafr))
  }
  return(m);
}

#
# bootstrap ks test implemented.  Fast version
#

ks.boot  <- function(Tr, Co, nboots=1000, alternative = c("two.sided", "less", "greater"), print.level=0)
  {
    alternative <- match.arg(alternative)
    tol <- sqrt(.Machine$double.eps)
    Tr <- Tr[!is.na(Tr)]
    Co <- Co[!is.na(Co)]

    w    <- c(Tr, Co)
    obs  <- length(w)
    n.x <- length(Tr)
    n.y <- length(Co)
    cutp <- n.x
    ks.boot.pval <- NULL
    bbcount <- 0

    if (nboots < 10)
      {
        nboots  <- 10
        warning("At least 10 'nboots' must be run; seting 'nboots' to 10")
      }

    if (nboots < 500)
      warning("For publication quality p-values it is recommended that 'nboots'\n be set equal to at least 500 (preferably 1000)") 
    
    fs.ks  <- Mks.test(Tr, Co, alternative=alternative)    

    if (alternative=="two.sided")
      {
        if (print.level > 0)
          cat("ks.boot: two.sided test\n")
        for (bb in 1:nboots)
          {
            if (print.level > 1)
              cat("s:", bb, "\n")
            
            sindx  <- sample(1:obs, obs, replace=TRUE)
            
            X1tmp <- w[sindx[1:cutp]]
            X2tmp <- w[sindx[(cutp+1):obs]]
            
            s.ks <- ks.fast(X1tmp, X2tmp, n.x=n.x, n.y=n.y, n=obs)
            
            if (s.ks >= (fs.ks$statistic - tol) )
              bbcount  <- bbcount + 1
          }
      } else {
        if (print.level > 0)
          cat("ks.boot:",alternative,"test\n")            
        for (bb in 1:nboots)
          {
            if (print.level > 1)
              cat("s:", bb, "\n")
            
            sindx  <- sample(1:obs, obs, replace=TRUE)
            
            X1tmp <- w[sindx[1:cutp]]
            X2tmp <- w[sindx[(cutp+1):obs]]
            
            s.ks <- Mks.test(X1tmp, X2tmp, alternative=alternative)$statistic
            
            if (s.ks >= (fs.ks$statistic - tol) )
              bbcount  <- bbcount + 1                                
          }        
      }
    ks.boot.pval  <- bbcount/nboots
    
    ret  <- list(ks.boot.pvalue=ks.boot.pval, ks=fs.ks, nboots=nboots)
    class(ret)  <- "ks.boot"
    
    return(ret)
  } #end of ks.boot

summary.ks.boot <- function(object, ..., digits=5)
  {
    if(!is.list(object)) {
      warning("object not a valid 'ks.boot' object")
      return()
    }

    if (class(object) != "ks.boot") {
      warning("Object not of class 'ks.boot'")
      return()
    }    
        
    cat("\n")
    cat("Bootstrap p-value:    ", format.pval(object$ks.boot.pvalue, digits=digits), "\n")
    cat("Naive p-value:        ", format(object$ks$p.value, digits=digits), "\n")
    cat("Full Sample Statistic:", format(object$ks$statistic, digits=digits), "\n")
#    cat("nboots completed      ", object$nboots, "\n")
    cat("\n")

    z <- list()
    class(z) <- "summary.ks.boot"
    return(invisible(z))
  } #end of summary.ks.boot

print.summary.ks.boot <- function(x, ...)
  {
    invisible(NULL)
  }


RmatchLoop <- function(Y, Tr, X, Z, V, All, M, BiasAdj, Weight, Weight.matrix, Var.calc, weight,
                       SAMPLE, ccc, cdd, ecaliper=NULL, exact=NULL, caliper=NULL, restrict=NULL,
                       MatchLoopC.indx=NULL, weights.flag, replace=TRUE, ties=TRUE,
                       version="standard", MatchbyAI=FALSE)
  {
    s1 <- MatchGenoudStage1caliper(Tr=Tr, X=X, All=All, M=M, weights=weight,
                                   exact=exact, caliper=caliper,
                                   distance.tolerance=cdd,
                                   tolerance=ccc)

    sum.caliper.drops <- 0
    X.orig <- X

#are we using the restriction matrix?
    if(is.matrix(restrict)) {
      restrict.trigger <- TRUE
    } else {
      restrict.trigger <- FALSE
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

# two slippage parameters that are used to determine whether distances are equal
# distances less than ccc or cdd are interpeted as zero.
# these are passed in.  ccc, cdd


# I. set up
# I.a. vector for which observations we want the average effect
# iot_t is the vector with weights in the average treatment effects
# iot_c is the vector of indicators for potential controls

    if (All==1)
      {
        iot.t <- weight;
        iot.c <- as.matrix(rep(1,length(Tr)))
      } else {
        iot.t <- Tr*weight;
        iot.c <- 1-Tr    
      }

# I.b. determine sample and covariate vector sizes
    N  <- nrow(X)
    Kx <- ncol(X)
    Kz <- ncol(Z)

# K covariates
# N observations
#    Nt <- sum(Tr)
#    Nc <- sum(1-Tr)
#    on <- as.matrix(rep(1,N))

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
# The matrix AA enables one to transform the user supplied weight matrix 
# to take account of this transformation.  BUT THIS IS NOT USED!!
# Mu_X and Sig_X keep track of the original mean and variances
#    AA    <- diag(Kx)
    Mu.X  <- matrix(0, Kx, 1)
    Sig.X <- matrix(0, Kx, 1)

    for (k in 1:Kx)
      {
        Mu.X[k,1] <- sum(X[,k]*weight)/sum(weight)
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weight)/sum(weight)-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        if(Sig.X[k,1] < ccc)
          Sig.X[k,1] <- ccc        
        X[,k]=eps/Sig.X[k,1]
#        AA[k,k]=Sig.X[k,1]
      } #end of k loop

#    Nv <- nrow(V)
    Mv <- ncol(V)
    Mu.V  <- matrix(0, Mv, 1)
    Sig.V <- matrix(0, Mv, 1)

    for (j in 1:Mv)
      {
        Mu.V[j,1]= ( t(V[,j])%*%weight ) /sum(weight)
#        dv <- V[,j]-Mu.V[j,1]
        sv <- sum(V[,j]*V[,j]*weight)/sum(weight)-Mu.V[j,1]^2
        if (sv > 0)
          {
            sv <- sqrt(sv)
          } else {
            sv <- 0
          }
        sv <- sv * sqrt(N/(N-1))
        Sig.V[j,1] <- sv
      } #end of j loop

# I.d. define weight matrix for metric, taking into account normalization of
# regressors.
# If the eigenvalues of the regressors are too close to zero the Mahalanobis metric
# is not used and we revert back to the default inverse of variances.
    if (Weight==1)
      {
        Weight.matrix=diag(Kx)
      } else if (Weight==2) {
        if (min (eigen( t(X)%*%X/N, only.values=TRUE)$values) > 0.0000001)
          {
            Weight.matrix= solve(t(X)%*%X/N)
          } else {
            Weight.matrix <- diag(Kx)
          }
      }
      # DO NOT RESCALE THE Weight.matrix!!
      #else if (Weight==3)
      #  {
      #    Weight.matrix <- AA %*% Weight.matrix %*% AA
      #  }

#    if (exact==1)
#      {
#        Weight.matrix <- cbind(Weight.matrix, matrix(0,nrow=Kx,ncol=Mv))
#        Weight.matrix <- rbind(Weight.matrix, cbind(matrix(0,nrow=Mv,ncol=Kx),
#                               1000*solve(diag(as.vector(Sig.V*Sig.V), nrow=length(Sig.V)))))
#        Weight.matrix <- as.matrix(Weight.matrix)
#        X <- cbind(X,V)
#        Mu.X  <- rbind(Mu.X, matrix(0, nrow(Mu.V), 1))
#        Sig.X <- rbind(Sig.X, matrix(1, nrow(Sig.V), 1))
#      } #end of exact

    if ( min(eigen(Weight.matrix, only.values=TRUE)$values) < ccc )
      Weight.matrix <- Weight.matrix + diag(Kx)*ccc

    ww <- chol(Weight.matrix) # so that ww*ww=w.m

    if(is.null(s1$ecaliper))
      {
        caliperflag <- 0
        use.ecaliper <- 0
      } else {
        caliperflag <- 1
        use.ecaliper <- s1$ecaliper
      }

    #if we have a diagonal matrix we can void cblas_dgemm
    if (Kx > 1)
      {
        DiagWeightMatrixFlag <- as.double(sum( (Weight.matrix!=diag(diag(Weight.matrix))) )==0)
      } else {
        DiagWeightMatrixFlag <- 1
      }

    if(is.null(MatchLoopC.indx))
      {
    #indx:
    # 1] I (unadjusted); 2] IM (unadjusted); 3] weight; 4] I (adjusted); 5] IM (adjusted)
        if(weights.flag==TRUE)
          {
            MatchLoopC.indx <- MatchLoopC(N=s1$N, xvars=Kx, All=s1$All, M=s1$M,
                                          cdd=cdd, caliperflag=caliperflag, replace=replace, ties=ties,
                                          ww=ww, Tr=s1$Tr, Xmod=s1$X,
                                          weights=weight,
                                          CaliperVec=use.ecaliper, Xorig=X.orig,
                                          restrict.trigger=restrict.trigger, restrict=restrict,
                                          DiagWeightMatrixFlag=DiagWeightMatrixFlag)
          } else {
            MatchLoopC.indx <- MatchLoopCfast(N=s1$N, xvars=Kx, All=s1$All, M=s1$M,
                                              cdd=cdd, caliperflag=caliperflag, replace=replace, ties=ties,
                                              ww=ww, Tr=s1$Tr, Xmod=s1$X,
                                              CaliperVec=use.ecaliper, Xorig=X.orig,
                                              restrict.trigger=restrict.trigger, restrict=restrict,
                                              DiagWeightMatrixFlag=DiagWeightMatrixFlag)
          }
      }

    indx <- MatchLoopC.indx

    if(indx[1,1]==0)
      {
        ret <- list()
        ret$valid <- 0
        if (caliperflag)
          {
            ret$sum.caliper.drops <- indx[1,6]
          } else {
            ret$sum.caliper.drops <- 0
          }
        return(ret)
      } 
    #we now keep going if we only have 1 valid match
    #else if (nrow(indx)< 2)
    #  {
    #    ret <- list()
    #    ret$valid <- 1
    #    if (caliperflag)
    #      {
    #        ret$sum.caliper.drops <- indx[1,6]
    #      } else {
    #        ret$sum.caliper.drops <- 0
    #      }
    #    return(ret)
    #  }
        
    if (All==2)
      {
        foo <- indx[,5]
        indx[,5] <- indx[,4]
        indx[,4] <- foo
      }

    if (caliperflag)
      {
        sum.caliper.drops <- indx[1,6]
      } else {
        sum.caliper.drops <- 0
      }

    #
    # Generate variables which we need later on
    #

    I <- indx[,1]
    IT <- Tr[indx[,1]]
    IM <- indx[,2]

#    IX <- X[indx[,1],]
#    Xt <- X[indx[,4],]
#    Xc <- X[indx[,5],]

#    IY <- Y[indx[,1]]
    Yt <- Y[indx[,4]]
    Yc <- Y[indx[,5]]

    W <- indx[,3]
    
    if(BiasAdj==1 & sum(W) < ncol(Z))
      {
        warning("Fewer (weighted) matches than variables in 'Z': BiasAdjust set to FALSE")
        BiasAdj=0
      }    

    if(BiasAdj==1)
      {
        if(sum(W) < ncol(Z))
          {
            warning("Fewer matches than variables for Bias Adjustment")
          }
        
        IZ <- Z[indx[,1],]
        Zt <- Z[indx[,4],]
        Zc <- Z[indx[,5],]
      }

    est.func <- function(N, All, Tr, indx, weight, BiasAdj, Kz)
      {
        Kcount <- as.matrix(rep(0, N))    
        KKcount <- as.matrix(rep(0, N))
        YCAUS <- matrix(0, nrow=N, ncol=1)

        if (BiasAdj==1)
          {
            ZCAUS <- matrix(0, nrow=N, ncol=Kz)
          } else {
            ZCAUS <- NULL
          }
        
        for (i in 1:N)
          {
            if ( ( Tr[i]==1 & All!=1) | All==1 )
              {

                foo.indx <- indx[,1]==i
                foo.indx2 <- foo.indx
                
                sum.foo <- sum(foo.indx)
                
                if (sum.foo < 1)
                  next;            
                
                foo <- rep(FALSE, N)
                foo.indx <- indx[foo.indx,2]
                foo[foo.indx] <- rep(TRUE,sum.foo)

                #   inner.func <- function(N, weight, indx, foo.indx2, Y, Tr, foo)
  
                Kcount <- Kcount + weight[i] * weight*foo/sum(foo*weight)
                
                KKcount <- KKcount + weight[i]*weight*weight*foo/
                  (sum(foo*weight)*sum(foo*weight))
                
                foo.indx2.2 <- indx[foo.indx2,2];
                foo.indx2.3 <- indx[foo.indx2,3];
                sum.foo.indx2.3 <- sum(foo.indx2.3)

                if(Tr[i]==1)
                  {
                    YCAUS[i] <- Y[i] - sum((Y[foo.indx2.2]*foo.indx2.3))/sum.foo.indx2.3
                  } else {
                    YCAUS[i] <- sum((Y[foo.indx2.2]*foo.indx2.3))/sum.foo.indx2.3 - Y[i]
                  }

                if (BiasAdj==1)
                  {
                    
                    if(Tr[i]==1)
                      {
                        if (sum.foo > 1)
                          {
                            ZCAUS[i,] <-
                              Z[i,] - t(Z[foo.indx2.2,]) %*% foo.indx2.3/sum.foo.indx2.3
                          } else {
                            ZCAUS[i,] <-
                              Z[i,] - Z[foo.indx2.2,]*foo.indx2.3/sum.foo.indx2.3
                          }
                      } else {
                        if (sum.foo > 1)
                          {
                            ZCAUS[i,] <- t(Z[foo.indx2.2,]) %*% foo.indx2.3/sum.foo.indx2.3 - Z[i,]
                          } else {
                            ZCAUS[i,] <- Z[foo.indx2.2,]*foo.indx2.3/sum.foo.indx2.3 - Z[i,]
                          }
                      }
                  } #endof BiasAdj
              }
          } #end of if
        return(list(YCAUS=YCAUS,ZCAUS=ZCAUS,Kcount=Kcount,KKcount=KKcount))
      } #end of est.func
    
    if(version=="standard" & BiasAdj==0)
      {
        ret <- .Call("EstFuncC", as.integer(N), as.integer(All), as.integer(nrow(indx)),
                     as.double(Y), as.double(Tr),
                     as.double(weight), as.double(indx),
                     PACKAGE="Matching")
        YCAUS <- ret[,1];
        Kcount <- ret[,2];
        KKcount <- ret[,3];
      } else if (version=="standard") {
        ret.est <- est.func(N=N, All=All, Tr=Tr, indx=indx, weight=weight, BiasAdj=BiasAdj, Kz=Kz)
        YCAUS <- ret.est$YCAUS
        ZCAUS <- ret.est$ZCAUS
        Kcount <- ret.est$Kcount
        KKcount <- ret.est$KKcount
      }

    if (All!=1)
      {
        I  <- as.matrix(I[IT==1])
        IT <- as.matrix(IT[IT==1])
        Yc <- as.matrix(Yc[IT==1])
        Yt <- as.matrix(Yt[IT==1])
        W  <- as.matrix(W[IT==1])
        if (BiasAdj==1)
          {
            if (Kz > 1)
              {
                Zc <- as.matrix(Zc[IT==1,])
                Zt <- as.matrix(Zt[IT==1,])
                IZ <- as.matrix(IZ[IT==1,])
              } else{
                Zc <- as.matrix(Zc[IT==1])
                Zt <- as.matrix(Zt[IT==1])
                IZ <- as.matrix(IZ[IT==1])            
              }
          }
#        IM <- as.matrix(IM[IT==1,])
#        IY <- as.matrix(IY[IT==1])
#        IX.u  <- as.matrix(IX.u[IT==1,])
#        Xc.u  <- as.matrix(Xc.u[IT==1,])
#        Xt.u  <- as.matrix(Xt.u[IT==1,])
#        Xc    <- as.matrix(Xc[IT==1,])
#        Xt <- as.matrix(Xt[IT==1,])
      } #end of if
    
    if (length(I) < 1)
      {
        return(list(sum.caliper.drops=N))
      }
    
    if (BiasAdj==1)
      {
        # III. Regression of outcome on covariates for matches
        if (All==1)
          {
            # number of observations
            NNt <- nrow(Z)
            # add intercept        
            ZZt <- cbind(rep(1, NNt), Z)
            # number of covariates
            Kx <- ncol(ZZt)
            xw <- ZZt*(sqrt(Tr*Kcount) %*% t(as.matrix(rep(1,Kx))))
            
            foo <- min(eigen(t(xw)%*%xw, only.values=TRUE)$values)
            foo <- as.double(foo<=ccc)
            foo2 <- apply(xw, 2, sd)

            options(show.error.messages = FALSE)
            wout <- NULL
            try(wout <- solve( t(xw) %*% xw + diag(Kx) * ccc * (foo) * max(foo2)) %*%
                (t(xw) %*% (Y*sqrt(Tr*Kcount))))
            if(is.null(wout))
              {
                wout2 <- NULL
                try(wout2 <- ginv( t(xw) %*% xw + diag(Kx) * ccc * (foo) * max(foo2)) %*%
                    (t(xw) %*% (Y*sqrt(Tr*Kcount))))
                if(!is.null(wout2))
                  {
                    wout <-wout2
                    warning("using generalized inverse to calculate Bias Adjustment probably because of singular 'Z'")
                  }
              }
            options(show.error.messages = TRUE)
            if(is.null(wout))
              {
                warning("unable to calculate Bias Adjustment probably because of singular 'Z'")
                BiasAdj <- 0
              } else {
                NW <- nrow(wout)
#                KW <- ncol(wout)
                Alphat <- wout[2:NW,1]
              }
          } else {
            Alphat <- matrix(0, nrow=Kz, ncol=1)
          } #end if ALL
      }

    if(BiasAdj==1)
      {
        # III.b.  Controls
        NNc <- nrow(Z)
        ZZc <- cbind(matrix(1, nrow=NNc, ncol=1),Z)
        Kx <- ncol(ZZc)
        
        xw <- ZZc*(sqrt((1-Tr)*Kcount) %*% matrix(1, nrow=1, ncol=Kx))
        
        foo <- min(eigen(t(xw)%*%xw, only.values=TRUE)$values)
        foo <- as.double(foo<=ccc)
        foo2 <- apply(xw, 2, sd)

        options(show.error.messages = FALSE)
        wout <- NULL        
        try(wout <- solve( t(xw) %*% xw + diag(Kx) * ccc * (foo) * max(foo2)) %*%
            (t(xw) %*% (Y*sqrt((1-Tr)*Kcount))))
        if(is.null(wout))
          {
            wout2 <- NULL
            try(wout2 <- ginv( t(xw) %*% xw + diag(Kx) * ccc * (foo) * max(foo2)) %*%
                (t(xw) %*% (Y*sqrt((1-Tr)*Kcount))))
            if(!is.null(wout2))
              {
                wout <-wout2
                warning("using generalized inverse to calculate Bias Adjustment probably because of singular 'Z'")
              }
          }        
        options(show.error.messages = TRUE)
        if(is.null(wout))
          {
            warning("unable to calculate Bias Adjustment probably because of singular 'Z'")
            BiasAdj <- 0
          } else {
            NW <- nrow(wout)
#            KW <- ncol(wout)
            Alphac <- as.matrix(wout[2:NW,1])
          }
      }

    if(BiasAdj==1)
      {
        # III.c. adjust matched outcomes using regression adjustment for bias adjusted matching estimator
        IZ <- as.matrix(IZ)
        Zc <- as.matrix(Zc)
        Zt <- as.matrix(Zt)
        Alphat <- as.matrix(Alphat)

        SCAUS <- YCAUS-Tr*(ZCAUS %*% Alphac)-(1-Tr)*(ZCAUS %*% Alphat)
        # adjusted control outcome
        Yc.adj <- Yc+BiasAdj * (IZ-Zc) %*% Alphac
        # adjusted treated outcome
        Yt.adj <- Yt+BiasAdj*(IZ-Zt) %*% Alphat
        Tau.i <- Yt.adj - Yc.adj
      } else {
        Yc.adj <- Yc
        Yt.adj <- Yt
        Tau.i <- Yt.adj - Yc.adj        
      }
       
    art.data <- cbind(I,IM)

    # III. If conditional variance is needed, initialize variance vector
    # and loop through all observations

    if (Var.calc>0)
      {
        #For R version of this function see Matching version < 4.5-0008
        Sigs <- VarCalcMatchC(N=N, xvars=ncol(X), Var.calc=Var.calc,
                              cdd=cdd, caliperflag=caliperflag, 
                              ww=ww, Tr=Tr, Xmod=s1$X,
                              CaliperVec=use.ecaliper, Xorig=X.orig,
                              restrict.trigger=restrict.trigger, restrict=restrict,
                              DiagWeightMatrixFlag=DiagWeightMatrixFlag,
                              Y=Y, weightFlag=weights.flag, weight=weight)            
      } #end of var.calc > 0

    est <- t(W) %*% Tau.i/sum(W) # matching estimator

    if(version=="standard")
      {
        if (Var.calc==0)
          {
            eps <- Tau.i - as.double(est)
            eps.sq <- eps*eps
            Sigs <- 0.5 * matrix(1, N, 1) %*% (t(eps.sq) %*% W)/sum(W) #sss <- sqrt(Sigs[1,1])
          } #end of Var.calc==0
        
        SN <- sum(iot.t)
        var.sample <- sum((Sigs*(iot.t+iot.c*Kcount)*(iot.t+iot.c*Kcount))/(SN*SN))
        
        if (All==1)
          {
            var.pop <- sum((Sigs*(iot.c*Kcount*Kcount+2*iot.c*Kcount-iot.c*KKcount))/(SN*SN))
          } else {
            var.pop=sum((Sigs*(iot.c*Kcount*Kcount-iot.c*KKcount))/(SN*SN))
          }
        
        if (BiasAdj==1)
          {
            dvar.pop <- sum(iot.t*(SCAUS-as.double(est))*(SCAUS-as.double(est)))/(SN*SN)
          } else {
            dvar.pop <- sum(iot.t*(YCAUS-as.double(est))*(YCAUS-as.double(est)))/(SN*SN)
          }
        
        var.pop <- var.pop + dvar.pop
        
        if (SAMPLE==1)
          {
            var <- var.sample
          } else {
            #var <- max(var.sample, var.pop)
            var <- var.pop
          }
        
        var.cond <- max(var.sample,var.pop)-var.sample
        se <- sqrt(var)
        se.cond <- sqrt(var.cond)
        #Sig <- sqrt(Sigs)        
      } else {
        se=NULL
        se.cond=NULL
      }
        
    if (All==2)
      est <- -1*est

#    if (exact==1)
#      {
#        Vt.u <- Xt.u[,(Kx-Mv+1):Kx]
#        Vc.u <- Xc.u[,(Kx-Mv+1):Kx]
#        Vdif <- abs(Vt.u-Vc.u)
#
#        if (Mv>1)
#          Vdif <- as.matrix(apply(t(Vdif), 2, sum))
#
#        em[1,1] <- length(Vdif)
#        em[2,1] <- sum(Vdif>0.000001)
#      }#end of exact==1

    if(!MatchbyAI)
      {
        return(list(est=est, se=se, se.cond=se.cond, W=W,
                    sum.caliper.drops=sum.caliper.drops,
                    art.data=art.data, 
                    MatchLoopC=MatchLoopC.indx))
      } else {
        if(Var.calc==0)
          Sigs <- NULL
        return(list(est=est, se=se, se.cond=se.cond, W=W,
                    sum.caliper.drops=sum.caliper.drops,
                    art.data=art.data, 
                    MatchLoopC=MatchLoopC.indx,
                    YCAUS=YCAUS, Kcount=Kcount, KKcount=KKcount,
                    Sigs=Sigs))
      }
  }# end of RmatchLoop

MatchLoopC <- function(N, xvars, All, M, cdd, caliperflag, replace, ties, ww, Tr, Xmod, weights, CaliperVec, Xorig,
                       restrict.trigger, restrict, DiagWeightMatrixFlag)
  {

    if(restrict.trigger)
      {
        restrict.nrow <- nrow(restrict)
      } else {
        restrict.nrow <- 0
      }    

    ret <- .Call("MatchLoopC", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                 as.double(cdd), as.integer(caliperflag), as.integer(replace), as.integer(ties),
                 as.double(ww), as.double(Tr),
                 as.double(Xmod), as.double(weights), as.double(CaliperVec), as.double(Xorig),
                 as.integer(restrict.trigger), as.integer(restrict.nrow), as.double(restrict),
                 as.double(DiagWeightMatrixFlag),
                 PACKAGE="Matching")
    return(ret)
  } #end of MatchLoopC

MatchLoopCfast <- function(N, xvars, All, M, cdd, caliperflag, replace, ties, ww, Tr, Xmod, CaliperVec,
                           Xorig, restrict.trigger, restrict, DiagWeightMatrixFlag)
  {

    if(restrict.trigger)
      {
        restrict.nrow <- nrow(restrict)
      } else {
        restrict.nrow <- 0
      }    
    
    ret <- .Call("MatchLoopCfast", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                 as.double(cdd), as.integer(caliperflag), as.integer(replace), as.integer(ties),
                 as.double(ww), as.double(Tr),
                 as.double(Xmod), as.double(CaliperVec), as.double(Xorig),
                 as.integer(restrict.trigger), as.integer(restrict.nrow), as.double(restrict),
                 as.double(DiagWeightMatrixFlag),
                 PACKAGE="Matching")
    return(ret)
  } #end of MatchLoopCfast

VarCalcMatchC <- function(N, xvars, Var.calc, cdd, caliperflag, ww, Tr, Xmod, CaliperVec,
                          Xorig, restrict.trigger, restrict, DiagWeightMatrixFlag,
                          Y, weightFlag, weight)
  {

    if(restrict.trigger)
      {
        restrict.nrow <- nrow(restrict)
      } else {
        restrict.nrow <- 0
      }    
    
    ret <- .Call("VarCalcMatchC", as.integer(N), as.integer(xvars), as.integer(Var.calc),
                 as.double(cdd), as.integer(caliperflag), 
                 as.double(ww), as.double(Tr),
                 as.double(Xmod), as.double(CaliperVec), as.double(Xorig),
                 as.integer(restrict.trigger), as.integer(restrict.nrow), as.double(restrict),
                 as.double(DiagWeightMatrixFlag),
                 as.double(Y),
                 as.integer(weightFlag), as.double(weight),
                 PACKAGE="Matching")
    return(ret)
  } #end of VarCalcMatchC



.onAttach <- function( ... )
{
  MatchLib <- dirname(system.file(package = "Matching"))
  version <- packageDescription("Matching", lib.loc = MatchLib)$Version
  BuildDate <- packageDescription("Matching", lib.loc = MatchLib)$Date

  foo <- paste("## \n##  Matching (Version ", version, ", Build Date: ", BuildDate, ")\n", 
               "##  See http://sekhon.berkeley.edu/matching for additional documentation.\n",
               "##  Please cite software as:\n",
               "##   Jasjeet S. Sekhon. 2011. ``Multivariate and Propensity Score Matching\n",
               "##   Software with Automated Balance Optimization: The Matching package for R.''\n",
               "##   Journal of Statistical Software, 42(7): 1-52. \n##\n",
               sep = "")
  packageStartupMessage(foo)
}


