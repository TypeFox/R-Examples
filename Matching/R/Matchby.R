Matchby <- function(Y=NULL, Tr, X, by, estimand="ATT", M=1, ties=FALSE, replace=TRUE,
                    exact = NULL, caliper = NULL,
                    AI=FALSE, Var.calc=0,
                    Weight = 1, Weight.matrix = NULL, 
                    distance.tolerance = 1e-05, tolerance = sqrt(.Machine$double.eps), 
                    print.level=1, version="Matchby", ...)
  {

    #index for raw obs
    Tr <- as.double(Tr)    
    nobs  <- length(Tr)
    orig.treated.nobs  <- sum(Tr==1)

    if(is.null(Y))
      {
        Y <- rep(0, nobs)
      } else {
        Y  <- as.double(Y)
        if (nobs != length(Y))
          {
            stop("length(Tr) != length(Y)")
          }
      }
    
    X  <- as.matrix(X)

    if( nobs != nrow(X))
      {
        stop("length(Tr) != nrow(X)")
      }        

    index.nobs <- 1:nobs
    t.index.nobs <- split(index.nobs, by, drop=TRUE)

    nindx  <- length(t.index.nobs)

    # *OUTPUT* observation weights
    weights  <- NULL

    index.treated <- NULL
    index.control <- NULL

    if (Var.calc < 0)
      {
        warning("User set 'Var.calc' to less than 0.  Resetting to the default which is 0.")
        Var.calc <- 0
      }

    if(Var.calc > 0)
      AI <- TRUE

    if(AI & estimand!="ATT")
      {
        warning("This function can only calculate Abadie-Imbens SEs for 'ATT'.  For AI SEs and other estimands, please use 'Match()'. Setting 'AI=FALSE'")
        AI <- FALSE
        Var.calc <- 0
      }
    if(AI & ties!=TRUE)
      {
        ties=TRUE
        warning("Abadie-Imbens SEs have been requested. Setting 'ties=TRUE'")
      }
    if(AI & replace!=TRUE)
      {
        replace=TRUE
        warning("Abadie-Imbens SEs have been requested. Setting 'replace=TRUE'")
      }
    if(AI)
      version <- "MatchbyAI"

    if (replace!=FALSE & replace!=TRUE)
      {
        warning("'replace' must be TRUE or FALSE.  Setting to TRUE")
        replace <- TRUE
      }
    if(replace==FALSE)
      {
        ties <- FALSE        
      }    
    if (ties!=FALSE & ties!=TRUE)
      {
        warning("'ties' must be TRUE or FALSE.  Setting to TRUE")
        ties <- TRUE
      }

    if(AI)
      {
        Kcount <- as.matrix(rep(0, nobs))    
        KKcount <- as.matrix(rep(0, nobs))
        YCAUS <- matrix(0, nrow=nobs, ncol=1)

        if(Var.calc>0)
          Sigs <- matrix(0, nrow=nobs, ncol=1)
      }

    for (i in 1:nindx)
      {
        if(print.level > 0)
          cat(i,"of", nindx, "groups\n")

        tmp.index.nobs <- t.index.nobs[[i]]
        
        f.Tr  <- Tr[tmp.index.nobs]
        #need at least 1 Tr and 1 Co
        if(length(f.Tr) < 2)
          {
            next;
          }
        #drop if we only have Tr or Co
        if( var(f.Tr)==0 )
          {
            next;
          }

        t1  <- Match(Y=Y[tmp.index.nobs], Tr=f.Tr, X=X[tmp.index.nobs,],
                     estimand=estimand, M=M,
                     Var.calc=Var.calc,
                     exact=exact, caliper=caliper, replace=replace, ties=ties,
                     Weight=Weight, Weight.matrix=Weight.matrix,
                     tolerance=tolerance, distance.tolerance=distance.tolerance,
                     version=version, ...)
        
        if(is.na(t1[1]))
          {
            if(!is.null(exact) | !is.null(caliper))
              {
                warning("no matches found in group ",i," (probably because of the exact or caliper option) continuing")
              } else {
                warning("no matches found in group ",i," continuing")}
            next
          }

        weights  <- c(weights, t1$weights)

        index.treated <- c(index.treated,(tmp.index.nobs)[t1$index.treated])
        index.control <- c(index.control,(tmp.index.nobs)[t1$index.control])

        if(AI)
          {
            YCAUS[tmp.index.nobs]   <- t1$YCAUS
            Kcount[tmp.index.nobs]  <- t1$Kcount
            KKcount[tmp.index.nobs] <- t1$KKcount

            if(Var.calc>0)
              {
                if(is.null(t1$Sigs))
                  {
                    pfoo <- paste("Var.calc=",Var.calc," cannot be calculated (group ",i,").  Var.calc is probably set to a number larger than the possible number of matches in a subgroup",sep="")
                    stop(pfoo)
                  } else {
                    Sigs[tmp.index.nobs]    <- t1$Sigs
                  }
              }#if(Var.calc>0)
          }#if(AI)
      }#i:nindx

    if(is.null(index.treated))
      {
        warning("no valid matches were found")
        z <- NA
        class(z)  <- "Matchby"    
        return(z)
      }

    Yt  <- Y[index.treated]
    Yc  <- Y[index.control]
    
    sum.tw  <- sum(weights)
    est  <- sum(Yt*weights)/sum.tw  - sum(Yc*weights)/sum.tw
    Tau  <- Yt - Yc
    varest  <- sum( ((Tau-est)^2)*weights)/(sum.tw*sum.tw)
    se.standard  <- sqrt(varest)
    
    ret  <- list(est=est, se=NULL, se.standard=se.standard,
                 index.treated=index.treated,
                 index.control=index.control,
                 weights=weights)
    ret$orig.nobs  <- nobs
    ret$orig.wnobs <- nobs
    ret$nobs  <- length(Yt)
    ret$wnobs  <- sum.tw
    ret$orig.treated.nobs  <- orig.treated.nobs
    ret$ndrops  <- orig.treated.nobs-ret$wnobs
    ret$estimand  <- estimand
    ret$version  <- version
    class(ret)  <- "Matchby"
    
    if(AI)
      {
        if(Var.calc==0)
          {
            eps <- Tau - est
            eps.sq <- eps*eps
            Sigs <- 0.5 * matrix(1, nobs, 1) %*% (t(eps.sq) %*% weights)/sum(weights)
          }
        
        SN <- orig.treated.nobs
        var.pop=sum((Sigs*((1-Tr)*Kcount*Kcount-(1-Tr)*KKcount))/(SN*SN))
        
        dvar.pop <- sum(Tr*(YCAUS-est)*(YCAUS-est))/(SN*SN)
        
        var <- var.pop + dvar.pop
        
        se <- sqrt(var)
        ret$se <- se
      }
    
    invisible(return(ret))
  }
  

summary.Matchby  <- function(object, ..., digits=5)
  {
    if(!is.list(object)) {
      warning("'Matchby' object contains less than two valid matches.  Cannot proceed.")
      return(invisible(NULL))
    }
    
    if (class(object) != "Matchby") {
      warning("Object not of class 'Matchby'")
      return(invisible(NULL))
    }

    cat("\n")
    cat("Estimate... ",format(object$est,digits=digits),"\n")

    if(!is.null(object$se))
      {
        cat("AI SE...... ",format(object$se,digits=digits),"\n")
        cat("AI T-stat.. ",format(object$est/object$se,digits=digits),"\n")
        cat("AI p.val... ",format.pval((1-pnorm(abs(object$est/object$se)))*2,digits=digits),"\n\n")
      }
    cat("SE......... ",format(object$se.standard,digits=digits),"\n")
    cat("T-stat..... ",format(object$est/object$se.standard,digits=digits),"\n")
    cat("p.val...... ",format.pval((1-pnorm(abs(object$est/object$se.standard)))*2,digits=digits),"\n")
    cat("\n")        

    cat("Original number of observations.............. ", object$orig.nobs,"\n")
    if(object$estimand!="ATC")
      {
        cat("Original number of treated obs............... ", object$orig.treated.nobs,"\n")
      } else {
        cat("Original number of control obs............... ",
            object$orig.nobs- object$orig.treated.nobs,"\n")        
      }
    cat("Matched number of observations............... ", round(object$wnobs, 3),"\n")
    cat("Matched number of observations  (unweighted). ", object$nobs,"\n")
    cat("\n")
        
    if ( object$ndrops > 0)
        {
          cat("Number of treated observations dropped....... ", object$ndrops, "\n")
          cat("\n")
        }

    z <- list()
    class(z) <- "summary.Matchby"
    return(invisible(z))
  } #end of summary.Matchby

print.summary.Matchby <- function(x, ...)
  {
    invisible(NULL)
  }

