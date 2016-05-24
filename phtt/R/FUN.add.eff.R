#########################################################################################
## FUN.add.eff
## Aim: Calculates Additive effects and Intercept(function) of a
## panel-data model with unobserved factor structure.
## The factor structure can be:
## a)empirical:        estimated by fAFactMod() given to the function via "fAFactMod.obj" 
## b)hypothetical:     given to the function via "g.fun"
## 
## Takes:
## 1)PF.obj        (return-object of FUN.Pformula(): Original data and Transformed Data)
## 2)fAFactMod.obj (return-object of fAFactMod())
## 3)beta.hat      ((re-)estimated beta coeficients)
##
## Gives a list with:
## mu     (Overallmean)
## tau    (individal effects)
## beta.0 (intercept-function)
##########################################################################################



FUN.add.eff <- function(PF.obj, fpca.fit.obj=NULL, g.fun=NULL, beta.hat)
  {
    if(!is.null(fpca.fit.obj) & !is.null(g.fun)){
      stop("Only one of the arguments >>fpca.fit.obj<< or >>g.fun<< is allowed to be used.")
    }
    P         <- length(PF.obj)-1
    y.in.list <- PF.obj[[1]]
    T         <- nrow(PF.obj[[1]]$ODM)
    N         <- ncol(PF.obj[[1]]$ODM)

    given.d   <- fpca.fit.obj$given.d

      
    ##=========================================================================================================
    YInC  <- y.in.list$TRm$InC                                 ## *Y**In*dividual *C*onstants
    YTiVC <- y.in.list$TRm$TiVC                                ## *Y**Ti*me V*arying *C*onstants
    YOVc  <- y.in.list$TRm$OVc                                 ## *Y**OV*erall *c*onstant          
    XInC  <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$InC)   ## *X**In*dividual *C*onstants
    XTiVC <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$TiVC)  ## *X**Ti*me V*arying *C*onstants
    XOVc  <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$OVc)   ## *X**OV*erall *c*onstant      
    ##=========================================================================================================

    ## mu:        overall mean effect
    mu           <- ifelse(y.in.list$I, YOVc - XOVc %*% beta.hat, 0)
    
    if(y.in.list$Tr=="individual"|y.in.list$Tr=="twoway"){
      ## tau:       individual effects
      tau        <- c((YInC  - YOVc) - (XInC  - matrix(rep(XOVc, each=N), N, P)) %*% beta.hat) ## dim(tau): Nx1  
    }else{
      ## diese else()-Abfrage könnte weg, da bei "none"/"time" die elemente in YInC gleich wie YOVc sind 
      tau    <- 0
    }
    if(y.in.list$Tr=="time"|y.in.list$Tr=="twoway"){
      ## see section 3.1, paper KSS-2009:
      tmp        <- (YTiVC - YOVc) - (XTiVC - matrix(rep(XOVc, each=T), T, P)) %*% beta.hat     ## dim(tmp): Tx1
      if(is.null(g.fun)){# if empirical factor structure
        ## theta.bar: scores regarding to TiVC
        theta.bar  <-  qr.solve(fpca.fit.obj$factors[,1:given.d, drop= FALSE], tmp)
        ## beta.0:    functional time effects
        beta.0     <-  fpca.fit.obj$factors[,1:given.d, drop= FALSE] %*% theta.bar         ## dim(tmp): Tx1
      }else{# if hypothetical factor structure
        theta.bar  <-  qr.solve(g.fun, tmp)
        beta.0     <-  g.fun %*% theta.bar
      }      
    }else{beta.0 <- 0}
    result    <- list(mu = mu, tau = tau, beta.0 = beta.0)
    return(result)
  }
