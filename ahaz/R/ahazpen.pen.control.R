"ahaz.format.init"<-function(initsol, name)
  {
    ## --------------------------------------------------------------------------
    ## Purpose: Format "initial solutions" (vectors/functions) supplied
    ##          to lasso or sscad
    ## --------------------------------------------------------------------------
    
    if(is.numeric(initsol)){
        
        ##### OUTPUT
        f <-function(surv, X, weights)
          {
            if(any(is.nan(initsol)) || length(initsol) != ncol(X))
              stop(paste("'",name,"' has incorrect dimensions",sep=""))
              return(initsol)
          }
         #####
          
      } else if(is.function(initsol)){
            # Check that we use the correct formal args
          if(!identical(names(formals(initsol)),c("surv","X","weights")))
            stop(paste("arguments of function '", name, "' should be precisely 'surv', 'X' and 'weights'",sep=""))
                 
          ##### OUTPUT
          f<-function(surv, X, weights)
            {
              out<-initsol(surv = surv,X = X,weights = weights)
              if(!is.numeric(out) || any(is.nan(out)) || length(out) != ncol(cbind(X)))
                stop("output of 'init.sol'/'ada.wgt' has incorrect dimensions/class")
              return(out)
            }
           #####
          
        } else {
          stop(paste("'",name,"' has the wrong class",sep=""))
        }
    return(f)
  }

"lasso.control"<-function(alpha = 1, ada.wgt=NULL)
  {
    ## Purpose: Control function for lasso
    ##          to be used in ahazpen calls
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   alpha         : elasticnet penalty parameter, 0<alpha<1. 
    ##   initsol       : initial solution. Tries to default to ahaz estimates
    ##                   otherwise defaults to zero
    ##   nsteps        : How many iterations (1 = one-step SCAD; > 1 ~ full SCAD)
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    if(!is.numeric(alpha) || alpha>1 || alpha<=0 ){ stop("'alpha' incorrectly specified") }

    # Weights for adaptive lasso
    if(!is.null(ada.wgt)){
      ada.wgt<-ahaz.format.init(ada.wgt,"ada.wgt")
    }
    
    out<-list("type"="lasso","init.sol"=NULL,"alpha"=alpha,"nsteps" = 1,"a"=1,"c"=NULL,
              "ada.wgt"=ada.wgt)
    class(out)<-"ahazpen.pen.control"
    return(out)
  }

"sscad.control"<-function(a = 3.7,nsteps = 1,init.sol=NULL, c=NULL)
  {
    ## Purpose: Control function for stepwise SCAD
    ##          to be used in ahazpen calls
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   a             : constant in penalty function. Should be > 2, defaults
    ##                   to 3.7 as per recommendations in literature
    ##   initsol       : initial solution. Tries to default to ahaz estimates
    ##                   otherwise defaults to zero
    ##   nsteps        : How many iterations (1 = one-step SCAD; > 1 ~ full SCAD)
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
     uspec<-FALSE
    if(!is.numeric(nsteps) || nsteps < 1){ stop("'nsteps' must be positive") }
    if(!is.numeric(a) || a < 2){ stop("'a' must be > 2") }
    if(!is.null(c) && (!is.numeric(c) || c < 0) ){ stop("'c' must be > 2") }

    # Initial solution for sscad
    if(is.null(init.sol)){
      init.sol<-function(surv,X,weights)
        {
          X<-cbind(X)
          if(nrow(X)<=ncol(X)){
            return(rep(0,ncol(X)))
          } else {
             beta<-coef(ahaz(surv=surv,X=X,weights=weights))
             if(any(is.na(beta)))
               stop("regression coefficients not defined because of singularities")
             return(beta)
           }
        }
    } else {
      init.sol<-ahaz.format.init(init.sol,"init.sol")
    }
    
    out<-list(type="sscad","init.sol"=init.sol,"alpha"=1,"nsteps" = nsteps,"a"=a,"c"=c,"ada.wgt"=NULL)
    class(out)<-"ahazpen.pen.control"
    return(out)
  }
