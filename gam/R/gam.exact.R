"gam.exact" <-
function(gam.obj)
### -----------------------------------------------------------------------------------
### gam.exact is a method for the gam class.
###  
### Computes the asymptotically exact variance-covariance matrix for the linear
### terms in the model (except for the intercept).
###
### Note: Use of lo in the model formula is not allowed.  
###
### Author: Aidan McDermott (AMcD)
### Date:   Mar 5, 2003
###
###         Mar 28, 2003
###         Fixed single linear term models -- thanks to Tim Ramsay
###         April 17, 2006
###         Modified to work in R by Trevor Hastie
###
### See:
### 
### [1] Issues in Semiparametric Regression: A Case Study of Time Series
###     Models in Air Pollution and Mortality,
###     Dominici F., McDermott A., Hastie T.J.,
###     Technical Report, Department of Biostatistics, Johns Hopkins University,
###     Baltimore, MD, USA. 
###  
### -----------------------------------------------------------------------------------
  {

    if ( is.na(match("gam",class(gam.obj))) ) {
      stop("not a gam object")
    }
    
    nl.df    <- gam.obj$nl.df
    terms    <- terms(gam.obj)
    at.terms <- attributes(terms)

    coef <- coef(gam.obj)
    
    w   <- gam.obj$weights
    mu  <- gam.obj$fitted.values
    eta <- gam.obj$additive.predictors
    y   <- as.matrix(gam.obj$y)

    family   <- family(gam.obj)
    mu.eta.val <- family$mu.eta(eta)
    z <- eta + (y - mu)/mu.eta.val

    
### Don't want lo in gam formula.
###    if ( length((at.terms$specials)$lo) > 0 ) {
###      stop("lo found in gam formula.")
###    }

    X   <- model.matrix(gam.obj)
    Y   <- as.matrix(gam.obj$y)

### only take terms that survived the original gam call
    names.coef <- names(coef)
    has.intercept <- match("(Intercept)",names.coef)
    if ( !is.na(has.intercept) ) names.coef <- names.coef[-has.intercept]
    X   <- X[,names.coef]
    tnames <- dimnames(X)[[2]]
    form   <- "y~"
    special.list <- c()
### Replace the df with the actual df returned by gam.
### Rewrite fromula to match names in X
    for ( k in 1:length(tnames) ) {
      if ( substring(tnames[k],1,2) == "s(" ) {
        s.call     <- match.call(s,parse(text=tnames[k]))
        this.name  <- as.name(paste("x",k,sep=""))

        which      <- match(tnames[k],names(nl.df))
        if ( is.na(which) ) stop(paste("can't find df for term",tnames[k]))
        this.df    <- nl.df[which]+1

        form <- paste(form,
                      "+s(",this.name,",df =",this.df,")")
        special.list <- c(special.list,k)
      }
      else if ( substring(tnames[k],1,3) == "lo(" ) {
        lname <- length(tnames[k])
        if ( substring(tnames[k],lname,lname) == "1" ) tnames[k] <- substring(tnames[k],1,(lname-1))
        if ( substring(tnames[k],lname,lname) == ")" ) {
        lo.call    <- match.call(lo,parse(text=tnames[k]))
        this.name  <- as.name(paste("x",k,sep=""))

        lo.call[[2]] <- this.name
        lo.call <- deparse(lo.call)
        form <- paste(form,"+",lo.call)
      }
        special.list <- c(special.list,k)
      }
      else form <- paste(form,"+x",k,sep="")
    }
    mydat <- data.frame(cbind(Y,X))
    names(mydat) <- c("y",paste("x",1:ncol(X),sep=""))

    XX <- X
    mydat[,"w"] <- w

    Control <- gam.obj$call$control
    if ( is.null(Control) ) {
      call      <- gam.obj$call
      call[[1]] <- as.name("gam.control")
      Control   <- eval(call,sys.parent())
    }

    for ( k in 1:length(tnames) ) {
      if ( substring(tnames[k],1,2) != "s("  & substring(tnames[k],1,3) != "lo(" ) {
        this.var <- paste("x",k,sep="")
        upd.form <- update(as.formula(form),paste(this.var,"~. -",this.var))

        XX[,k] <- gam(formula=upd.form,data=mydat,family=gaussian,weights=w,
                      control=eval(Control))$fitted
      }
    }

### Need to test we get some data
    if ( length(X) == 0 ) stop("nothing to do")
    
    X   <- X[,-special.list,drop=FALSE]
    sx  <- XX[,-special.list,drop=FALSE]
    swx <- w*sx
    
    if ( length(X) == 0 ) stop("no linear terms in the model -- nothing to do")
    
    A <- t(X) %*% ( w * X ) - t(X) %*% ( w * sx )
    B <- t(X*w) - t(swx)
    H <- solve(A) %*% B

    beta    <- H %*% z
    varbeta <- (H * (1/w)) %*% t(H) * as.vector(summary(gam.obj)$dispersion)
    se      <- sqrt(diag(varbeta))

    coef <- cbind(summary.glm(gam.obj)$coef,NA,NA,NA)
    tab <- cbind(beta,se,beta/se,2*(1-pnorm(beta/se)))
    coef[dimnames(tab)[[1]],c(5,6,7)] <- tab[,c(2,3,4)]

    dimnames(coef) <- list(dimnames(coef)[[1]],
                           c(dimnames(coef)[[2]][1:4],
                             "A-exact SE","A-exact Z","A-exact P")) 

    out.object <- list(coefficients=coef,covariance=varbeta)
    class(out.object) <- c("gamex")
    
    return(out.object)
  }

