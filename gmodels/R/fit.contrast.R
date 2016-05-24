# $Id: fit.contrast.R 2060 2015-07-19 03:22:30Z warnes $

fit.contrast.lm <- function(model, varname, coeff, showall=FALSE,
                            conf.int=NULL, df=FALSE, ...)
{
  # check class of model
  if( !(any(class(model) %in% c("lm", "aov", "lme") ) ))
    stop("contrast.lm can only be applied to objects inheriting from 'lm'",
         "and 'lme' (eg: lm,glm,aov,lme).")

  # make sure we have the NAME of the variable
  if(!is.character(varname))
     varname <- deparse(substitute(varname))

  # make coeff into a matrix
  if(!is.matrix(coeff))
    {
       coeff <- matrix(coeff, nrow=1)
     }

  # make sure columns are labeled
  if (is.null(rownames(coeff)))
     {
       rn <- vector(length=nrow(coeff))
       for(i in 1:nrow(coeff))
          rn[i] <- paste(" c=(",paste(coeff[i,],collapse=" "), ")")
       rownames(coeff) <- rn
     }

  # now convert into the proper form for the contrast matrix
  cmat <- make.contrasts(coeff, ncol(coeff) )
  cn <- paste(" C",1:ncol(cmat),sep="")
  cn[1:nrow(coeff)] <- rownames(coeff)
  colnames(cmat) <- cn

  # recall fitting method with the specified contrast
  m <- model$call

  if(is.null(m$contrasts))
    m$contrasts <- list()
  m$contrasts[[varname]] <- cmat

  if(is.R())
    r <- eval(m, parent.frame())
  else
    r <- eval(m)

  # now return the correct elements ....
  if( 'lme' %in% class(model) )
    {
      est <- r$coefficients$fixed
      se  <- sqrt(diag(r$varFix))
      tval <- est/se
      df.lme   <- r$fixDF$X
      retval <- cbind(
                      "Estimate"= est,
                      "Std. Error"= se,
                      "t-value"= tval,
                      "Pr(>|t|)"=  2 * (1 - pt(abs(tval), df.lme)),
                      "DF"=df.lme
                      )

    }
  else if('glm' %in% class(model))
    {
      smodel <- summary.glm(r)
      retval <- cbind(coef(smodel), "DF"=smodel$df[2])
    }
  else # lm, aov
    {
      smodel <- summary.lm(r)
      retval <- cbind(coef(smodel), "DF"=smodel$df[2])
    }

  if( !showall )
    {
      if( !is.R() && ncol(cmat)==1 )
        {
          retval <- retval[varname,,drop=FALSE]
          rownames(retval) <- rn
        }
      else
        {
          rn <- paste(varname,rownames(coeff),sep="")
          ind <- match(rn,rownames(retval))
          retval <- retval[ind,,drop=FALSE]
        }

    }

  if(!is.null(conf.int)) # add confidence intervals
    {
      alpha <- 1-conf.int
      retval <- cbind( retval,
                      "lower CI"=retval[,1] -
                      qt(1-alpha/2,retval[,5])*retval[,2],
                      "upper CI"=retval[,1] +
                      qt(1-alpha/2,retval[,5])*retval[,2] )
    }

  if(!df)
    return(retval[,-5,drop=FALSE])
  else
    return(retval)
}

# fit.contrast.lme and fit.contrast.mer are necessary because
# 'lme' and 'mer' objects do not inherit from 'lm'.
#
# **Make sure that the argument list *exactly* matches the one
# for fit.contrast.lm() above.**
#
fit.contrast.lme <- function(model, varname, coeff, showall=FALSE,
                            conf.int=NULL, df=FALSE, ...)
  {
    fit.contrast.lm(model, varname, coeff, showall, conf.int, df)
  }

## # I made rather dramatic changes here and do all calculations in fit.contrast.mer rather than
## # fit.contrast.lm because of the simulation extras ... added sim.mer and n.sim to the parameter list
## fit.contrast.mer <- function(model, varname, coeff, showall=FALSE,
##                             conf.int=NULL, sim.mer=TRUE, n.sim=1000, ...)
## {
##   require(lme4)

##   # make sure we have the NAME of the variable
##   if(!is.character(varname))
##      varname <- deparse(substitute(varname))

##   # make coeff into a matrix
##   if(!is.matrix(coeff))
##     {
##        coeff <- matrix(coeff, nrow=1)
##      }

##   # make sure columns are labeled
##   if (is.null(rownames(coeff)))
##      {
##        rn <- vector(length=nrow(coeff))
##        for(i in 1:nrow(coeff))
##           rn[i] <- paste(" c=(",paste(coeff[i,],collapse=" "), ")")
##        rownames(coeff) <- rn
##      }

##   # now convert into the proper form for the contrast matrix
##   cmat <- make.contrasts(coeff, ncol(coeff) )
##   cn <- paste(" C",1:ncol(cmat),sep="")
##   cn[1:nrow(coeff)] <- rownames(coeff)
##   colnames(cmat) <- cn

##   m <- model@call

##   if(is.null(m$contrasts))
##     m$contrasts <- list()
##   m$contrasts[[varname]] <- cmat

##   if(is.R())
##     r <- eval(m, parent.frame())
##   else
##     r <- eval(m)
##   # now return the correct elements ....
##   r.effects <- fixef(r)
##   n <- length(r.effects)

##   if(sim.mer)
##   {
##     retval <- est.mer(obj = r, cm = diag(n), beta0 = rep(0, n),
##                        conf.int = conf.int, show.beta0 = FALSE,
##                        n.sim=n.sim)
##     rownames(retval) <- names(r.effects)
##   }else{
##     if(!is.null(conf.int))
##       warning("Confidence interval calculation for mer objects requires simulation -- use sim.mer = TRUE")

##     est <- fixef(r)
##     se  <- sqrt(diag(as.matrix(vcov(r))))
##     tval <- est/se
##     retval <- cbind(
##                     "Estimate"= est,
##                     "Std. Error"= se,
##                     "t-value"= tval
##                     )
##   }

##   if( !showall )
##   {
##     if( !is.R() && ncol(cmat)==1 )
##     {
##       retval <- retval[varname,,drop=FALSE]
##       rownames(retval) <- rn
##     }else{
##       rn <- paste(varname,rownames(coeff),sep="")
##       ind <- match(rn,rownames(retval))
##       retval <- retval[ind,,drop=FALSE]
##     }
##   }

##   return(retval)
## }


fit.contrast <- function(model, varname, coeff, ...)
  UseMethod("fit.contrast")

