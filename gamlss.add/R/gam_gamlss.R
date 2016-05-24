# This is the interface to use the gam() function of Simon Wood
# within gamlss()
# fit smoothing terms using the gam function of mgcv 
# which is used in the backfitting 
# Author Mikis Stasinopoulos
# latest change 21--8-12 MS
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
ga <-function(formula, control=ga.control(...),...) 
{ 
#------------------------------------------
# function starts here
#------------------------------------------
    scall <- deparse(sys.call(), width.cutoff = 500L) # 14-oct-2013 DS 
if (!is(formula, "formula")) 
  stop("formula argument in ga() needs a formula starting with ~")
# get where "gamlss" is in system call
# it can be in gamlss() or predict.gamlss()
    rexpr <- grepl("gamlss",sys.calls()) ## 
for (i in length(rexpr):1)
   { 
 position <- i # get the position
 if (rexpr[i]==TRUE) break
   }
  #
gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
##---
## get the data
## this has been modified on the 12-12-14 to make sure that 
##  if model.frame.gamlss() is used as for example in term.plot() the
## function does not fails (It need to implemented to all smoother using formulea?)
if (sys.call(position)[1]=="predict.gamlss()")
     { # if predict is used 
      Data <- get("data", envir=gamlss.env)
     }
else if (sys.call(position)[1]=="gamlss()") 
     { # if gamlss() is used
     if (is.null(get("gamlsscall", envir=gamlss.env)$data)) 
         { # if no data argument but the formula can be interpreted
     Data <- model.frame(formula)	
         }
     else
         {# data argument in gamlss 
     Data <- get("gamlsscall", envir=gamlss.env)$data
         }
     }
else  {Data <- get("data", envir=gamlss.env)}
     Data <- data.frame(eval(substitute(Data)))
     #===== 
      len <- dim(Data)[1] # get the lenth of the data
## out
     xvar <- rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
   attr(xvar,"formula")     <- formula
   attr(xvar,"control")    <- control
   attr(xvar, "gamlss.env") <- gamlss.env
   attr(xvar, "data")       <- as.data.frame(Data)
   attr(xvar, "call")       <- substitute(gamlss.ga(data[[scall]], z, w, ...)) 
   attr(xvar, "class")      <- "smooth"
   xvar
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
ga.control <-  function (offset=NULL, method="GCV.Cp",
    optimizer=c("outer","newton"), control=list(),
    select=FALSE, knots=NULL, sp=NULL, min.sp=NULL, H=NULL, gamma=1,
    fit=TRUE, paraPen=NULL, G=NULL, in.out=NULL,... ) 
{
	    control <- gam.control(...)
       list(offset=offset, method=method, optimizer=optimizer, control=control, select=select, knots=knots, sp=sp, min.sp=min.sp, H=H, gamma=gamma, fit=fit, paraPen=paraPen, G=G, in.out=in.out,...)
}

#method :"GCV.Cp"  "GACV.Cp"  "REML" "P-REML" for REML estimation, but using a Pearson estimate of the scale. "ML" and "P-ML" 
# optimizer: "perf" "outer" f "outer" can use several optimizer: "newton"  "bfgs", "optim", "nlm" and "nlm.fd"
# select: if this is TRUE then gam can add an extra penalty to each term so that it can be penalized to zero. This means that the smoothing parameter estimation that is part of fitting can completely remove terms from the model. 
# knots :  user specified knot values to be used (note that the number of knots is not always just k). See tprs for what happens in the "tp"/"ts" case. Different terms can use different numbers of knots, unless they share a covariate.
# sp : A vector of smoothing parameters can be provided here
# min.sp : lower band for lambdas
# H: A user supplied fixed quadratic penalty on the parameters of the GAM can be supplied, with this as its coefficient matrix.
# gamma : inflate the model degrees of freedom in the GCV or UBRE/AIC
# fit : TRUE
#paraPen Penanlits in the parameters
#G usally NULL
# in.out list(sp, scale)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.ga <-function(x, y, w, xeval = NULL, ...)
{     
   formula <- attr(x,"formula")
   formula <- as.formula(paste("Y.var", deparse(formula, width.cutoff = 500L), sep=""))
  control <- as.list(attr(x, "control"))  
#gamlss.env <- as.environment(attr(x, "gamlss.env"))
      OData <- attr(x,"data") 
       Data <-  if (is.null(xeval)) OData #the trick is for prediction
                else  OData[seq(1,length(y)),]
      Y.var <- y
      W.var <- w   
       Data <- data.frame(eval(substitute(Data)),Y.var,W.var)  
       fit <- gam(formula, data=Data, weights=W.var, 
                  offset=control$offset, method=control$method,
                  optimizer=control$optimizer, control=control$control, 
                  select=control$select, knots=control$knots,sp=control$sp,
                   min.sp=control$min.sp, H=control$H, gamma=control$gamma,
                   fit=TRUE, paraPen=control$paraPen, G=control$G, 
                   in.out=control$in.out)
        df <- sum(fit$edf)-1 
        fv <- fitted(fit) 
 residuals <- y-fv
  #     var <- predict(fit, se=T)$se 
   # the se of the fit can be obtained  using predict(getSmo(m1))$se
  if (is.null(xeval))
    {
   list(fitted.values=fv, residuals=residuals,
     nl.df = df, lambda=fit$sp[1], #
     coefSmo = fit, var=NA)    # var=fv has to fixed
    }
else 
    {
   ll<-dim(OData)[1]
   pred <- predict(fit,newdata = OData[seq(length(y)+1,ll),])
    }         

}
      
