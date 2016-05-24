#-------------------------------------------------------------------------------
# This is an the interface to use the lme() function of 
# Pinheiro and Bates (2000)  within gamlss()
# fit random effect  terms using the lme() function of nlme 
# which is used in the backfitting 
#-------------------------------------------------------------------------------
# Authors Mikis Stasinopoulos, Marco Enea
# created  3-6-14 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# TO DO 
# i) the function gives the same results as the re()
# ii) we need to check
#      a) the correlation componet OK
#      b) whether the control elements are functioning OK
#      c) can we use all the martices at table 4.3 page 158?
#      d) multilevel models p 167?
#      e) I am assuming that we can not model the variance function for modelling
#        Hetetoscedasticity page 205 since we have those facilities within gamlss
#      f) corStruct Classes page 233
#          look at table 5.3 page 234 can we fit any of those??
# formula, random = NULL, correlation = NULL,
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
re <-function(fixed=~1, random = NULL, correlation = NULL, method = "ML", ...) 
{ 
#------------------------------------------
# function starts here
#------------------------------------------
    scall <- deparse(sys.call(), width.cutoff = 500L) # 
if (!is(fixed, "formula")) stop("fixed argument in lme() needs a formula starting with ~")
#if (!is(random, "formula")) stop("formula argument in lme() needs a formula starting with ~")
# we have to do somehing with corelation
# if (!is.null(correlation)) {
#   cor.for <- attr(correlation, "formula")
#   if (!is.null(cor.for)) 
#     cor.vars <- all.vars(cor.for)
# }
# else cor.vars <- NULL
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
##--- get the lme control values
control  <- nlme::lmeControl(...)
## get the data
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
   Data <- if (any(attributes(eval(substitute(Data)))$class=="groupedData")) eval(substitute(Data))
else data.frame(eval(substitute(Data))) 
     #===== 
      len <- dim(Data)[1] # get the lenth of the data
## out
     xvar <- rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
   attr(xvar,"fixed")       <- fixed
   attr(xvar,"random")      <- random
   attr(xvar,"method")      <- method  
   attr(xvar,"correlation") <- correlation
   attr(xvar,"control")     <- control
   attr(xvar, "gamlss.env") <- gamlss.env
   if (any(attributes(Data)$class=="groupedData")) {
         attr(xvar, "data") <- Data } else {
         attr(xvar, "data") <- as.data.frame(Data)
   }
   attr(xvar, "call")       <- substitute(gamlss.re(data[[scall]], z, w, ...)) 
   attr(xvar, "class")      <- "smooth"
   xvar
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# re.control <-  function (offset=NULL, method="GCV.Cp",
#     optimizer=c("outer","newton"), control=list(),
#     select=FALSE, knots=NULL, sp=NULL, min.sp=NULL, H=NULL, gamma=1,
#     fit=TRUE, paraPen=NULL, G=NULL, in.out=NULL,... ) 
# {
# 	    control <- gam.control(...)
#        list(offset=offset, method=method, optimizer=optimizer, control=control, select=select, knots=knots, sp=sp, min.sp=min.sp, H=H, gamma=gamma, fit=fit, paraPen=paraPen, G=G, in.out=in.out,...)
# }

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
gamlss.re <-function(x, y, w, xeval = NULL, ...)
{    
      fixed <- attr(x,"fixed")
     random <- attr(x,"random")
correlation <- attr(x,"correlation")
     method <- attr(x,"method")
fix.formula <- as.formula(paste("Y.var", deparse(fixed, width.cutoff = 500L), sep=""))
    control <- as.list(attr(x, "control"))  
#gamlss.env <- as.environment(attr(x, "gamlss.env"))
      OData <- attr(x,"data") 
       Data <-  if (is.null(xeval)) OData #the trick is for prediction
                else  OData[seq(1,length(y)),]
  if (any(attributes(Data)$class=="groupedData")) {
    Data$W.var <- 1/w 
    Data$Y.var <- y 
  } else {
    Y.var <- y
    W.var <- 1/w   
    Data <- data.frame(eval(substitute(Data)),Y.var,W.var) 
  }     
#       Data <- data.frame(eval(substitute(Data)),y,wei=1/w)
# fit  <-  lme(All$fixed, data = Data, random=All$random, weights=varFixed(~wei),  method="ML")
#       
#              (fixed, data = sys.frame(sys.parent()), random, correlation = NULL, 
#           weights = NULL, subset, method = c("REML", "ML"), na.action = na.fail, 
#           control = list(), contrasts = NULL, keep.data = TRUE) 
# lme(fixed = fixed, random = random, data = data, 
#    correlation = correlation, control = control, weights = varFixed(w.formula), 
#    method = "ML", ...)        
       fit <- lme(fixed=fix.formula, data=Data, random=random,  weights=varFixed(~W.var), 
                  correlation = correlation, control=control,  method = method) 
        fv <- fitted(fit)  
 residuals <- y-fv
         N <- sum(w!=0)
        df <-  N-(sum(w*(y-fv)^2))/(fit$sigma^2)                 
  if (is.null(xeval))
    {
   list(fitted.values=fv, residuals=residuals,
     nl.df = df, lambda=fit$sigma, #
     coefSmo = fit, var=NA)    # var=fv has to fixed
    }
else 
    {
   ll<-dim(OData)[1]
   pred <- predict(fit,newdata = OData[seq(length(y)+1,ll),])
    }         

}
#-------------------------------------------------------------------------------

