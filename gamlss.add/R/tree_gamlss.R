#require(tree)
#require(rpart)
# this functions fits regression trees  
# using backfitting 
# TO DO
# i)  check why for not identical link it is produce more fitted values than exist
# ii) needs prediction
#----------------------------------------------------------------------------------------
tr <-function(formula, method=c("rpart"), control=rpart.control(...), ...) #  ,"tree" 
{
#---------------------------------
    method <- match.arg(method)
     scall <-  deparse(sys.call(), width.cutoff = 200L)
# check the formula
if (!is(formula, "formula")) stop("formula argument in tr() needs a formula starting with ~")
# get where "gamlss" is in system call
# it can be in gamlss() or predict.gamlss()       
rexpr <- grepl("gamlss",sys.calls()) ## 
for (i in length(rexpr):1)
{ 
  position <- i # get the position, we are geting the fist from the last
  if (rexpr[i]==TRUE) break
}
gamlss.env <- sys.frame(position) #gamlss or predict.gamlss# 
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
      xvar <- rep(0, len)
   attr(xvar, "formula") <- formula
   attr(xvar, "method")  <- method
   attr(xvar,"control")  <- control

attr(xvar, "gamlss.env") <- gamlss.env
    attr(xvar, "data")   <- as.data.frame(Data)
   attr(xvar, "call")    <- substitute(gamlss.tr(data[[scall]], z, w, ...)) 
   attr(xvar, "class")   <- "smooth"
   xvar
}
#----------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.tr <-function(x, y, w, xeval = NULL, ...)
{
     formula <- attr(x,"formula")
     formula <- as.formula(paste("y",deparse(formula, width.cutoff = 500L), sep=""))
      method <- attr(x, "method")
     control <- as.list(attr(x, "control"))
  gamlss.env <- as.environment(attr(x, "gamlss.env"))
       OData <- attr(x,"data") 
        Data <-  if (is.null(xeval)) OData #this is for prediction
                 else  OData[seq(1,length(y)),]
        Data <- data.frame(eval(substitute(Data)),y,w) 
          rexpr<-regexpr("gamlss",sys.calls())
         fit <- if (method=="rpart")  
           {
           rpart(formula, data=Data, weights=w, cp = control$cp, minsplit = control$minsplit, minbucket = control$minbucket,  maxcompete = control$maxcompete,  maxsurrogate =  control$maxsurrogate, usesurrogate=control$usesurrogate, xval=control$xval, surrogatestyle=control$surrogatestyle,maxdepth = control$maxdepth)
         }
       #         else  tree(formula, data=Data, weights=w)   
        frame <- fit$frame
       leaves <- frame$var == "<leaf>"
         size <- sum(leaves)
           df <-  size #frame$n[1] - size
    residuals <- resid(fit)
           fv <- predict(fit)  # this can be a matrix ??
  if (is.null(xeval))
    {
   list(fitted.values=fv, residuals=residuals,
     nl.df = df, lambda=NA, ## we nead df's here 
     coefSmo = fit, var=NA)    # var=fv has to fixed
    }
else 
    {
    #stop("not prediction for the tree() function exist yet")
    ndata <-subset(OData, source=="newdata")
    pred <- predict(fit, newdata=ndata)
    }         

}
      
