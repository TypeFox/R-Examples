
# fit neural networks terms using the nnet function of nnet 
# which is used in the backfitting 
#----------------------------------------------------------------------------------------
# this function need several argument which can not pass throught ...
#x, y, weights, size, Wts, mask,
#     maxit = 100, Hess = FALSE, trace = TRUE, MaxNWts = 1000,
#     abstol = 1.0e-4, reltol = 1.0e-8, 
nn <-function(formula, control=nn.control(...), ...)           
{ 
#------------------------------------------
# function starts here
#------------------------------------------
    scall <-  deparse(sys.call(), width.cutoff = 500L)
# check the formula
 if (!is(formula, "formula")) stop("formula argument in nn() needs a formula starting with ~")
# get where "gamlss" is in system call
# it can be in gamlss() or predict.gamlss()       
    rexpr <- grepl("gamlss",sys.calls()) ## 
for (i in length(rexpr):1)
   { 
 position <- i # get the position, we are geting the fist from the last
 if (rexpr[i]==TRUE) break
   }
gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
##---
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
Data <- data.frame(eval(substitute(Data)))
#===== 
len <- dim(Data)[1] # get the lenth of the data
## what we are doing with starting values 
## here I save only a logical variable in the envir of gamlss
## in gamlss.nn the Wts (coefficients for nnet) are saved 
               sl <-sample(letters, 4)
      fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
    ifStartedName <- paste("if.started",fourLetters, sep=".")
   assign(ifStartedName, FALSE, envir=gamlss.env)  
## out
                          xvar <- rep(0,  len) 
   attr(xvar,"formula")        <- formula
   attr(xvar,"control")        <- control
   attr(xvar, "gamlss.env")    <- gamlss.env
   attr(xvar, "ifStartedName") <- ifStartedName
   attr(xvar, "data")          <- as.data.frame(Data)
   attr(xvar, "call")          <- substitute(gamlss.nn(data[[scall]], z, w, ...)) 
   attr(xvar, "class")         <- "smooth"
   xvar
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# mask logical vector indicating which parameters should be optimized (default all).  
# size number of units in the hidden layer. Can be zero if there are skip-layer units.  
# linout switch for linear output units. Default logistic output units.  
# entropy switch for entropy (= maximum conditional likelihood) fitting. Default by least-squares.  
# softmax switch for softmax (log-linear model) and maximum conditional likelihood fitting. 
#          linout, entropy, softmax and censored are mutually exclusive.  
# censored A variant on softmax, in which non-zero targets mean possible classes. 
#          Thus for softmax a row of (0, 1, 1) means one example each of classes 2 and 3, 
#          but for censored it means one example whose class is only known to be 2 or 3.  
# skip switch to add skip-layer connections from input to output.  
# rang  Initial random weights on [-rang, rang]. Value about 0.5 unless the inputs are large, 
#       in which case it should be chosen so that rang * max(|x|) is about 1.  
nn.control <- function(   size = 3,  linout = TRUE, entropy = FALSE, softmax = FALSE, #Wts=NULL,
                      censored = FALSE, skip = FALSE, rang = 0.7, decay = 0,
                         maxit = 100, Hess = FALSE, trace = FALSE, MaxNWts = 1000,
                        abstol = 1.0e-4, reltol = 1.0e-8)                        
{ 
        if(size < 0) {
warning("the value of size supplied is negative the default value of 3 was used instead")
                order <- 3}
        if(rang <= 0) {
warning("the value of rang supplied is less than zero the default value of 0.7 was used instead")
                rang <- 0.7}
        if(decay < 0) {
warning("the value of decay supplied is negative the default value of 0 was used instead")
                decay <- 0}
        if(maxit < 0) {
warning("the value of maxit supplied is negative the default value of 100 was used instead")
                maxit <- 100}
        if(MaxNWts < 0) {
warning("the value of MaxNWts supplied is negative the default value of 1000 was used instead")
                MaxNWts <- 1000} 
        if(abstol < 0) {
warning("the value of abstol supplied is negative the default value of  1.0e-4 was used instead")
                abstol <- 1.0e-4}
        if(reltol < 0) {
warning("the value of reltol supplied is negative the default value of  1.0e-8 was used instead")
                reltol <- 1.0e-8}                                                          
        list(   size = size, linout =  as.logical(linout)[1], entropy= as.logical(entropy)[1], 
             softmax = as.logical(softmax)[1], #mask=as.logical(mask)[1],
            censored = as.logical(censored)[1], skip = as.logical(skip)[1], 
                rang = rang, decay = decay, maxit = maxit, Hess = as.logical(Hess)[1],
                trace = as.logical(trace)[1], MaxNWts= MaxNWts,  abstol = abstol, reltol = reltol
               #, Wts=Wts
                )
}
#----------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.nn <-function(x, y, w, xeval = NULL, ...)
{       
        formula <- attr(x,"formula")
        formula <- as.formula(paste("y",deparse(formula, width.cutoff = 500L), sep=""))
        control <- as.list(attr(x, "control"))
  ifStartedName <- as.character(attr(x, "ifStartedName")) # get the of the logical 
  theOldWtsName <- paste("theOldWts",substr(ifStartedName, 11, 15),sep="") # create a name for the Wts
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
      ifStarted <-  if (is.null(xeval))  get(ifStartedName, envir=gamlss.env) # 
                    else  TRUE
          OData <- attr(x,"data") 
           Data <-  if (is.null(xeval)) OData #this is for prediction
                   else  OData[seq(1,length(y)),]
           Data <- data.frame(eval(substitute(Data)),y,w)
## whether a proper fit or prediction
if (is.null(xeval)) # proper fit
  {
   if (!ifStarted) # depends on whether we have saved or not Wts
    { 
       fit <- nnet(formula, data=Data, weights=w, 
                   size=control$size,  linout=control$linout, 
                   entropy = control$entropy, softmax = control$softmax,
                   censored = control$censored, skip = control$skip, rang = control$rang, decay = control$decay,
                   maxit = control$maxit, Hess = control$Hess, trace = control$trace, MaxNWts = control$MaxNWts,
                  abstol = control$abstol, reltol = control$reltol )
      assign(ifStartedName, TRUE, envir=gamlss.env)
      assign(theOldWtsName, fit$wts, envir=gamlss.env)
     }
     else # here we have save the Wts in the first fit
     {
     theWts <-  get(theOldWtsName, envir=gamlss.env)
     fit <- nnet(formula, data=Data, weights=w, Wts=theWts, 
                   size=control$size,  linout=control$linout, 
                   entropy = control$entropy, softmax = control$softmax,
                   censored = control$censored, skip = control$skip, rang = control$rang, decay = control$decay,
                   maxit = control$maxit, Hess = control$Hess, trace = control$trace, MaxNWts = control$MaxNWts,
                  abstol = control$abstol, reltol = control$reltol )
      assign(theOldWtsName, fit$wts, envir=gamlss.env)
     }     
            df <- length(coef(fit))# ??? I am not sure about that what happents if decay is not zer0?? 
            fv <- fitted(fit) 
     residuals <- resid(fit)
fit$call$formula <- formula
   list(fitted.values=fv, residuals=residuals,
     nl.df = df, lambda=NA,  coefSmo = fit, var=NA)    #
  }
else # predict 
  { 
       obj <- get("object", envir=gamlss.env ) # get the object from predict
      # now the problem is to pick the right parameter
      # and also the right smoother within the parameter, that is,
      # obj$??.coefSmo[[??]] 
      # paste("obj$", get("what", envir=gamlss.env), ".coefSmoo[[1]]", sep="")
      # match(TT,SL) will finf the number in [[??]]
       #paste("obj$", get("what", envir=gamlss.env), ".coefSmoo[[",as.character(match(TT,SL)), "]]", sep="")
        TT <- get("TT", envir=gamlss.env ) # get wich position is now
        SL <- get("smooth.labels", envir=gamlss.env) # all the labels of the smoother
        NN <- eval(parse(text=paste("obj$", get("what", envir=gamlss.env), ".coefSmo[[",as.character(match(TT,SL)), "]]", sep="")))
    theWts <- coef(NN) # get the Wts
       fit <- nnet(formula, data=Data, weights=w, Wts=theWts, # fit using only one iteration  
                   size=control$size,  linout=control$linout, 
                   entropy = control$entropy, softmax = control$softmax,
                   censored = control$censored, skip = control$skip, rang = control$rang, decay = control$decay,
                   maxit = 1, Hess = control$Hess, trace = control$trace, MaxNWts = control$MaxNWts,
                  abstol = control$abstol, reltol = control$reltol )
   ll<-dim(OData)[1]
   newdata <- OData[seq(length(y)+1,ll),]
   pred <- predict(fit,newdata=newdata) # now predict
  }         

}
#----------------------------------------------------------------------------------------
