## TO DO clear and put in gamlss
##   i) check the term option 
##  ii) NA coefficients should be checked
##  iii) what happends if data is not declared   
##-------------------------------------------------------------------------------------------
## last change 12-04-11 DS
## this corrects a bug the reported by Gillian Heller where offset was used to fit the model 
## the following comment is old 
# there is a problem with predict and  lpred
# the type = "response" is not working in
# the link is not the default. The reason is that in 
#  gamlss.family(family(object)[1])[[paste(what,"linkinv",sep=".")]](pred)
# family(object)[1] do not pick up if the link is not the default   
# maybe this will do
#paste(family(m1)[1],"(",what,".link=",eval(parse(text=(paste("family$",what,".link", sep="")))),")", sep="")
#--------------------------------------------------------------------------------------------- 
### The predict.gamlss function 
### Thursday, June 17, 2004 at 10:23 
### it is based on the safe.predict.gam() of Trevor Hastie
### author: Mikis Stasinopoulos
### last change Thursday,  12-04-11 DS
### BUGS
### 1) the type = "terms" is not working : fixed
### 2) offset is not working : this is fixed but needs checking
### 3) se.fit is not supported for new data at all additves terms even for just linear 
### 4) NA coefficients are not working yet : should be OK
predict.gamlss <- function(object, 
                           what = c("mu", "sigma", "nu", "tau"), 
                            parameter = NULL,
                           newdata = NULL, 
                           type = c("link", "response", "terms"), # terms not working 
                           terms = NULL, 
                           se.fit = FALSE, 
                           data = NULL, ...)                                                                  
{
## this little function put data frames together 
##  originated from an the R-help reply by B. Ripley
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##-------- concat starts here
concat <- function(..., names=NULL) 
  { 
    tmp <- list(...) 
    if(is.null(names)) names <- names(tmp) 
    if(is.null(names)) names <- sapply( as.list(match.call()), deparse)[-1] 
    if( any( 
            sapply(tmp, is.matrix) 
            | 
            sapply(tmp, is.data.frame) ) ) 
      { 
        len <- sapply(tmp, function(x) c(dim(x),1)[1] ) 
        len[is.null(len)] <- 1 
        data <- rbind( ... ) 
      } 
    else 
      { 
        len <- sapply(tmp,length) 
        data <- unlist(tmp)     
      } 
    namelist <- factor(rep(names, len), levels=names)          
    return( data.frame( data, source=namelist) ) 
  } 
##----------concat finish here
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##   main function starts here
##----------------------------
## If no new data just use lpred() and finish
if (is.null(newdata))  # 
    {
    predictor<- lpred(object, what = what, type = type, terms = terms, se.fit = se.fit, ... )
    return(predictor)
    }
## at the moment se.fit is not supported for new data
if (se.fit) 
    warning(" se.fit = TRUE is not supported for new data values at the moment \n")
##  stop if newdata is not data frame
## note that atomic is not working here so better to take it out Mikis 23-10-13 
## if (!(is.atomic(newdata) | inherits(newdata, "data.frame")))
  if (!(inherits(newdata, "data.frame")))
    stop("newdata must be a data frame ") # or a frame mumber
## getting which parameter and type   
       what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
       type <- match.arg(type)
## get the original call 
       Call <- object$call
## we need both the old and the new data sets
## the argument data can be provided by predict
data<- data1 <- if (is.null(data))
{        ## if it is not provided then get it from the original call
            if (!is.null(Call$data)) eval(Call$data) 
            else stop("define the original data using the option data") 
           }
        else data # if it provide get it 
## keep only the same variables 
## this assumes that all the relevant variables will be in newdata
## what happens if not?
      data <- data[match(names(newdata),names(data))]  
## merge the two data together
      data <- concat(data,newdata)
## get the formula 
   parform <- formula(object, what)# object[[paste(what, "formula", sep=".")]]
## put response to NULL
if (length(parform)==3)
  parform[2] <- NULL
## define the terms 
     Terms <- terms(parform)
## get the offset
 offsetVar <- if (!is.null(off.num <- attr(Terms, "offset"))) # new 
             eval(attr(Terms, "variables")[[off.num + 1]], data) 
## model frame 
         m <- model.frame(Terms, data, xlev = object[[paste(what,"xlevels",sep=".")]])           
             ## model design matrix y and w 
         X <- model.matrix(Terms, data, contrasts = object$contrasts)
         y <- object[[paste(what,"lp",sep=".")]] 
         w <- object[[paste(what,"wt",sep=".")]] 
## leave for future checks
       #  aN <- dim(newdata)[1]
      #zeros <- rep(0,aN)
       #ones <- rep(1,aN)
       #yaug <- as.vector(c(y,zeros))
       #waug <- as.vector(c(w,zeros))
## for keeping only the original data
  onlydata <- data$source == "data" # TRUE or FALSE
## whether additive terms are involved in the fitting 
   smo.mat <- object[[paste(what,"s",sep=".")]]
## if offset take it out from fitting
if (!is.null(off.num))
         y <- (y - offsetVar[onlydata])
## if smoothing 
if (!is.null(smo.mat))
      {
     n.smooths <- dim(smo.mat)[2]
             y <- (y - smo.mat %*% rep(1, n.smooths))
      }
## refit the model
    refit <- lm.wfit(X[onlydata,  , drop = FALSE], y, w)
## ckeck the residuals if they are zero
##if (any(abs(resid(refit))>1e-005)) 
if (abs(sum(resid(refit)))>1e-001||abs(sum(coef(object, what=what)-coef(refit), na.rm=TRUE))>1e-005)
warning(paste("There is a discrepancy  between the original and the re-fit",
               " \n used to achieve 'safe' predictions \n ", sep = "" ))  
## this is disturbing fit and refit have different coefficients  why?
## fit <- lm.wfit(X, yaug, waug)
## get the coefficients
       coef <- refit$coef         ## save the coefficints
         nX <- dimnames(X)        ## the names of rows and columns 
   rownames <- nX[[1]][!onlydata] ## only the newdata rows
      nrows <- sum(!onlydata)     ## the number of rows in the new data
        nac <- is.na(coef)        ## whether they are NA in coefficients
assign.coef <- attr(X, "assign")  ## X is a matrix
   collapse <- type != "terms"## !collapse is for whether type is not "terms"     
      Xpred <- X[!onlydata,]
      Xpred <- matrix(Xpred, nrow=nrows) # I think this probably is not needed sinse allready a matrix
# I will check this later      
  if (!collapse)       ## whether type=="terms" 
    { 
          aa <- attr(X, "assign")
          ll <- attr(Terms, "term.labels")
     if (attr(Terms, "intercept") > 0)  ll <- c("(Intercept)", ll)
         aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
hasintercept <- attr(Terms, "intercept") > 0
           p <- refit$qr$rank  
          p1 <- seq(len = p)
         piv <- refit$qr$pivot[p1]
     if (hasintercept) 
        {
     asgn$"(Intercept)" <- NULL
                    avx <- colMeans(X[onlydata, ])
             termsconst <- sum(avx[piv] * coef[piv])
        }
      # TT <- sum(onlydata)
      # xbar <- drop(array(1/TT, c(1, TT)) %*% X[onlydata, !nac])
      nterms <- length(asgn)
  #    if (nterms > 0) 
  # define the prediction matrix
        pred <- matrix(ncol = nterms, nrow = nrows)
dimnames(pred) <- list(rownames(newdata), names(asgn))
  #          if (se.fit ) 
  #          {
  #              ip <- matrix(ncol = nterms, nrow = NROW(X))
  #              dimnames(ip) <- list(rownames(X), names(asgn))
  #              Rinv <- qr.solve(qr.R(obj[[paste(what,"qr",sep=".")]])[p1, p1])
  #          }
       if (hasintercept) 
        Xpred <- sweep(Xpred, 2, avx)
        unpiv <- rep.int(0, NCOL(Xpred))
   unpiv[piv] <- p1
            for (i in seq(1, nterms, length = nterms)) 
            {
             iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
    iipiv[ii == 0] <- 0
         pred[, i] <- if (any(iipiv > 0)) # ms Thursday, May 1, 2008 at 10:12
                      Xpred[, iipiv, drop = FALSE] %*% coef[iipiv]
                      else 0
  #              if (se.fit ) 
  #                ip[, i] <- if (any(iipiv > 0)) 
  #                  as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii, 
  #                    , drop = FALSE])^2 %*% rep.int(1,p)
  #                else 0
            }
  attr(pred, "constant") <- if (hasintercept)  termsconst
                            else 0     
     #   Xpred <- Xpred - outer(rep(1, sum(!onlydata)), avx)
     if (!is.null(terms)) 
            {   
                pred <- pred[, terms, drop = FALSE]
             #   if (se.fit) 
             #     ip <- ip[, terms, drop = FALSE]
            }
    } ## end if for terms 
    else  ## if type is not terms but "link" or "response" 
    { 
     pred <- drop(Xpred[, !nac, drop = FALSE] %*% coef[!nac])
    if (!is.null(off.num) && collapse)
     pred <- pred + offsetVar[!onlydata]   
    }
## 
## now the smoothing part
##
if (!is.null(smo.mat))
      {
      cat("new prediction", "\n")
       smooth.labels <- dimnames(smo.mat)[[2]]       ## getting the labels i.e. "pb(Fl)" "pb(A)"
              pred.s <- array(0, c(nrows, n.smooths), list(names(pred), 
                            dimnames(smo.mat)[[2]])) ## creating the prediction matrix 
                                 # smooth.labels[smooth.labels%in%colnames(X)]    
                                 # smooth.wanted <- smooth.labels[match(smooth.labels, colnames(X), 0) > 0] 
                                                               ## getting the smoothing call
        smooth.calls <- lapply(m[smooth.labels], attr, "call") # i.e $`pb(Fl)`
                                                               #     gamlss.pb(data[["pb(Fl)"]], z, w)
                data <- subset(m, onlydata, drop=FALSE)        ## get the  original data
 attr(data, "class") <- NULL                                   ## note that m is the data.frame with all data 
               new.m <- subset(m, !onlydata, drop=FALSE)       ## get the new data
attr(new.m, "class") <- NULL
           residuals <-  if (!is.null(off.num)) object[[paste(what,"wv",sep=".")]] - object[[paste(what,"lp",sep=".")]]+offsetVar[onlydata]
                         else object[[paste(what,"wv",sep=".")]] - object[[paste(what,"lp",sep=".")]]
       for(TT in smooth.labels)
         { 
            if (is.matrix(m[[TT]])) # the problem is that for some smoother the m[[TT]] is a matrix (for example pvc())
             { # MS 27-6-11         # in this case  we have to protect the dim attributes of data[[tt]]
              nm <- names(attributes(m[[TT]])) # first we get the names of all attributes 
              attributes(data[[TT]]) <- attributes(m[[TT]])[nm[-c(1,2)]]# then we pass all but
             }                                 # 1 and 2 i.e. dim and names
            else   attributes(data[[TT]]) <- attributes(m[[TT]])
                  Call <- smooth.calls[[TT]] # 
            Call$xeval <- substitute(new.m[[TT]], list(TT = TT))
                     z <- residuals + smo.mat[, TT]
     # debug(gamlss.pvc)
          pred.s[, TT] <- eval(Call)
         }
        if(type == "terms")
            {
            # pred[, smooth.wanted] <- pred[, smooth.wanted] + pred.s[, smooth.wanted]
            pred[, smooth.labels] <- pred[, smooth.labels] + pred.s[, smooth.labels]
            }  
        else pred <- drop(pred + pred.s %*% rep(1, n.smooths)) 
      }
if(type == "response") 
   {
   ## Saturday, April 14, 2007  change to work with trun and cens
    # pred <- try(gamlss.family(family(object)[1])[[paste(what,"linkinv",sep=".")]](pred) , silent = TRUE) 
    #         if (any(class(pred)%in%"try-error"))
    #             { 
    # this is the latest change DS: Monday, March 10, 2008 at 10:00
     if (is(eval(parse(text=object$family[[1]])),"gamlss.family"))
      {
      pred<- eval(parse(text=object$family[[1]]))[[ paste(what,"linkinv",sep=".")]](pred)
      }
     else
      { 
      pred <- gamlss.family(eval(parse(text=paste(family(object)[1],"(",what,".link=",# ms 
                 eval(parse(text=(paste("object$",what,".link", sep="")))),")", sep=""))
                                           ))[[paste(what,"linkinv",sep=".")]](pred) 
      }            
     #pred <- gamlss.family(family(object)[1])[[paste(what,"linkinv",sep=".")]](pred) 
   }       
pred
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# this is to for later release of GAMLSS
# allows the user to get all the parameters using the  predict.gamlss().
# creates a list containing y if exist in the newdata 
# and predicted mu sigma nu and tau 
#----------------------------------------------------------------------------------------
predictAll <-function(object, 
                    newdata = NULL, 
                       type = c("response", "link", "terms"),# note that default is "response" 
                      terms = NULL,   
                     se.fit = FALSE, 
                             ...)       
 {
  type <- match.arg(type)
## if no new data then give all the fitted from the old
 if (is.null(newdata))  # 
    {
    out <- list(y=object$y)
      if ("mu" %in% object$par)  
         out$mu <- lpred(object, what = "mu", type = type, terms = terms, se.fit = se.fit, ... )
     if ("sigma" %in% object$par)  
      out$sigma <- lpred(object, what = "sigma", type = type, terms = terms, se.fit = se.fit, ... )
     if (  "nu" %in% object$par)     
         out$nu <- lpred(object, what = "nu", type = type, terms = terms, se.fit = se.fit, ... )
     if ( "tau" %in% object$par)    
        out$tau <- lpred(object, what = "tau", type = type, terms = terms, se.fit = se.fit, ... )
        attr(out, "family") <- object$family
    return(out)
    }
  else
    {
      out <- list()
      if ("mu" %in% object$par) #
         out$mu <- predict(object,newdata=newdata, what = "mu", type = type, terms = terms, se.fit = se.fit, ... )
     if ("sigma" %in% object$par)  
      out$sigma <- predict(object, newdata=newdata, what = "sigma", type = type, terms = terms, se.fit = se.fit, ... )
     if ("nu" %in% object$par)  
         out$nu <- predict(object, newdata=newdata, what = "nu", type = type, terms = terms, se.fit = se.fit, ... )
     if ("tau" %in% object$par)  
        out$tau <- predict(object, newdata=newdata, what = "tau", type = type, terms = terms, se.fit = se.fit, ... )
     if (as.character(object$mu.formula[[2]])%in%names(newdata)) 
          out$y <-  newdata[,as.character(object$mu.formula[[2]])]
     attr(out, "family") <- object$family
     #out<- list(out,  family=object$family, parameters=object$parameters,  call=object$call, 
     #         weights=object$weights, G.deviance=object$G.deviance, N=object$N, type=object$type, 
     #         #residuals=object , 
     #         noObs=object$noObs,df.fit=object$df.fit, df.residual=object$df.residuals)
     #      class(out) <- c("gamlssPredict", "gamlss")
       return(out)
     }  
 }
#----------------------------------------------------------------------------------------    
  
