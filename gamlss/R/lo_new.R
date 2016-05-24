##---------------------------------------------------------------------------------------
## the lo(), lo.control() and gamlss.lo() functions
## are based on R loess() and S-plus lo() 
## author Mikis Stasinopoulos
##------------------------------------------------------------------------------
## created by MS : 13-08-12
## based on Brian Ripley impementation of loess in R function
## This replace the previous version 2002 version
#------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
lo <-function(formula, control=lo.control(...), ...) 
{ 
#------------------------------------------
# function starts here
#------------------------------------------
    scall <- deparse(sys.call(), width.cutoff = 500L)
if (!is(formula, "formula")) stop("lo() needs a formula starting with ~")
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

      len <- dim(Data)[1] # get the lenth of the data
    #  get the free arguments 
    alist <- list(...)
    # check if df are defined
    if (!is.null(alist$df))  control$enp.target <- alist$df 
## out
     xvar <- rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
   attr(xvar,"formula")     <- formula
    attr(xvar,"control")    <- control 
   attr(xvar, "gamlss.env") <- gamlss.env
   attr(xvar, "data")       <- as.data.frame(Data)
   attr(xvar, "call")       <- substitute(gamlss.lo(data[[scall]], z, w, ...)) 
   attr(xvar, "class")      <- "smooth"
   xvar
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
lo.control <-  function (span = 0.75, enp.target=NULL, degree = 2,
      parametric = FALSE, drop.square = FALSE, normalize = TRUE,
      family = c("gaussian", "symmetric"),
      method = c("loess", "model.frame"),
      surface = c("interpolate", "direct"), 
      statistics = c("approximate", "exact", "none"), 
      trace.hat = c("exact", "approximate"), 
      cell = 0.2, 
      iterations = 4,  iterTrace = FALSE, ...) 
{
       list(span = span, enp.target=enp.target, degree=degree,
         parametric = parametric, drop.square = drop.square, normalize = normalize, family=family,
         method=method, surface = surface, statistics = statistics, trace.hat = trace.hat, 
         cell = cell, iterations = iterations, iterTrace = FALSE)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.lo <-function(x, y, w, xeval = NULL, ...)
{           
    formula <- attr(x,"formula")
    formula <- as.formula(paste("Y.var",deparse(formula), sep=""))
    control <- as.list(attr(x, "control"))
#gamlss.env <- as.environment(attr(x, "gamlss.env"))
      OData <- attr(x,"data") 
       Data <-  if (is.null(xeval)) OData #the trick is for prediction
               else  OData[seq(1,length(y)),]
       Y.var <- y
       W.var <- w     
      Data <- data.frame(eval(substitute(Data)),Y.var,W.var)   
       if (is.null(control$enp.target))
        { 
       fit <- loess(formula, data=Data, weights=W.var, 
                     span=control$span, 
                     degree=control$degree, normalize=control$normalize,
                     family=control$family, 
                     control=loess.control(surface=control$surface,
                     statistics=control$statistics,
                     trace.hat=control$trace.hat,
                     iterations= control$iterations, 
                     iterTrace = control$iterTrace)) 
         }
       else 
         { 
       	fit <-loess(formula, data=Data, weights=W.var,
                     enp.target=control$enp.target, 
                     degree=control$degree, normalize=control$normalize,
                     family=control$family, 
                     control=loess.control(surface=control$surface,
                     statistics=control$statistics,
                     trace.hat=control$trace.hat,
                     iterations= control$iterations,
                     iterTrace = control$iterTrace)) 
         } 
        df <- fit$trace.hat-1 
        fv <- fitted(fit) 
 residuals <- Y.var-fv
       var <- NA  #if(control$se) (predict(fit, se=T)$se.fit)^2 else NA
  if (is.null(xeval))
    {
   list(fitted.values=fv, residuals=residuals,
     nl.df = df, lambda=fit$pars[["span"]], ## we nead df's here 
     coefSmo = fit, var=var) # for more than one  
    }
else 
    {
   ll<-dim(OData)[1]
   pred <- predict(fit,newdata = OData[seq(length(y)+1,ll),])
    }         
}
#------------------------------------------------------------------------------
# to do 
# i) can we put points in 3-dim plots?
# ii) maybe we should put SE in 3-d
#--------------------------------------------------------------------------------      
vis.lo<- function(obj,
                  se = -1,
                  rug = FALSE,
                  partial.resid = FALSE,
                  col.term = "darkred",
                  col.shaded = "gray", 
                  col.res = "lightblue", 
                  col.rug = "gray",
                  lwd.term = 1.5,
                  cex.res = 1, 
                  pch.res = par("pch"),
                  type = c("persp", "contour"), 
                  col.surface = "gray", 
                  nlevels = 30, 
                  n.grid = 30, 
                  image = TRUE, 
                  ...)
{
  #-------------------------------------------------------------------------------  
  se.shaded <- function(x, ff = 2) 
  {
    se <- ff * x$se.fit[oo]
    yy <- x$fit[oo]
    xx <- obj$x[oo]
    rx <- c(xx,rev(xx))
    ry <- c(yy - se, rev(yy + se))
    polygon(rx, ry, col = col.shaded, border = col.shaded)
  }
  #-------------------------------------------------------------------------------
  # checking whether loess
  if  (!is(obj,"loess")) stop("The objects should be loess class")
  # getting the number of x-variables
  xnames <- obj$xnames
  lnames <- length(xnames)
  # stop if more than three
  if (lnames>=3) stop("currently only up to two x-variables are plotted")
  # if only one use the usual term.plot
  if (lnames==1)
  {
    #   if (se>=0) # always gives se in this case 
    #   {
    oo <- order(obj$x)
    xx <- obj$x[oo]
    objPred <- predict(obj, term=obj$xnames[[1]], se=TRUE)
    plot(obj$x[oo], obj$y[oo], type = "n", xlab = xnames[[1]], 
         ylab = "fitted values", xlim = range(obj$x), ylim = range(obj$y), 
         main = "loess fit", col = "black", lwd = lwd.term, ...) 
    se.shaded(objPred)  
    lines(obj$x[oo], obj$fitted[oo], col = col.term,lwd = lwd.term,...) 
    if (partial.resid) 
    {
      points(obj$x[oo], obj$y[oo], cex = cex.res, pch = pch.res, 
             col = col.res)
    }
    if (rug) 
    {
      ylims = range(obj$y)
      n <- length(xx)
      lines(rep.int(jitter(xx), rep.int(3, n)), rep.int(ylims[1] + 
                                                          c(0, 0.05, NA) * diff(ylims), n), col=col.rug)
    }  
    #   }
  } # finish the one variable plot
  if (lnames==2) # three dimensional plots here
  {
    # ckeck the type of 3-d plot
    type  <- match.arg(type)
    rx1 <- range(obj$x[,1]) # range of x1
    rx2 <- range(obj$x[,2]) # range of x2
    rz <- range( obj$y)    # range of y
    # get agrid        
    EG <- expand.grid(x1= seq(rx1[1], rx1[2], length=n.grid), 
                      x2= seq(rx2[1], rx2[2], length=n.grid))
    colnames(EG) <- xnames # replace the real names 
    newData <- data.frame(EG) # the data frame
    predValues <- predict(obj, newdata=newData)   # we do not dneed se here yet                
    X1 <- seq(rx1[1], rx1[2], length=n.grid)
    X2 <- seq(rx2[1], rx2[2], length=n.grid)
    z <-  matrix(predValues, nrow=length(X1))
    # rz <- range(z)
    if (type=="contour")  # contour 
    {
      if (image)
      {
        image(X1, X2, z, xlab=xnames[[1]], ylab=xnames[[2]], ...) 
        contour(X1, X2, z, add=TRUE, 
                nlevels=nlevels, xlab=xnames[[1]], ylab=xnames[[2]])
      } else 
        contour(X1, X2, z, nlevels=nlevels, xlab=xnames[[1]], ylab=xnames[[2]])
    } else # if 
    {
      if (se>=0)       # if persp and need se plors
      {
        predValues <- predict(obj, newdata=newData, se=TRUE) #
        # lower se plot
        z1 <-  matrix(predValues$fit - se*predValues$se.fit , nrow=length(X1))
        persp(X1, X2, z1, col=NA,
              zlab="partial fit",   zlim=c(rz[1],rz[2]), 
              xlab=xnames[[1]], ylab=xnames[[2]], border="gray",...)
        # the fitted values 
        par(new = TRUE)
        persp(X1, X2, z, col=NA, border="black",
              zlab="partial fit",   zlim=c(rz[1],rz[2]), 
              xlab=xnames[[1]], ylab=xnames[[2]])
        # the upper se plot
        z2 <-  matrix(predValues$fit + se*predValues$se.fit , nrow=length(X1))
        par(new = TRUE)
        res<-persp(X1, X2, z2, col=NA, border="gray",
                   zlab="partial fit",   zlim=c(rz[1],rz[2]), 
                   xlab=xnames[[1]], ylab=xnames[[2]],...)
        if (partial.resid)  # if partial plot the data
        {
          suppressWarnings(
          points(trans3d(EG[,1], EG[,2], obj$y, pmat=res ), col =col.res, 
                 pch = pch.res, cex = cex.res))
        }
      } else # if persp and no se's plots'
      {
        res<-  persp(X1, X2, matrix(predValues, nrow=length(X1)), 
                     zlab="fitted values",  
                     zlim=c(rz[1],rz[2]), xlab=xnames[[1]], ylab=xnames[[2]],
                     col=col.surface,
                     ...)
        if (partial.resid) # if partial plot the data
        {
          suppressWarnings(
          points(trans3d(EG[,1], EG[,2], obj$y, pmat=res ), col =col.res, 
                 pch = pch.res, cex = cex.res))
        }  
      } # end  if persp and no se's plots'
    }  
    
  }
}
