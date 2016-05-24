#' @title Linear Regression

#' @description Regression analysis (one numerical predictor variable) with simplified output.
#'   Wrapper function for \code{lm} in package \code{stats}.
#' 
#' @rdname lmGC
#' @usage lmGC(form,data=parent.frame(),graph=FALSE,check=FALSE)
#' @param form formula of form y~x, both variables numeric
#' @param data dataframe supplying y and x above.  If one or more of the variables is not in data, then
#' they will be searched for in the parent environment.
#' @param graph Produce scatterplot with fitted ploynomial, together with prediction standard error bands
#' @param check Asks to produce a lowess or gam curve with approximate 95%-confidence band.  If the
#' fitted line wanders outside the band, then perhaps a linear fit is not appropriate.
#' @return A list of class "GClm".  Elements that may be queried include "slope", "intercept",
#' "s" (residual standard error), "R^2" (unadjusted).
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #To study the relationship between two numerical variables:
#' lmGC(fastest~GPA,data=m111survey,graph=TRUE)
lmGC <-function(form,data=parent.frame(),graph=FALSE,check=FALSE)  {
  
  prsd <- ParseFormula(form)
  respname <- as.character(prsd$lhs)
  expname <- as.character(prsd$rhs)
  
  if (length(expname)>1) stop("Only one predictor variable permitted")
  
  resp <- simpleFind(varName=respname,data=data)
  exp <- simpleFind(varName=expname,data=data)
  
  
  if (!is(resp,"numeric")) stop("Response variable must be numerical")
  if (!is(exp,"numeric")) stop("Predictor variable must be numerical")
  
  #get the numbers from stats::lm
  df <- data.frame(exp,resp)
  names(df) <- c(expname,respname)
  form <- as.formula(paste0(respname,"~",expname))
  resultslm <- lm(form, data=df)
  results <- summary(resultslm)
  
  
  residse <- results$sigma
  
  
  #Collect what we need for our print function:
  results2 <- list(expname=expname,
                   respname=respname,
                   exp=exp,
                   resp=resp,
                   coefficients=results$coefficients,
                   r.squared=results$r.squared,
                   resid.sterr=residse,
                   graph=graph,
                   check=check,
                   mod=resultslm)
  
  class(results2) <- "GClm"
  return(results2)
  
}#end GClm

#' @title Diagnostic Plots for GC Linear Regression

#' @description Used by generic plot function
#' 
#' @rdname plot.GClm
#' @method plot GClm
#' @usage 
#' \S3method{plot}{GClm}(x,...)
#' @param x An object of class GClm
#' @param \ldots ignored
#' @return two diagmostic plots
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' SpeedModel <- lmGC(fastest~GPA,data=m111survey)
#' plot(SpeedModel)
plot.GClm <-function(x,...)  {
  
  GClm <- x
  mod <- GClm$mod
  residuals <- mod$residuals
  fitted.values <- mod$fitted.values
  
  p1 <- lattice::densityplot(~residuals,xlab="residuals",main="Residuals")
  p2 <- lattice::xyplot(residuals~fitted.values,xlab="predicted y values",
                        ylab="residuals",main="Residuals vs. Fits",pch=19,
                        panel=function(...){
                          lattice::panel.xyplot(...)
                          lattice::panel.abline(h=0)
                        })   
  print(p1,split=c(1,1,1,2), more=TRUE)
  print(p2,split=c(1,2,1,2))
  
}#end plot.GClm

#' @title Prediction Function for GC Linear Regression

#' @description Used by generic predict function
#' 
#' @rdname predict.GClm
#' @method predict GClm
#' @usage 
#' \S3method{predict}{GClm}(object,x,level=NULL,...)
#' @param object An object of class GClm
#' @param x value of the predictor variable
#' @param level desired level of prediction interval
#' @param \ldots ignored
#' @return numeric prediction
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #predict fastest speed driven, for person with GPA=3.0:
#' SpeedModel <- lmGC(fastest~GPA,data=m111survey)
#' predict(SpeedModel,x=3.0)
#' #include prediction interval:
#' predict(SpeedModel,x=3.0,level=0.95)
predict.GClm <-function(object,x,level=NULL,...)  {
  
  expname <- object$expname
  respname <- object$respname
  exp <- object$exp
  model <- object$mod
  
  if (!is.null(level)) {
    if (level <=0 || level >= 1) {
      stop("Level must be a number between 0 and 1")
    }
  }
  
  residse <- object$resid.sterr
  
  newdf <- data.frame(x)
  names(newdf) <- expname
  
  prediction1 <- predict(model,newdata=newdf,se.fit=TRUE)
  predVal <- prediction1$fit
  sepred <- sqrt(residse^2+(prediction1$se.fit)^2)
  
  cat(paste0("Predict ",respname," is about ",signif(predVal,4),
             ",\ngive or take ",signif(sepred,4)," or so for chance variation.\n\n"))
  
  if (!is.null(level)) {
    
    prediction2 <- suppressWarnings(predict(model,newdata=newdf,interval="prediction",level=level))
    lower <- prediction2[2]
    upper <- prediction2[3]
    cat(paste0(100*level,"%-prediction interval:\n"))
    int <- c(lower,upper)
    cat(sprintf("%-10s%-20s%-20s","","lower.bound","upper.bound"),"\n")
    cat(sprintf("%-10s%-20f%-20f","",int[1],int[2]),"\n\n")
    
  }
  
}#end predict.lmGC

#' @title Print Function for GC Linear Regression

#' @description Utility print function
#' @keywords internal
#' 
#' @rdname print.GClm
#' @method print GClm
#' @usage 
#' \S3method{print}{GClm}(x,...)
#' @param x an object of class GClm
#' @param \ldots ignored
#' @return graphical output and output to console
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
print.GClm <-function(x,...)  {
  GClm <- x
  coefs <- GClm$coefficients
  exp <- GClm$exp
  resp <- GClm$resp
  respname <- GClm$respname
  expname <- GClm$expname
  cat("\n") 
  
  cat("\tLinear Regression\n\n")
  cat("Correlation coefficient r = ",signif(cor(GClm$exp,GClm$resp,use="na.or.complete"),4),"\n\n")
  cat("Equation of Regression Line:\n\n")
  cat("\t",respname,"=",round(coefs[1],4),"+",round(coefs[2],4),"*",
      expname,"\n")
  cat("\n")
  
  cat("Residual Standard Error:\ts   =",round(GClm$resid.sterr,4),"\n")
  cat("R^2 (unadjusted):\t\tR^2 =",round(GClm$r.squared,4),"\n")
  
  
  #make data frame with complete cases to suppress warnings in ggplot RE missing data
  df <- data.frame(GClm$exp,GClm$resp)
  names(df) <- c(expname,respname)
  df <- df[complete.cases(df),]
  
  
  if (GClm$graph && !GClm$check) {
    
    title <- paste0("Scatterplot with linear fit")
    
    
    p1 <- ggplot2::ggplot(df,ggplot2::aes_string(x=expname,y=respname))+
      ggplot2::ggtitle(title)+
      ggplot2::geom_point()+
      ggplot2::stat_smooth(method = "lm",size = 1,se=FALSE)+
      ggplot2::xlab(expname)+ggplot2::ylab(respname)
    
    suppressWarnings(print(p1)) #suppress just in case Hadley has more friendly advice
    
  }
  
  
  if (GClm$check) {
    
    if (length(GClm$exp) < 1000) method <- "loess" else method <- "gam"
    
    title <- paste0("Checking the Model Fit\n(Model is blue; ",method,
                    " curve is red;\n95%-confidence band for curve included)")
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x=expname,y=respname))+
      ggplot2::ggtitle(title)+
      ggplot2::geom_point()+
      ggplot2::stat_smooth(method = "lm", size = 1,se=FALSE)+
      ggplot2::xlab(expname)+ggplot2::ylab(respname) + 
      ggplot2::stat_smooth(method=method,color="red",size=1,se=TRUE)
    
    
    suppressWarnings(print(p1))
    
  }
  
  
}#end print.GClm
