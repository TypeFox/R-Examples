#' @title Polynomial Regression

#' @description Regression analysis (one numerical predictor variable) with simplified output.
#'   Wrapper function for \code{lm} in package \code{stats}.
#' 
#' @rdname polyfitGC
#' @usage polyfitGC(form,data=parent.frame(),degree=2,graph=TRUE,check=FALSE)
#' @param form formula of form y~x, both variables numeric
#' @param data dataframe supplying y and x above.  If one or more of the variables is not in data, then
#' they will be searched for in the parent environment.
#' @param degree desired degree of polynomial (for degree 1 use lmgC)
#' @param graph Produce scatterplot with fitted ploynomial.
#' @param check Asks to produce a lowess or gam curve with approximate 95%-confidence band.  If the
#' fitted line wanders outside the band, then perhaps a linear fit is not appropriate.
#' @return A list of class "polyGC".  Elements that may be queried include
#' "s" (residual standard error) and "R^2" (unadjusted).
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #To study the relationship between two numerical variables:
#' polyfitGC(mpg~wt,data=mtcars,degree=2,graph=TRUE)
#' #check the second-fdegree fit:
#' polyfitGC(mpg~wt,data=mtcars,degree=2,check=TRUE)
polyfitGC <-function(form,data=parent.frame(),degree=2,graph=TRUE,check=FALSE)  {
  
  if (degree==1) stop("For linear fit, use lmGC().")
  
  prsd <- ParseFormula(form)
  respname <- as.character(prsd$lhs)
  expname <- as.character(prsd$rhs)
  
  if (length(expname)>1) stop("Only one predictor variable permitted")
  
  resp <- simpleFind(varName=respname,data=data)
  exp <- simpleFind(varName=expname,data=data)
  
  
  if (!is(resp,"numeric")) stop("Response variable must be numerical")
  if (!is(exp,"numeric")) stop("Predictor variable must be numerical")
  
  #get the numbers from stats::lm
  
  # first center the x-values to minimize numerical issues
  
  meanX <- mean(exp,na.rm=TRUE)
  sdX <- sd(exp,na.rm=TRUE)
  expCent <- (exp - meanX)/sdX
  
  df <- data.frame(expCent,resp)
  names(df) <- c(expname,respname)
  form <- as.formula(paste0(respname,"~ poly(",expname,",",degree,", raw=TRUE)"))
  polyMod <- lm(form, data=df)
  results <- summary(polyMod)
  

  residse <- results$sigma
  
  #Collect what we need for our print function:
  results2 <- list(expname=expname,
                   respname=respname,
                   exp=exp,
                   resp=resp,
                   residuals=polyMod$residuals,
                   coefficients=results$coefficients,
                   r.squared=results$r.squared,
                   resid.sterr=residse,
                   graph=graph,
                   check=check,
                   degree=degree,
                   meanX=meanX,
                   sdX=sdX,
                   mod=polyMod)
  
  class(results2) <- "polyGC"
  return(results2)
  
}#end polyfitGC


#' @title Prediction Function for GC Polynomial Regression

#' @description Used by generic predict function
#' 
#' @rdname predict.polyGC
#' @method predict polyGC
#' @usage 
#' \S3method{predict}{polyGC}(object,x,level=NULL,...)
#' @param object An object of class polyGC
#' @param x value of the predictor variable
#' @param level desired level of prediction interval
#' @param \ldots ignored
#' @return numeric prediction
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #predict mpg for a car weighing 3 tons:
#' mpgModel <- polyfitGC(mpg~wt,data=mtcars,degree=2)
#' predict(mpgModel,x=3.0)
#' #include prediction interval:
#' predict(mpgModel,x=3.0,level=0.95)
predict.polyGC <-function(object,x,level=NULL,...)  {
  
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
  
  xCent <- (x - object$meanX)/object$sdX
  
  newdf <- data.frame(xCent)
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
  
}#end predict.polyGC

#' @title Diagnostic Plots for GC Polynomial Regression

#' @description Used by generic plot function
#' 
#' @rdname plot.polyGC
#' @method plot polyGC
#' @usage 
#' \S3method{plot}{polyGC}(x,...)
#' @param x An object of class polyGC
#' @param \ldots ignored
#' @return two diagmostic plots
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' mpgModel <- polyfitGC(mpg~wt,data=mtcars)
#' plot(mpgModel)
plot.polyGC <-function(x,...)  {
  
  polyGC <- x
  mod <- polyGC$mod
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
  
}#end plot.polyGC

#' @title Print Function for GC Polynomial Regression

#' @description Utility print function
#' @keywords internal
#' 
#' @rdname print.polyGC
#' @method print polyGC
#' @usage 
#' \S3method{print}{polyGC}(x,...)
#' @param x an object of class polyGC
#' @param \ldots ignored
#' @return graphical output and output to console
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
print.polyGC <-function(x,...)  {
  
  length_out <- 1000 #for plotting fitted curve
  
  polyGC <- x
  degree <- polyGC$degree
  
  respname <- polyGC$respname
  expname <- polyGC$expname
  
  #          if (degree==2) {
  #            cat("\tQuadratic Regression\n\n")
  #            cat("Equation of Fitted Parabola:\n\n")
  #            cat(paste0("\t",respname,"=",signif(coefs[1],4)," + ",
  #                signif(coefs[2],4),"*",expname," + " ,
  #                signif(coefs[3],4),"*",expname,"^2","\n"))
  #            cat("\n")
  #          }
  # 
  #          if (degree==3) {
  #            cat("\tCubic Regression\n\n")
  #            cat("Equation of Fitted Cubic:\n\n")
  #            cat(paste0("\t",respname," = ",signif(coefs[1],4)," + ",
  #                     signif(coefs[2],4),"*",expname," + ",
  #                     signif(coefs[3],4),"*",expname,"^2"," + ",
  #                     signif(coefs[4],4),"*",expname,"^3","\n"))
  #            cat("\n")
  #          }
  #          
  #          if (degree > 3) {
  #            cat("\tPolynomial Regression\n\n")
  #            cat("Equation of Fitted Polynomial:\n\n")
  #            eq <- paste0("\t",respname," = ",signif(coefs[1],4)," + ",
  #                           signif(coefs[2],4),"*",expname,"\n")
  #            for (i in 3:(degree+1)) {
  #              eq <- c(eq,paste0("\t\t\t\t + ",
  #                               signif(coefs[i],4),"*",expname,"^",i-1,"\n"))
  #            }
  #            cat(eq)
  #            cat("\n")
  #          }
  
  cat(paste0("Polynomial Regression, Degree = ",degree,"\n\n"))
  cat("Residual Standard Error:\ts   =",round(polyGC$resid.sterr,4),"\n")
  cat("R^2 (unadjusted):\t\tR^2 =",round(polyGC$r.squared,4),"\n")
  
  
  #make data frame with complete cases to suppress warnings in ggplot RE missing data
  df <- data.frame(polyGC$exp,polyGC$resp)
  names(df) <- c(expname,respname)
  df <- df[complete.cases(df),]
  
  exp <- df[,expname]
  resp <- df[,respname]
  
  xFill <- seq(min(exp),max(exp),length.out=length_out)
  meanX <- polyGC$meanX
  sdX <- polyGC$sdX
  
  xFillCent <- (xFill - meanX)/sdX
  newdf <- data.frame(xFillCent)
  names(newdf) <- expname
  
  mod <- polyGC$mod
  fitsFill <- predict(mod,newdata=newdf)
  
  predfr <- data.frame(xFill,fitsFill)
  names(predfr) <- c(expname,respname)
  
  
  if (polyGC$graph && !polyGC$check) {
    
    
    
    
    if (degree == 2) {
      title <- paste0("Scatterplot with quadratic fit")
    }
    
    if (degree == 3) {
      title <- paste0("Scatterplot with cubic fit")
    }
    
    if (degree >= 4) {
      title <- paste0("Scatterplot with degree-",degree," fit")
    }
    
    p1 <- ggplot2::ggplot(df,ggplot2::aes_string(x=expname,y=respname))+
      ggplot2::ggtitle(title)+
      ggplot2::geom_point()+
      ggplot2::geom_point(data=predfr,ggplot2::aes_string(x=expname,y=respname),col="blue",size=1)+
      ggplot2::xlab(expname)+ggplot2::ylab(respname)
    
    suppressWarnings(print(p1))
    
  }
  
  
  if (polyGC$check) {
    
    if (length(polyGC$exp) < 1000) method <- "loess" else method <- "gam"
    
    title <- paste0("Checking the Model Fit\n(Model is blue; ",method,
                    " curve is red;\n95%-confidence band for ",method," curve is included)")
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x=expname,y=respname))+
      ggplot2::ggtitle(title)+
      ggplot2::geom_point()+
      ggplot2::stat_smooth(method=method,color="red",size=1,se=TRUE)+
      ggplot2::geom_point(data=predfr,ggplot2::aes_string(x=expname,y=respname),col="blue",size=1)+
      ggplot2::xlab(expname)+ggplot2::ylab(respname)
    
    
    suppressWarnings(print(p1))
    
  }
  
  
}#end print.polyGC
