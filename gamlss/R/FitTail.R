# library(gamlss)
# library(gamlss.tr)
# those are function for checking the tails of a distribution
# created by Mikis Stasinopoulos Bob Rigby and Vlasios Voudouris Dec 2011 
# TO DO
# i) logSurv : which parameters we should save? and what they represent? save? 
# ii)
# it contains
#-------------------------------------------------------------------------------
# log log Survival plots
# there are three functions for that plus one which connect them al
# i)   loglogSurv1()  : for type I tails
# ii)  loglogSurv2()  : for type II tails
# iii) loglogSurv3()  : for type III tails
# iv)  loglogSurv     : combines the above functions and selects the best
#-------------------------------------------------------------------------------
# log Survival plot
# v)    logSurv() : plots the empirical suvuival function log(1-ecdf) against log(y) OK
#-------------------------------------------------------------------------------
# vi)     fitTail : fits a truncated gamlss.family distribution to the tail of the data
# vii) fitTailAll : fits a (Hill type) series of Fit using the fitTail function 
#------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# TYPE I
#-------------------------------------------------------------------------------
loglogSurv1 <- function(y, 
            percentage = 10, 
               howmany = NULL, 
                  type = c("right", "left"), 
                  plot = TRUE, 
                 print = TRUE,
                  save = FALSE)
{
#-----------------
#  local function  
  tailFun <- function(y, percentage, howmany, type )
 {
        ly <- length(y)
   howmany <- if (is.null(howmany)) floor(ly*(percentage/100)) else howmany
         Y <- if (type=="right") tail(y[order(y)], howmany)	
              else abs(head(y[order(y)], howmany))
         Y
 }
# body of the function starts here  
  type <- match.arg(type)
   cdF <- ecdf(y)
     Y <- tailFun(y, percentage=percentage, howmany=howmany, type=type )
    mY <- min(Y)
if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1")
# newY <- log(-log(1-cdF(Y)+0.0001))
  newY <- log(-log(1-cdF(Y)+(1/(2*length(y)))))# different versions to define the ecdf
 #newY <- log(-log(1-cdF(Y)*(length(y)/(length(y)+1))))
     x <- log(log(Y))
    m1 <- gamlss(newY ~ x, trace=FALSE)
   ess <- sum((newY-fitted(m1))^2)
 Xlab1 <- paste("log(log(", deparse(substitute(y)), "))", sep="")
  Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")  
if (print)
  {
   cat("coefficients",  coef(m1), "\n")
  cat("error sum of squares", ess, "\n")
  }
if (plot)
  {
  plot(newY~x, xlab=Xlab1, ylab=Ylab)
  lines(fitted(m1)~x, col="red")  
  if (is.null(howmany)) 
    title(paste("Log Log Survival plot (Type I) for", paste(paste(percentage,"%",sep=""), "of the", type, "tail", "of",  deparse(substitute(y)))))
  else
 title(paste("Log Log Survival plot (Type I) for", paste( howmany,"observations", "of the", type, "tail", "of",  deparse(substitute(y)))))
  }
if (save)
  {
return(m1)  
  }  
}
#--------------------------------------------------------------------------------
# Type II
#--------------------------------------------------------------------------------
loglogSurv2 <- function(y, 
            percentage = 10, 
               howmany = NULL, 
                  type = c("right", "left"), 
                  plot = TRUE, 
                 print = TRUE,
                  save = FALSE)
{
#-----------------
#  local function  
  tailFun <- function(y, percentage, howmany, type )
 {
        ly <- length(y)
   howmany <- if (is.null(howmany)) floor(ly*(percentage/100)) else howmany
         Y <- if (type=="right") abs(tail(y[order(y)], howmany))  
              else abs(head(y[order(y)], howmany))
         Y
 }
# body of the function starts here
require(gamlss)
  type <- match.arg(type)
   cdF <- ecdf(y)
     Y <- tailFun(y, percentage=percentage, howmany=howmany, type=type )
    mY <- min(Y)
if (mY<0) stop("For the method for Type II to work the minimum value of the tail must be greater than 0")
# newY <- log(-log(1-cdF(Y)+0.0001))
  newY <- log(-log(1-cdF(Y)+(1/(2*length(y)))))# ??????????
# newY <- log(-log(1-cdF(Y)*(length(y)/(length(y)+1))))
     x <- log(Y)
    m2 <- gamlss(newY ~ x, trace=FALSE)
   ess <- sum((newY-fitted(m2))^2)
 Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
  Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")
if (print)
   {
   cat("coefficients",  coef(m2), "\n")
   cat("error sum of squares", ess, "\n")
   }
if (plot)
   {
  plot(newY~x, xlab=Xlab1, ylab=Ylab)
  lines(fitted(m2)~x, col="red")  
  if (is.null(howmany)) 
    title(paste("Log Log Survival plot (Type II) for", paste(paste(percentage,"%",sep=""), "of the", type, "tail", "of",  deparse(substitute(y)))))
  else
 title(paste("Log Log Survival plot (Type II) for", paste( howmany,"observations", "of the", type, "tail", "of",  deparse(substitute(y)))))
   }
if (save)
  {
return(m2)  
  }  
}
#--------------------------------------------------------------------------------
# Type III
#--------------------------------------------------------------------------------
loglogSurv3 <- function(y, 
            percentage = 10, 
               howmany = NULL, 
                  type = c("right", "left"), 
                  plot = TRUE, 
                 print = TRUE,
                  save = FALSE)
{
#-----------------
#  local function  
  tailFun <- function(y, percentage, howmany, type )
 {
        ly <- length(y)
   howmany <- if (is.null(howmany)) floor(ly*(percentage/100)) else howmany
         Y <- if (type=="right") abs(tail(y[order(y)], howmany))  
              else abs(head(y[order(y)], howmany))
         Y
 }
# body of the function starts here
require(gamlss)
  type <- match.arg(type)
   cdF <- ecdf(y)
     Y <- tailFun(y, percentage=percentage, howmany=howmany, type=type )
    mY <- min(Y)
# newY <- log(-log(1-cdF(Y)+0.0001))
  newY <- log(-log(1-cdF(Y)+(1/(2*length(y)))))# ??????????
# newY <- log(-log(1-cdF(Y)*(length(y)/(length(y)+1))))
     x <- Y
    m3 <- gamlss(newY ~ x, trace=FALSE)
   ess <- sum((newY-fitted(m3))^2)
  Xlab <- paste(deparse(substitute(y)))
  Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")
if (print)
  {
   cat("coefficients",  coef(m3), "\n")
  cat("error sum of squares", ess, "\n")
  }
if (plot)
  {
  plot(newY~x, xlab=Xlab, ylab=Ylab)
  lines(fitted(m3)~x, col="red")  
  if (is.null(howmany)) 
    title(paste("Log Log Survival plot (Type III) for", paste(paste(percentage,"%",sep=""), "of the", type, "tail", "of",  deparse(substitute(y)))))
  else
 title(paste("Log Log Survival plot (Type III) for", paste( howmany,"observations", "of the", type, "tail", "of",  deparse(substitute(y)))))
  }
if (save)
  {
return(m3)  
  }  
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# Select the best from type I II and III
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
loglogSurv <- function(y, 
            percentage = 10, 
               howmany = NULL, 
                  type = c("right", "left"), 
                  plot = TRUE, 
                 print = TRUE,
                  save = FALSE)
{
#-----------------
#  local function  
  tailFun <- function(y, percentage, howmany, type )
 {
        ly <- length(y)
   howmany <- if (is.null(howmany)) floor(ly*(percentage/100)) else howmany
         Y <- if (type=="right") abs(tail(y[order(y)], howmany))  
              else abs(head(y[order(y)], howmany))
         Y
 }
# body of the function starts here
  type <- match.arg(type)
   cdF <- ecdf(y)
     Y <- tailFun(y, percentage=percentage, howmany=howmany, type=type )
    mY <- min(Y)
   if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1")
 # newY <- log(-log(1-cdF(Y)+0.0001))
  newY <-  log(-log(1-cdF(Y)+(1/(2*length(y)))))# ??????????
 #  newY <-  log(-log(1-cdF(Y)*(length(y)/(length(y)+1))))
    x1 <- log(log(Y))
    m1 <- gamlss(newY ~ x1, trace=FALSE)
  ess1 <- sum((newY-fitted(m1))^2)
    x2 <- log(Y)
    m2 <- gamlss(newY ~ x2, trace=FALSE)
  ess2 <- sum((newY-fitted(m2))^2)
    x3 <- Y
    m3 <- gamlss(newY ~ x3, trace=FALSE)
   ess3 <- sum((newY-fitted(m3))^2) 
    ess <- c(ess1, ess2, ess3)
   num <- which.min(ess)
if (print)
 {
matcoef <- rbind(coef(m1), coef(m2), coef(m3))
matcoef <- cbind(matcoef, ess)
 dimnames(matcoef) <- list(c("type I", "type II", "type III"), c(" Intercept",
" slope", " Error SS"))
cat("Linear regression coefficients", "\n")
printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
matk <- matcoef[,-3]
matk[,1] <-exp(matk[,1])
cat("Estimates for parameters k", "\n")
printCoefmat(matk, digits = 3, signif.stars = TRUE)
 }
if (plot)
 {
  Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")
  switch(num, 
         { Xlab1 <- paste("log(log(", deparse(substitute(y)), "))", sep="")
           plot(newY~x1, xlab=Xlab1, ylab=Ylab)
           lines(fitted(m1)~x1, col="red")  
           if (is.null(howmany)) 
             title(paste("Log Log Survival plot (Type I) for", paste(paste(percentage,"%",sep=""), "of the", type, "tail", "of",  deparse(substitute(y)))))
           else
             title(paste("Log Log Survival plot (Type I) for", paste( howmany,"observations", "of the", type, "tail", "of",  deparse(substitute(y)))))},
         {
           Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
           plot(newY~x2, xlab=Xlab1, ylab=Ylab)
           lines(fitted(m2)~x2, col="red")  
           if (is.null(howmany)) 
             title(paste("Log Log Survival plot (Type II) for", paste(paste(percentage,"%",sep=""), "of the", type, "tail", "of",  deparse(substitute(y)))))
           else
             title(paste("Log Log Survival plot (Type II) for", paste( howmany,"observations", "of the", type, "tail", "of",  deparse(substitute(y))))) 
         },
         {
           Xlab1 <- paste(deparse(substitute(y)), sep="")
           plot(newY~x3, xlab=Xlab1, ylab=Ylab)
           lines(fitted(m3)~x3, col="red")  
           if (is.null(howmany)) 
             title(paste("Log Log Survival plot (Type III) for", paste(paste(percentage,"%",sep=""), "of the", type, "tail", "of",  deparse(substitute(y)))))
           else
             title(paste("Log Log Survival plot (Type III) for", paste( howmany,"observations", "of the", type, "tail", "of",  deparse(substitute(y))))) 
         }
         )   
 }
if (save)
 {
model <- switch(num, m1,m2,m3)   
return(model) 
 }
}

#--------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this plots the empirical log(1-ecdf) agains log(y)
# The complementary cumulative distribution function (CCDF) plot
# fits also a linear and a quadraitic fit to the resulting plot yo help interpetation 
logSurv <- function(y, 
            percentage = 10, 
               howmany = NULL, 
                  type = c("right", "left"), 
                  plot = TRUE, 
                 print = TRUE,
                  save = FALSE)
{
#-----------------
#  local function  
  tailFun <- function(y, percentage, howmany, type )
 {
        ly <- length(y)
   howmany <- if (is.null(howmany)) floor(ly*(percentage/100)) else howmany
         Y <- if (type=="right") abs(tail(y[order(y)], howmany))  
              else abs(head(y[order(y)], howmany))
         Y
 }
#----------------
# body of the function starts here
#  require(gamlss)
   type <- match.arg(type)
     #y <- as.numeric(y)
    cdF <- ecdf(y) # get ecdf
      Y <- tailFun(y, percentage=percentage, howmany=howmany, type=type ) # ger %
     mY <- min(Y)
   newY <- log(1-cdF(Y)+(1/(2*length(y))))# 
 # newY <-  log(1-cdF(Y)*(length(y)/(length(y)+1)))
     m1 <- gamlss(newY ~ log(Y), trace=FALSE) # fitting k1=1 model
 #   m2 <- gamlss(newY ~ I(log(Y)^2), trace=FALSE)
     m3 <- gamlss(newY ~ log(Y)+I(log(Y)^2), trace=FALSE)  # fitting k1=2 model
  newY2 <- log(-log(1-cdF(Y)+0.0001))
#     x <- log(Y)
     m4 <- gamlss(newY2 ~ log(Y), trace=FALSE) # fitting k
#  fv4  <- fitted(m4)
   fv4  <- -exp(fitted(m4))
   Xlab <- paste("log(", deparse(substitute(y)), ")", sep="")
   Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
if (plot)
 {
     plot(newY~log(Y), xlab=Xlab, ylab=Ylab)
    lines(fitted(m1)~log(Y), col="darkgreen", lty=3, lwd=2)
#    lines(fitted(m2)~log(Y), col=3, lty=2)
    lines(fitted(m3)~log(Y), col="blue", lty=2 , lwd=2)
    lines(fv4~log(Y), col="red", lty=1, lwd=2 )
   if (is.null(howmany)) title(paste(paste(percentage,"%",sep=""), "of the", type, "tail", "of", deparse(substitute(y))))
   else title(paste(howmany,"observations", "of the", type, "tail", "of", deparse(substitute(y))))
  legend("bottomleft", legend=c("linear", "quadratic", "exponential"), col=c("darkgreen", "blue", "red"), lty=c(3,2,1), lwd=2 )
 }
if (save)
{
  invisible(rbind(coef(m1), coef(m3)[-2], coef(m4))) 
}
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
