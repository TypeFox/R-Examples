#' @title Find the Regression Line

#' @description The regression minimizes the residual sum of squares (RSS).  In this game, the player chooses slope and y-intercept
#' of a line so as to approximate the regression line.  The moveable line is set initially as a horizontal line with height equal to
#' the mean of the y-coordinates of the scatterplot, so initially the residual sum of squares equals the total sum of squares (TSS).
#' The player's score is the sum of the number of turns taken and the difference between the current RSS
#' and the regression line's RSS (as a percentage of  TSS-RSS for regression line).  The aim is to lower one's score.
#' 
#' @rdname FindRegLine
#' @usage FindRegLine  
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note Requires package \code{manipulate}, available only in RStudio.
#' @examples
#' \dontrun{
#' if (require(manipulate)) FindRegLine()
#' }
FindRegLine <-
function()  {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  lowa <- -5
  higha <- 5
  lowb <- -2
  highb <- 2
  sigma <- 3
  ta <- round(runif(1,min=lowa,max=higha),1)
  tb <- round(runif(1,min=lowb,max=highb),1)
  n <- 10  #number of points
  x <- 1:n
  y <- ta+tb*x+rnorm(n,mean=0,sd=sigma)
  
  #SS for the regression line:
  mod <- lm(y~x)
  ess <- sum((resid(mod))^2)
  
  #determine nice limits for plot (and a slider):
  reg.slope <- coef(mod)[2]
  reg.int <- coef(mod)[1]
  
  #Find range of y-intercepts of lines through
  #points on scatterplot having slope = reg.slope
  int.min <- min(y-reg.slope*x)
  int.max <- max(y-reg.slope*x)
  int.band <- (int.max-int.min)/2
  
  #Expand this range, and make sure it includes 0:
  int.mid <- (int.max+int.min)/2
  lowa.slider <- floor(min(c(int.mid-1.2*int.band,-1,min(y)-1)))
  higha.slider <- ceiling(max(c(int.mid+1.2*int.band,1,max(y)+1)))
  
  #plot limits reflect this range, too:
  ymin <- lowa.slider
  ymax <- higha.slider
  y.mean <- mean(y)
  
  #SS for the line initially placed (a=0,b=0):
  total.ss <- sum((y-mean(y))^2)
  your.ss.init <- total.ss
  #You start at line with slope 0, intercept = mean(y)
  
  #Stuff to help keep score:
  turns <- 0
  aprev <- y.mean
  bprev <- 0
  a <- y.mean #Probably not needed
  b <- 0  #probably not needed
  
  #Score is number of turns taken
  #plus a relative measure of your closeness to optimal:
  close <- 100*(your.ss.init-ess)/(total.ss-ess)
  score <- turns+close
  #Initial score will alwauys be 100, because you
  #start at line with slope 0, intercept = mean(y)
  
  manipulate(
    a=slider(lowa.slider,higha.slider,step=0.01,initial=y.mean,label="Intercept"),
    b=slider(2*lowb,2*highb,step=0.01,initial=0,label="Slope"),
    displayscore=checkbox(FALSE,"Show my score"),
    enditall=checkbox(FALSE,"I quit -- show the reg line!"),
{
 plot(x,y,pch=16,col="blue",ylim=c(ymin,ymax),
      xlim=c(0,n))
 points(0,0,cex=0.8,pch=16,col="green")
 abline(a,b)
 abline(0,0,lty=2,col="green")
 lines(x=c(0,0),y=c(ymin,ymax),lty=2,col="green")
 your.y <- a+b*x
 fits <- fitted(mod)
 your.ss <- sum((y-your.y)^2)
 for(i in 1:n)  {
   lines(x=c(x[i],x[i]),y=c(your.y[i],y[i]))
 }
 mtext(paste("Your Sum of Squares =",
             round(your.ss,2)),line=2)
 mtext(paste("Minimal Sum of Squares =",
             round(ess,2)),line=0.5)
 
 if (!identical(a,aprev) || !identical(b,bprev)){
   #Update:
   turns <<- turns+1
   close <<- 100*(your.ss-ess)/(total.ss-ess)
   score <<- turns+close
   aprev <<- a
   bprev <<- b
 }
 
 if (displayscore)  {
   cat("\n")
   cat("You have taken",turns,"turns.\n")
   cat("Your current closeness (0 is perfect) is:",close,"\n")
   cat("So your current score is",score,"\n")
   cat("Think carefully about your next move!\n")
 }
 
 if(enditall)  {
   coefs <- coef(mod)
   abline(coefs,col="red",lwd=3)
   coefs <- round(coefs,2)
   cat("Equation of the regression line is:\n")
   cat(paste("y = ",coefs[1]," + ",coefs[2],"x\n",sep=""))
   cat("\n")
   cat("Your final score is",score,"\n")
   return(cat("Thanks for playing!"))
 }
 
}#end manipulate actions
  )#end manipulate
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("displayscore","enditall"))
