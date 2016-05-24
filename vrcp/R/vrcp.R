##' Estimation for Varying Regression Segments and Change Point in Heteroscedastic Data
##' 
##' @description Estimation of two segments and a change point in 2-segment regression models with varying variances and varying types of regression segments, with or without a smoothness constraint at the change point.
##' @param dataset either a data frame, or a matrix, containing 2 columns, where the first column contains covariate values and the second column contains response values.
##' @param lo,hi lower and upper bounds for the x value of the change point, to be set by a user.
##' @param smooth smoothness constraint of the regression function at the change point. 
##' Default constraint is "c0," continuous at the change point. "c1" indicates that the first derivative at the change point is continuous, while "u" indicates that there is constraint at the change point.
#' @param segment1,segment2 regression model used to compute parameters in segment1(or segment2) with additive Gaussian errors. Currently allowable models are:
##'
##' L: Linear model. y ~ a0 + a1 * x 
##'
##' Q: Quadratic model. y ~ a0 + a1 * x + a2 * x^2
##'
##' Exp: Linearizable Exponential model. y ~ a0 + a1 * exp(x)
##'
##' Log: Linearizable Logarithm model. y ~ a0 + a1 * log(x) 
##'
##' NLExp: Nonlinearizable Exponential model.
##' y ~ a0 + a1 * exp(a2 * (x - k)), where k is the change point in x.
##' 
##' These 5 models lead to the following 15 allowable combinations of varying types of segments, reflecting reasonable models we have seen from Degredation Science: 
##' 
##' "L"-"L", "L"-"Q", "L"-"Exp", "L"-"Log", "L"-"NLExp",
##' 
##' "Q"-"L", "Q"-"Q", "Q"-"Exp", "Q"-"Log", and "Q"-"NLExp",
##' 
##' "Exp"-"L", "Exp"-"Exp",
##' 
##' "Log"-"L", "Log"-"Log",
##' 
##' and "NLExp"-"L".
##' 
##' @param variance variance type of the data set, either "Common" that requires the variances at two segments to be the same, or "Diff" that does NOT require them to be the same.
##' @param spline a "TRUE" or "FALSE" logical argument, where "TRUE" shows a B-spline fit with the knot at the change point, as an extra option only available for segment combinations of "L"-"L" and "Q"-"Q". The default is "FALSE".
##' @param start a named numeric vector or a logic status "FALSE". Default is "start=FALSE", which will compute initial values automatically based on data.
##' "start=a vector" specifies initial values of vector parameters, typically in the case with a nonlinearizable segment as specified below, where a0, a1, a2 are regression parameters for segment1, and b0, b1, b2 are regression parameters for segment2, etc:
##' 
##' "L"-"NLExp": 3 initial values for (a0, a1; b2) may be specified.
##'
##' "Q"-"NLExp": 4 initial values for (a0, a1, a2; b2) may be specified.
##' 
##' "NLExp"-"L": 3 initial values for (b0, b1; a2) may be specified.
##'
##' @export
##' @importFrom 'ggplot2' 'aes' 'scale_colour_manual'
##' @importFrom 'stats' 'lm' 'nls' 'predict'
##' @return
##' maxloglik: maximum log-likelihood value.
##' 
##' sigma2: estimated variance(s) for two segments
##'
##' coe: coefficients for two regression segments,
##' beta = (a0,a1,a2,b0,b1,b2). No a2, b2 output for linear/linearizable segment.
##' 
##' changepoint: change point in x value.
##' 
##' @references
##' Stephen J. Ganocy and Jiayang Sun (2015), 
##' "Heteroscedastic Change Point Analysis and Application to Footprint Data",
##' J of Data Science, v.13.
##' @examples 
##' library(ggplot2)
##' 
##' # Test the vrcp() using simulated data sets
##'
##' # Example 1: L-L model with "c0", continuity at change point and common variance
##' # Simulate the data
##' x1<-seq(0,2,by=0.05)
##' x2<-seq(2.05,5,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 2+0.5*x1
##' yt2 <- -1+2*x2
##'
##' # Add noises
##' y1<-yt1+rnorm(length(x1),0,0.25)
##' y2<- yt2+rnorm(length(x2),0,0.25)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##' 
##' # It looks like a L-L regression with a change point between 1.5 and 2.5
##' # Fit with vrcp with L-L segments and "c0" constraint
##' ans <- vrcp(z,1.5,2.5,"c0","L","L","Common", spline = "TRUE") # Fit with common variance
##' ans
##' 
##' # The fitted L-L regression and spline are superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt), color = c("blue"), size=1) +
##' ggtitle("LL-c0-com model: Estimates vs. true model (in blue)")
##'
##' ans <- vrcp(z,1.5,2.5,"c0","L","L","Diff",spline = "TRUE") # Fit with different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt), color = "blue", size=1) +
##' ggtitle("LL-c0-diff model: Estimates vs. true model (in blue)") # compare
##'
##' 
##' \dontrun{
##' # Example 2: L-Log model with "c1" change point and common variance.
##' # Simulate the data
##' x1<-seq(0.05,2.05,by=0.05)
##' x2<-seq(2.1,5.05,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 3+1*x1
##' yt2 <- 3.61+2*log(x2)
##'
##' # Add noises
##' y1<- yt1+rnorm(length(x1),0,0.5)
##' y2<- yt2+rnorm(length(x2),0,0.5)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)<-c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##'
##' # It looks like a L-Log regression with a change point between 1.9 and 2.2
##' # Fit with vrcp with specification of L-Log segments and "c1" options with 
##' # and without common variance restriction
##' ans <- vrcp(z,1.9,2.2,"c1","L","Log","Common")
##' ans
##' 
##' # The fitted L-Log regression is superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LLog-c1-com model: Estimate (in magenta) vs. true model") # Fit with common variance
##'
##' ans <- vrcp(z,1.9,2.5,"c1","L","Log","Diff") 
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LLog-c1-diff model: Estimate (in magenta) vs. true model") # Fit with different variance
##' 
##' # both fits look good
##'
##' # Check what would look like with misspecification of smoothness at change point
##' ans <- vrcp(z,1.9,2.2,"c0","L","Log","Common") # Fit with common variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LLog-c0-com fit to LLog-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.5,2.5,"c0","L","Log","Diff") # Fit with different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LLog-c0-diff fit to LLog-c1-com model: Estimate (in magenta) vs. true model")
##'
##' # both look lack of fit, especially at the change point. 
##' # Hence, the correct specification of model is important 
##'
##'
##' # Example 3: Log-L - Simulated data set is "c1", smooth.
##' # Simulate the data
##' x1<-seq(2,4,by=0.05)
##' x2<-seq(4,7,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 1.6+0.5*log(x1)
##' yt2 <- 1.89+0.1*x2
##'
##' # Add noises
##' y1<- yt1+rnorm(length(x1),0,0.1)
##' y2<- yt2+rnorm(length(x2),0,0.1)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##' 
##' # It looks like a Log-L regression with a change point between 3.9 and 4.5
##' # Fit with vrcp with specification of Log-L segments and "c1" options with 
##' # and without common variance restriction
##' ans <- vrcp(z,3.9,4.5,"c1","Log","L","Common") # Fit with common variance
##' ans
##'
##' # The fitted Log-L regression is superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") + 
##' ggtitle("LogL-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,3.8,4.2,"c1","Log","L","Diff") # Fit with different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LogL-c1-diff model: Estimate (in magenta) vs. true model")
##' 
##' # results look similar, not bad.
##'
##' # Fit with Log-L segments and "c0" options with and without common variance restriction
##' ans <- vrcp(z,3.5,4.5,"c0","Log","L","Common") # Common variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LogL-c0-com fit to LLog-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,3.5,4.5,"c0","Log","L","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LogL-c0-diff fit to LLog-c1-com model: Estimate (in magenta) vs. true model")
##' 
##' # Little worse than the one with c1 constraint
##'
##' # Fit with Log-L segments and u" options with and without common variance restriction
##' ans <- vrcp(z,3.5,4.5,"u","Log","L","Common") # Common variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LLog-u-com fit to LLog-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,3.5,4.5,"u","Log","L","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LLog-u-diff fit to LLog-c1-diff model: Estimate (in magenta) vs. true model")
##'
##' # Clearly shows lack of fit at the change point.   
##' # Again, the correct specification of the model, or use of available information is important.
##'
##'
##' # Example 4: QL-c1-com model, fitted by Q-L and Exp-L models, 
##' # with and without a common variance constraint, respectively.
##'
##' # Simulate Q-L data
##' x1<-seq(0,2,by=0.05)
##' x2<-seq(2,5,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 2+2*x1+2*x1^2
##' yt2 <- -6+10*x2
##'
##' # Add noises
##' y1<- yt1+rnorm(length(x1),0,3)
##' y2<- yt2+rnorm(length(x2),0,3)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##'
##' # It looks like a Q-L regression with a change point between 1.8 and 2.5
##' # Fit with vrcp with specification of Q-L segments and "c1" options with 
##' # and without common variance restriction
##' ans <- vrcp(z,1.8,2.5,"c1","Q","L","Common")  # Common variance
##' ans
##'
##' # The fitted Q-L regression is superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("QL-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.8,2.5,"c1","Q","L","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("QL-c1-diff model: Estimate (in magenta) vs. true model")
##'
##' # Fit with vrcp with specification of Exp-L segments and "c1" options with 
##' # and without common variance restriction
##'
##' ans <- vrcp(z,1.5,2.5,"c1","Exp","L","Common") # Common variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpL-c1-com fit to QL-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.5,2.5,"c1","Exp","L","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpL-c1-diff fit to QL-c1-com model: Estimate (in magenta) vs. true model")
##'
##' # Fit with vrcp with specification of Exp-L segments and "c0" options with 
##' # and without common variance restriction
##'
##' ans <- vrcp(z,1.5,2.5,"c0","Exp","L","Common") # Common variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpL-c0-com fit to QL-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.5,2.5,"c0","Exp","L","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpL-c0-com fit to QL-c1-com model: Estimate (in magenta) vs. true model")
##'
##' # Fit with vrcp with specification of Exp-L segments and "u" options with 
##' # and without common variance restriction
##'
##' ans <- vrcp(z,1.5,2.5,"u","Exp","L","Common") # Common variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpL-u-com fit to QL-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.5,2.5,"u","Exp","L","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpL-u-com fit to QL-c1-com model: Estimate (in magenta) vs. true model")
##' 
##' # Exp-L fits surprisingly well in this case.
##' 
##'
##' # Example 5: Exp-Exp with "c0" change point and common variance. - No option of smoothness
##' # Simulate the data
##' x1<-seq(0,2,by=0.05)
##' x2<-seq(2.05,5,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 0.916+2*exp(x1)
##' yt2 <- 12+0.5*exp(x2)
##'
##' # Add noises
##' y1<-yt1+rnorm(length(x1),0,5)
##' y2<-yt2+rnorm(length(x2),0,5)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##'
##' # It looks like a Exp-Exp regression with a change point between 1.5 and 2.5
##' # Fit with vrcp with specification of Exp-Exp segments and "c0" options with 
##' # and without common variance restriction
##'
##' ans <- vrcp(z,1.5,2.5,"c0","Exp","Exp","Common") # Common variance ## simulation of smooth
##' ans
##'
##' # The fitted Exp-Exp regression is superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpExp-c0-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.5,2.5,"c0","Exp","Exp","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpExp-c0-diff model: Estimate (in magenta) vs. true model")
##'
##' # Fit with vrcp with specification of Exp-Exp segments and "u" options with 
##' # and without common variance restriction
##'
##' ans <- vrcp(z,1.5,2.2,"u","Exp","Exp","Common") # Common variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpExp-u-com fit to ExpExp-c0-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.5,2.5,"u","Exp","Exp","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("ExpExp-u-diff fit to ExpExp-c0-com model: Estimate (in magenta) vs. true model")
##'
##' # Unconstraint fits okay, considering information of continuity was not used.
##'
##'
##' # Example 6: Log-Log with "c0" change point and common variance. - No option of smoothness
##' # Simulate the data
##' x1<-seq(2,4,by=0.05)
##' x2<-seq(4.05,7,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 2 - 2*log(x1)
##' yt2 <- 13.1 - 10*log(x2)
##'
##' # Add noises
##' y1<- yt1+rnorm(length(x1),0,0.5)
##' y2<- yt2+rnorm(length(x2),0,0.5)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##'
##' # It looks like a Log-Log regression with a change point between 3.5 and 4.5
##' # Fit with vrcp with specification of Log-Log segments and "c0" options with 
##' # and without common variance restriction
##'
##' ans <- vrcp(z,3.5,4.5,"c0","Log","Log","Common") # Common variance
##' ans
##'
##' # The fitted Log-Log regression is superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LogLog-c0-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,3.5,4.5,"c0","Log","Log","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LogLog-c0-diff model: Estimate (in magenta) vs. true model")
##'
##' # Fit with vrcp with specification of Log-Log segments and "u" options with 
##' # and without common variance restriction
##'
##' ans <- vrcp(z,3.7,4.5,"u","Log","Log","Common") # Common variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LogLog-u-com fit to LogLog-c0-com model: Estimate (in magenta) vs. true model")
##'
##'
##' ans <- vrcp(z,3.7,4.5,"u","Log","Log","Diff") # Different variance
##' ans

##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LogLog-u-diff fit to LogLog-c0-com model: Estimate (in magenta) vs. true model")
##'
##'
##' # Example 7: Q-Exp with "c1" change point and common variance. 
##' # Simulate the data
##' x1<-seq(0,2,by=0.05)
##' x2<-seq(2,5,by=0.05)
##'
##' # The true regression functions
##' yt1 <-  .2+.2*x1+.5*x1^2
##' yt2 <- .3832+.3*exp(x2)
##'
##' # Add noises
##' y1<- yt1+rnorm(length(x1),0,3)
##' y2<- yt2+rnorm(length(x2),0,3)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##'
##' # It looks like a Q-Exp regression with a change point between 1.5 and 2.2
##' # Fit with vrcp with specification of Q-Exp segments and "c1" options with 
##' # and without common variance restriction
##'
##' ans <- vrcp(z,1.5,2.2,"c1","Q","Exp","Common") # Common variance
##' ans
##'
##' # The fitted Q-Exp regression is superimposed on the data
##' # Let's compare it with the true regression 
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("QExp-c0-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.5,2.2,"c1","Q","Exp","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("QExp-c0-diff model: Estimate (in magenta) vs. true model")
##'
##'
##' # Example 8: Q-Log with "c1" change point and common variance.
##' # Simulate the data
##' x1<-seq(0.05,2.05,by=0.05)
##' x2<-seq(2.1,5.05,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 2+1*x1+5*x1^2
##' yt2 <- 0+35*log(x2)
##'
##' # Add noises
##' y1<-yt1+rnorm(length(x1),0,4)
##' y2<-yt2+rnorm(length(x2),0,4)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##'
##' # It looks like a Q-Log regression with a change point between 1.5 and 2.5
##' # Fit with vrcp with specification of Q-Log segments and "c1" options with 
##' # and without common variance restriction
##' ans <- vrcp(z,1.5,2.5,"c1","Q","Log","Common") # Common variance
##' ans
##'
##' # The fitted Q-Log regression is superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("QLog-c1-com model: Estimate (in magenta) vs. true model")
##'
##' ans <- vrcp(z,1.5,2.5,"c1","Q","Log","Diff") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("QLog-c1-diff model: Estimate (in magenta) vs. true model")
##'
##'
##' # Example 9: Q-Q with "c1" change point and common variance.
##' # Simulate the data
##' x1<-seq(0.05,2,by=0.05)
##' x2<-seq(2.05,5,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 2+10*x1-5*x1^2
##' yt2 <- 30-21*x2+3.5*x2^2
##'
##' # Add noises
##' y1<-yt1+rnorm(length(x1),0,2)
##' y2<-yt2+rnorm(length(x2),0,2)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##'
##' # It looks like a Q-Q regression with a change point between 1.5 and 2.5
##' # Fit with vrcp with specification of Q-Q segments and "c1" options with 
##' # and without common variance restriction
##' ans <- vrcp(z,1.5,2.5,"c1","Q","Q","Common",spline="TRUE") # Common variance
##' ans
##'
##' # The fitted Q-Q regression and spline are superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt), color = "blue", size=1) +
##' ggtitle("QQ-c1-com model: Estimates vs. true model (in blue)")
##'
##' ans <- vrcp(z,1.5,2.5,"c1","Q","Q","Diff",spline="TRUE") # Different variance
##' ans
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt), color = "blue", size=1) +
##' ggtitle("QQ-c1-diff model: Estimates vs. true model (in blue)")
##'
##' # vrcp fits better than splines fitting.
##'
##'
##' # Example 10: L-NLExp with "c1" change point and common variance.
##' # Simulate the data
##' x1<- seq(0,2,by=0.05)
##' x2<- seq(2.05,5,by=0.05)
##'
##' # The true regression functions
##' yt1 <- 10-0.8*x1
##' yt2 <- 8.4 - (-0.78/0.5) + (-0.78/0.5)*(exp(0.5*(x2-2)))
##'
##' # Add noises
##' y1<- yt1+rnorm(length(x1),0,0.5)
##' y2 <- yt2+rnorm(length(x2),0,0.5)
##' z<-data.frame(c(x1,x2),c(y1,y2))
##' names(z)=c("x","y")
##'
##' # z is the simulated data in data frame. Let's visualize it
##' plot(z)
##'
##' # It looks like a L-NLExp regression with a change point between 1.8 and 2.05
##' # Fit with smooth L-NLExp, common or different variances
##'
##' tryCatch(vrcp(z,1.8,2.05,"c1","L","NLExp","Common",start="FALSE"), error=function(e)
##' {return("Try different starting values. If this still fails, try a different nonlinear model 
##' that might be more suitable to data.")})
##' 
##' ans <- vrcp(z,1.8,2.05,"c1","L","NLExp","Common",start="FALSE")
##' ans
##'
##' # The fitted L-NLExp regression is superimposed on the data
##' # Let's compare it with the true regression
##' x<-z$x
##' yt<-c(yt1,yt2)
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) +
##' scale_colour_grey(name = "Model") +
##' ggtitle("LNLExp-c1-com model: Estimate (in magenta) vs. true model")
##'
##' tryCatch(vrcp(z,1.8,2.05,"c1","L","NLExp","Diff",start="FALSE"), error=function(e){
##' return("Try different starting values. If this still fails, try a different nonlinear model 
##' that might be more suitable to data.")})
##' ans <- vrcp(z,1.8,2.05,"c1","L","NLExp","Diff",start="FALSE")
##' ans$plot + ggplot2::geom_line(aes(x = x, y = yt, colour = c("true")), size=1) + 
##' scale_colour_grey(name = "Model") +
##' ggtitle("LNLExp-c1-diff model: Estimate (in magenta) vs. true model")
##' }

vrcp <-
  function(dataset,lo,hi, smooth = c("c0","c1","u"), 
           segment1 = c("L","Q","Log","Exp","NLExp"),
           segment2 = c("L","Q","Log","Exp","NLExp"),
           variance = c("Common", "Diff"), spline=c("FALSE","TRUE"), start){
	lm1 <- NULL
	lm2 <- NULL
	lm3 <- NULL
	start1 <- NULL
	start2 <- NULL
	exp1 <- NULL
	exp2 <- NULL
	log1 <- NULL
    options(warn=-1)
    low<-lo
    high<-hi
    z<-dataset[order(dataset[,1]),]
    x<-z[,1]
    y<-z[,2]
    n<-length(x)
    lo<-length(x[x<=low])
    hi<-length(x[x<=high])
    spline <- match.arg(spline)
    
    smooth <- match.arg(smooth)
    if (smooth == "c0") {
    continuity <- "Cont"
    smooth <- "FALSE"
    }
    else if (smooth == "c1"){
      continuity <- "Cont"
      smooth <- "TRUE"
    }
    else if (smooth == "u"){
      continuity <- "Uncont"
      smooth <- "FALSE"
    }
      
      
    segment1 <- match.arg(segment1)
    segment2 <- match.arg(segment2)
    
    if (segment1 == "L" | segment1 == "Q"){ 
    }
    else if (segment1 == "Log"){
      segment1 <-"E2"
    } 
    else if (segment1 == "Exp"){
      segment1 <-"E1"
    }
    else if (segment1 == "NLExp"){
      segment1 <-"E3"
    }
    else {
      stop("Segment 1 can only be 'L', 'Q', 'Log', 'Exp', or 'NLExp'.")
    }
    if (segment2 == "L" | segment2 == "Q"){ 
    }
    else if (segment2 == "Log"){
      segment2 <-"E2"
    } 
    else if (segment2 == "Exp"){
      segment2 <-"E1"
    }
    else if (segment2 == "NLExp"){
      segment2 <-"E3"
    }
    else {
      stop("Segment 2 can only be 'L', 'Q', 'Log','Exp', or 'NLExp'.")
    }
    choice <- paste(segment1,segment2,sep="")
    met <- c("LL","LQ","QL","QQ","E1E1","E2E2","LE1","LE2","LE3","E1L","E2L","E3L",
             "QE1","QE2","QE3")
    if (is.na(match(choice,met)) == TRUE){
      stop("The method combination is not recognized.")
    }
    
    method <- choice
    variance <- match.arg(variance)

if (method=="QQ") {
      
      if (continuity == "Uncont"){
        
        if (variance == "Common"){
          hats<-llsearch.QQ.C(x,y,n,lo,hi)
          est<-p.est.QQ.C(x,y,n,hats$jhat)
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] +
            est$a2 * z$x[z$x <= est$xj]^2
          z$lm2 = NA
          z$lm2[z$x >= est$xj] = est$b0 + est$b1 * z$x[z$x >= est$xj] +
            est$b2 * z$x[z$x >= est$xj]^2
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,changepoint=est$xj,plot=plot))
        }
        
        if (variance == "Diff"){
          hats<-llsearch.QQ.D(x,y,n,lo,hi)
          est<-p.est.QQ.D(x,y,n,hats$jhat)
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] +
            est$a2 * z$x[z$x <= est$xj]^2
          z$lm2 = NA
          z$lm2[z$x >= est$xj] = est$b0 + est$b1 * z$x[z$x >= est$xj] +
            est$b2 * z$x[z$x >= est$xj]^2
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1,est$b2),
                      sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
        }
      }
      
      if (continuity == "Cont"){
        
        if (smooth == "TRUE"){
          
          if (variance == "Diff"){
            hats<-llsearch.QQ.CDS(x,y,n,lo,hi)
            est<-p.est.QQ.CDS(x,y,n,hats$jhat)
            z$lm1 = NA
            z$lm1[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2*(z$x[z$x <= est$xj])^2
            z$lm2 = NA
            z$lm2[z$x >= est$xj] = est$b0 + est$b1*z$x[z$x >= est$xj] + est$b2*(z$x[z$x >= est$xj])^2
           
            if (spline == "TRUE"){
              z$lm3 = NA
              z$lm3 = predict(lm(y~splines::bs(x, df = 1, knots = est$xj)))
              plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
                ggplot2::geom_line(aes(x = x, y = lm1,color = "vrcp"),size=1) + 
                ggplot2::geom_line(aes(x = x, y = lm2,color = "vrcp"),size=1)+
                ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)+
                ggplot2::geom_line(aes(x = x, y = lm3, color = "spline"),size=1)+
                scale_colour_manual(name ="", values=c(vrcp = "deeppink", spline="orange"))
            }
            else{
            plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
              ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
              ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
              ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
            }
            return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1,est$b2),
                        sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))  
            
            }
          
          if (variance == "Common"){
            hats<-llsearch.QQ.CCS(x,y,n,lo,hi)
            est<-p.est.QQ.CCS(x,y,n,hats$jhat)
            z$lm1 = NA
            z$lm1[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2*(z$x[z$x <= est$xj])^2
            z$lm2 = NA
            z$lm2[z$x >= est$xj] = est$b0 + est$b1*z$x[z$x >= est$xj] + est$b2*(z$x[z$x >= est$xj])^2
            if (spline == "TRUE"){
              z$lm3 = NA
              z$lm3 = predict(lm(y~splines::bs(x, df = 1, knots = est$xj)))
              plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
                ggplot2::geom_line(aes(x = x, y = lm1,color = "vrcp"),size=1) + 
                ggplot2::geom_line(aes(x = x, y = lm2,color = "vrcp"),size=1)+
                ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)+
                ggplot2::geom_line(aes(x = x, y = lm3, color = "spline"),size=1)+
                scale_colour_manual(name ="", values=c(vrcp = "deeppink", spline="orange"))
            }
            else{
            plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
              ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
              ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
              ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
            }
            return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                        sigma2=est$sigma2,changepoint=est$xj,plot=plot))    
            
          }
        }
        
        else{
          if (variance =="Common") {
            hats<-llsearch.QQ.CC(x,y,n,lo,hi)
            est<-p.est.QQ.CC(x,y,n,hats$jhat)
            z$lm1 = NA
            z$lm1[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2*(z$x[z$x <= est$xj])^2
            z$lm2 = NA
            z$lm2[z$x >= est$xj] = est$b0 + est$b1*z$x[z$x >= est$xj] + est$b2*(z$x[z$x >= est$xj])^2
            if (spline == "TRUE"){
              z$lm3 = NA
              z$lm3 = predict(lm(y~splines::bs(x, df = 1, knots = est$xj)))
              plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
                ggplot2::geom_line(aes(x = x, y = lm1,color = "vrcp"),size=1) + 
                ggplot2::geom_line(aes(x = x, y = lm2,color = "vrcp"),size=1)+
                ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)+
                ggplot2::geom_line(aes(x = x, y = lm3, color = "spline"),size=1)+
                scale_colour_manual(name ="", values=c(vrcp = "deeppink", spline="orange"))
            }
            else{
            plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
              ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
              ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
              ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
            }
            return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1,est$b2),
                        sigma2=est$sigma2,changepoint=est$xj,plot=plot))  
          }
          
          if (variance =="Diff") {
            hats<-llsearch.QQ.CD(x,y,n,lo,hi)
            est<-p.est.QQ.CD(x,y,n,hats$jhat)
            z$lm1 = NA
            z$lm1[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2*(z$x[z$x <= est$xj])^2
            z$lm2 = NA
            z$lm2[z$x >= est$xj] = est$b0 + est$b1*z$x[z$x >= est$xj] + est$b2*(z$x[z$x >= est$xj])^2
            if (spline == "TRUE"){
              z$lm3 = NA
              z$lm3 = predict(lm(y~splines::bs(x, df = 1, knots = est$xj)))
              plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
                ggplot2::geom_line(aes(x = x, y = lm1,color = "vrcp"),size=1) + 
                ggplot2::geom_line(aes(x = x, y = lm2,color = "vrcp"),size=1)+
                ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)+
                ggplot2::geom_line(aes(x = x, y = lm3, color = "spline"),size=1)+
                scale_colour_manual(name ="", values=c(vrcp = "deeppink", spline="orange"))
            }
            else{
            plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
              ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
              ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
              ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
            }
            return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                        coe=est$beta,changepoint=est$tau,plot=plot))
          }
        }
      }
    }    
    
if (method=="LL") {
      
      if (continuity == "Uncont"){ 
      
        if (variance == "Common"){
        hats<-llsearch.LL.C(x,y,n,lo,hi)
        est<-p.est.LL.C(x,y,n,hats$jhat)
        z$lm1 = NA
        z$lm1[z$x <= est$xj] = est$a0 + est$a1* z$x[z$x <= est$xj]
        z$lm2 = NA
        z$lm2[z$x >= est$xj] = est$b0 + est$b1 * z$x[z$x >= est$xj]
        if (spline == "TRUE"){
          z$lm3 = NA
          z$lm3 = predict(lm(y~splines::bs(x, df = 1, knots = est$xj)))
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm1,color = "vrcp"),size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm2,color = "vrcp"),size=1)+
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)+
            ggplot2::geom_line(aes(x = x, y = lm3, color = "spline"),size=1)+
            scale_colour_manual(name ="", values=c(vrcp = "deeppink", spline="orange"))
        }
        else{
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        }
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot = plot))
      }
       
      if (variance == "Diff"){
        hats<-llsearch.LL.D(x,y,n,lo,hi)
        est<-p.est.LL.D(x,y,n,hats$jhat)
        z$lm1 = NA
        z$lm1[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$lm2 = NA
        z$lm2[z$x >= est$xj] = est$b0 + est$b1 * z$x[z$x >= est$xj]
        if (spline == "TRUE"){
          z$lm3 = NA
          z$lm3 = predict(lm(y~splines::bs(x, df = 1, knots = est$xj)))
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm1,color = "vrcp"),size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm2,color = "vrcp"),size=1)+
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)+
            ggplot2::geom_line(aes(x = x, y = lm3, color = "spline"),size=1)+
            scale_colour_manual(name ="", values=c(vrcp = "deeppink", spline="orange"))
        }
        else{
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        }
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }
    } 
      
      if (continuity == "Cont"){
        
        if (variance =="Common") {
          hats<-con.search.LL.C(x,y,n,lo,hi)
          est<-con.vals.LL.C(x,y,n,hats$jhat)
          z$lm1 = NA
          z$lm1[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau]
          z$lm2 = NA
          z$lm2[z$x >= est$tau] = est$beta[3] + est$beta[4] * z$x[z$x >= est$tau]
          
          if (spline == "TRUE"){    
            z$lm3 = NA
            z$lm3 = predict(lm(y~splines::bs(x, df = 1, knots = est$tau)))
              
            plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
              ggplot2::geom_line(aes(x = x, y = lm1,color = "vrcp"),size=1) + 
              ggplot2::geom_line(aes(x = x, y = lm2,color = "vrcp"),size=1)+
              ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)+
              ggplot2::geom_line(aes(x = x, y = lm3, color = "spline"),size=1)+
              scale_colour_manual(name ="", values=c(vrcp = "deeppink", spline="orange"))
          }
          else{
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
          }
          return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                      coe=est$beta,changepoint=est$tau,plot=plot))
        }
        
        if (variance =="Diff") {
          hats<-con.search.LL.D(x,y,n,lo,hi)
          est<-con.vals.LL.D(x,y,n,hats$jhat)
          z$lm1 = NA
          z$lm1[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau]
          z$lm2 = NA
          z$lm2[z$x >= est$tau] = est$beta[3] + est$beta[4] * z$x[z$x >= est$tau]
          if (spline == "TRUE"){
            z$lm3 = NA
            z$lm3 = predict(lm(y~splines::bs(x, df = 1, knots = est$tau)))
            plot<- ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
              ggplot2::geom_line(aes(x = x, y = lm1,color = "vrcp"),size=1) + 
              ggplot2::geom_line(aes(x = x, y = lm2,color = "vrcp"),size=1)+
              ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)+
              ggplot2::geom_line(aes(x = x, y = lm3, color = "spline"),size=1)+
              scale_colour_manual(name ="", values=c(vrcp = "deeppink", spline="orange"))
          }
          else{
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
          }
          return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                      coe=est$beta,changepoint=est$tau,plot=plot))
        }
      }
    }  

if (method=="LQ") {
     
  if (continuity == "Uncont"){
      
      if (variance == "Common"){
        hats<-llsearch.LQ.C(x,y,n,lo,hi)
        est<-p.est.LQ.C(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$lm2 = NA
        z$lm2[z$x >= est$xj] = est$b0 + est$b1 * z$x[z$x >= est$xj] +
          est$b2 * z$x[z$x >= est$xj]^2
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.LQ.D(x,y,n,lo,hi)
        est<-p.est.LQ.D(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$lm2 = NA
        z$lm2[z$x >= est$xj] = est$b0 + est$b1 * z$x[z$x >= est$xj] +
          est$b2 * z$x[z$x >= est$xj]^2
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }
     }
     
     if (continuity == "Cont"){
       
       if (smooth == "TRUE"){
         
         if (variance == "Diff"){
           hats<-llsearch.LQ.CDS(x,y,n,lo,hi)
           est<-p.est.LQ.CDS(x,y,n,hats$jhat)
           z$lm1 = NA
           z$lm1[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
           z$lm2 = NA
           z$lm2[z$x >= est$xj] = est$b0 + est$b1*z$x[z$x >= est$xj] + est$b2*(z$x[z$x >= est$xj])^2
           plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
             ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
             ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
             ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
           return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                       sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))  
         }
         
         if (variance == "Common"){
           hats<-llsearch.LQ.CCS(x,y,n,lo,hi)
           est<-p.est.LQ.CCS(x,y,n,hats$jhat)
           z$lm1 = NA
           z$lm1[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
           z$lm2 = NA
           z$lm2[z$x >= est$xj] = est$b0 + est$b1*z$x[z$x >= est$xj] + est$b2*(z$x[z$x >= est$xj])^2
           plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
             ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
             ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
             ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
           return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                       sigma2=est$sigma2,changepoint=est$xj,plot=plot))  
         
       }
       }
       
       else{
       if (variance =="Common") {
         hats<-con.search.LQ.C(x,y,n,lo,hi)
         est<-con.vals.LQ.C(x,y,n,hats$jhat)
         z$lm = NA
         z$lm[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau]
         z$lm2 = NA
         z$lm2[z$x >= est$tau] =  est$beta[3] + est$beta[4] * z$x[z$x >= est$tau] + 
           est$beta[5] * z$x[z$x >= est$tau]^2
         plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
           ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
           ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
           ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
         return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                     coe=est$beta,changepoint=est$tau,plot=plot))
       }
       
       if (variance =="Diff") {
         hats<-con.search.LQ.D(x,y,n,lo,hi)
         est<-con.vals.LQ.D(x,y,n,hats$jhat)
         z$lm = NA
         z$lm[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau]
         z$lm2 = NA
         z$lm2[z$x >= est$tau] =  est$beta[3] + est$beta[4] * z$x[z$x >= est$tau] + 
           est$beta[5] * z$x[z$x >= est$tau]^2
         plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
           ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
           ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
           ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
         return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                     coe=est$beta,changepoint=est$tau,plot=plot))
       }
     }
     }
    }  

if (method=="QL") {
  
  if (continuity == "Uncont"){
      
      if (variance == "Common"){
        hats<-llsearch.QL.C(x,y,n,lo,hi)
        est<-p.est.QL.C(x,y,n,hats$jhat)
        z$lm2 = NA
        z$lm2[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] +
          est$a2 * z$x[z$x <= est$xj]^2
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1)+
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.QL.D(x,y,n,lo,hi)
        est<-p.est.QL.D(x,y,n,hats$jhat)
        z$lm2 = NA
        z$lm2[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] +
          est$a2 * z$x[z$x <= est$xj]^2
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1)+
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }
  }
  
  if (continuity == "Cont"){
    if (smooth == "TRUE"){
      
      if (variance =="Common") {
        hats<-llsearch.QL.CCS(x,y,n,lo,hi)
        est<- p.est.QL.CCS(x,y,n,hats$jhat)
        z$lm1 = NA
        z$lm1[z$x <= est$xj] = est$a0 + est$a1*z$x[z$x<=est$xj] + est$a2 * (z$x[z$x<=est$xj])^2
        z$lm2 = NA
        z$lm2[z$x >= est$xj] = est$b0 + est$b1 * (z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                    sigma2=est$sigma2, changepoint=est$xj,plot=plot))
        
      }
      
    if (variance == "Diff"){
      hats<-llsearch.QL.CDS(x,y,n,lo,hi)
      est<- p.est.QL.CDS(x,y,n,hats$jhat)
      z$lm1 = NA
      z$lm1[z$x <= est$xj] = est$a0 + est$a1*z$x[z$x<=est$xj] + est$a2 * z$x[z$x<=est$xj]^2
      z$lm2 = NA
      z$lm2[z$x >= est$xj] = est$b0 + est$b1 * z$x[z$x >= est$xj]
      plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
        ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
        ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1)+
        ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
      return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                  sigma2=est$sigma2,tau2=est$tau2, changepoint=est$xj,plot=plot))
      
    }   
  }
    else{
    if (variance =="Common") {
      hats<-con.search.QL.C(x,y,n,lo,hi)
      est<-con.vals.QL.C(x,y,n,hats$jhat)
      z$lm2 = NA
      z$lm2[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau] +
        est$beta[3] * (z$x[z$x <= est$tau])^2
      z$lm = NA
      z$lm[z$x >= est$tau] =  est$beta[4] + est$beta[5] * z$x[z$x >= est$tau]
      plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
        ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
        ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1)+
        ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
      return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                  coe=est$beta,changepoint=est$tau,plot=plot))
    }
    
    if (variance =="Diff") {
      hats<-con.search.QL.D(x,y,n,lo,hi)
      est<-con.vals.QL.D(x,y,n,hats$jhat)
      z$lm2 = NA
      z$lm2[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau] +
        est$beta[3] * z$x[z$x <= est$tau]^2
      z$lm = NA
      z$lm[z$x >= est$tau] =  est$beta[4] + est$beta[5] * z$x[z$x >= est$tau]
      plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
        ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
        ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1)+
        ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
      return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                  coe=est$beta,changepoint=est$tau,plot=plot))
    }
  }  
}
}

if (method=="LE1") {
  
  if (continuity == "Uncont"){
      
      if (variance == "Common"){
        hats<-llsearch.LE1.C(x,y,n,lo,hi)
        est<-p.est.LE1.C(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$exp = NA
        z$exp[z$x >= est$xj] = est$b0 + est$b1 * exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.LE1.D(x,y,n,lo,hi)
        est<-p.est.LE1.D(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$exp = NA
        z$exp[z$x >= est$xj] = est$b0 + est$b1 * exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }
  }
  
  if (continuity == "Cont"){
    
    if (variance =="Common") {
      
      if (smooth == "TRUE"){
        hats<-llsearch.LE1.CCS(x,y,n,lo,hi)
        est<-p.est.LE1.CCS(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$exp = NA
        # ypred <- est$a0 + est$a1*est$xj
        z$exp[z$x >= est$xj] = est$b0 + est$b1*exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1)+
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))  
      }
      
      else{
        hats<-con.search.LE1.C(x,y,n,lo,hi)
        est<-con.vals.LE1.C(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$tau] = est$beta[1] + est$beta[2]* z$x[z$x <= est$tau]
        z$exp = NA
        z$exp[z$x >= est$tau] = est$beta[3] + est$beta[4] * exp(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
    }
    
    if (variance =="Diff") {
      
      if (smooth == "TRUE"){
        hats<-llsearch.LE1.CDS(x,y,n,lo,hi)
        est<-p.est.LE1.CDS(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$exp = NA
        ypred <- est$a0 + est$a1*est$xj
        z$exp[z$x >= est$xj] = est$b0 + est$b1*exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1)+
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))  
      }
      
      else{
        hats<-con.search.LE1.D(x,y,n,lo,hi)
        est<-con.vals.LE1.D(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$tau] = est$beta[1] + est$beta[2]* z$x[z$x <= est$tau]
        z$exp = NA
        z$exp[z$x >= est$tau] = est$beta[3] + est$beta[4] * exp(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
    }  
  }
  
}

if (method=="LE2") {
  
  if (continuity == "Uncont"){
      
      if (variance == "Common"){
        hats<-llsearch.LE2.C(x,y,n,lo,hi)
        est<-p.est.LE2.C(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$log = NA
        z$log[z$x >= est$xj] = est$b0 + est$b1 * log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.LE2.D(x,y,n,lo,hi)
        est<-p.est.LE2.D(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$log = NA
        z$log[z$x >= est$xj] = est$b0 + est$b1 * log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }    
}
  
  if (continuity == "Cont"){
    
    if (smooth == "TRUE"){
      
      if (variance == "Common"){
        hats<-llsearch.LE2.CCS(x,y,n,lo,hi)
        est<-p.est.LE2.CCS(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$log = NA
        z$log[z$x >= est$xj] = est$b0 + est$b1 * log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if(variance == "Diff"){
        hats<-llsearch.LE2.CDS(x,y,n,lo,hi)
        est<-p.est.LE2.CDS(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$log = NA
        z$log[z$x >= est$xj] = est$b0 + est$b1 * log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))
      }
      
    }
    
    else{
      
      if (variance =="Common") {
        hats<-con.search.LE2.C(x,y,n,lo,hi)
        est<-con.vals.LE2.C(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$tau] = est$beta[1] + est$beta[2]* z$x[z$x <= est$tau]
        z$log = NA
        z$log[z$x >= est$tau] = est$beta[3] + est$beta[4] * log(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=1/est$eta[1],coe=est$beta,
                    changepoint=est$tau,plot=plot))
      }
      
      if (variance =="Diff") {
        hats<-con.search.LE2.D(x,y,n,lo,hi)
        est<-con.vals.LE2.D(x,y,n,hats$jhat)
        z$lm = NA
        z$lm[z$x <= est$tau] = est$beta[1] + est$beta[2]* z$x[z$x <= est$tau]
        z$log = NA
        z$log[z$x >= est$tau] = est$beta[3] + est$beta[4] * log(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
    }
  }
}

if (method=="LE3") {
  
  if (variance == "Common"){
    
    if (continuity == "Uncont"){
      
      if(start == "FALSE"){
        est<-llsearch.LE3.C.WITHOUT.I(x,y,n,lo,hi)
      }
      
      else{ # with initial values
        hats<-llsearch.LE3.C(x,y,n,lo,hi,start[1],start[2],start[3])
        est<-p.est.LE3.C(x,y,n,hats$jhat,start[1],start[2],start[3])
      }    
    
      z$lm = NA
      z$lm[z$x < est$xj] = est$a0 + est$a1 * z$x[z$x < est$xj]
      z$exp = NA
      ypred <- est$a0 + est$a1*est$xj
      z$exp[z$x > est$xj] = est$b0 + est$b1*exp(est$b2 * (z$x[z$x > est$xj] - est$xj))

      plot <- ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) +
        ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
        ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
        ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
      
      return(list(maxloglik=est$maxloglik,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                  sigma2=est$sigma2,changepoint=est$xj,plot=plot))
    }
    
    if (continuity == "Cont"){
      if (smooth == "TRUE"){
        if (start == "FALSE"){
          initials <- llsearch.LE3.C.WITHOUT.I(x,y,n,lo,hi)
          hats<-llsearch.LE3.CCS.WITHOUT.I(x,y,n,lo,hi,initials$a0,initials$a1,initials$b0,initials$b1,initials$b2)
          est<-p.est.LE3.CCS.WITHOUT.I(x,y,n,hats$jhat,initials$a0,initials$a1,initials$b0,initials$b1,initials$b2)
        }       
        else{
          hats<-llsearch.LE3.CCS(x,y,n,lo,hi,start[1],start[2],start[3])
          est<-p.est.LE3.CCS(x,y,n,hats$jhat,start[1],start[2],start[3])
        }
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$exp = NA
        ypred <- est$a0 + est$a1*est$xj
        z$exp[z$x >= est$xj] = est$b0 + est$b1*exp(est$b2 * (z$x[z$x >= est$xj] - est$xj))
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))     
      }
      
      else{ # c0: continuous
        if (start == "FALSE"){
          est <- llsearch.LE3.CC.WITHOUT.I(x,y,n,lo,hi)
          
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
          z$exp = NA
          ypred <- est$a0 + est$a1*est$xj
          z$exp[z$x >= est$xj] = est$b0 + est$b1*exp(est$b2 * (z$x[z$x >= est$xj] - est$xj))
          
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          
          return(list(maxloglik=est$maxloglik,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,changepoint=est$xj,plot=plot))
        }
        else{
          hats<-llsearch.LE3.CC(x,y,n,lo,hi,start1,start2)
          est<-p.est.LE3.CC(x,y,n,hats$jhat,start1,start2)
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
          z$exp = NA
          ypred <- est$a0 + est$a1*est$xj
          z$exp[z$x >= est$xj] = ypred - est$b0 * (1 - exp(est$b1 * (z$x[z$x >= est$xj] - est$xj)) )
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,changepoint=est$xj,plot=plot))
        }
      }
    }
  }
  
  if(variance == "Diff"){
    
    if (continuity == "Uncont"){
      
      if(start == "FALSE"){
        est<-llsearch.LE3.D.WITHOUT.I(x,y,n,lo,hi)
      }
      else{ # with initial values
        hats<-llsearch.LE3.CD(x,y,n,lo,hi,start[1],start[2])
        est<-p.est.LE3.CD(x,y,n,hats$jhat,start[1],start[2])
      }    
      z$lm = NA
      z$lm[z$x < est$xj] = est$a0 + est$a1 * z$x[z$x < est$xj]
      z$exp = NA
      z$exp[z$x > est$xj] = est$b0 + est$b1*(exp(est$b2*(z$x[z$x > est$xj] - est$xj))) 
      plot <- ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) +
        ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
        ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
        ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
      
      return(list(maxloglik=est$maxloglik,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                  sigma2=est$sigma2,eau2=est$tau2,changepoint=est$xj,plot=plot))
    }
    
    if (continuity == "Cont"){
      if (smooth == "TRUE"){
        if (start == "FALSE"){
          initials <- llsearch.LE3.C.WITHOUT.I(x,y,n,lo,hi)
          hats<-llsearch.LE3.CDS.WITHOUT.I(x,y,n,lo,hi,initials$a0,initials$a1,initials$b0,initials$b1,initials$b2)
          est<-p.est.LE3.CDS.WITHOUT.I(x,y,n,hats$jhat,initials$a0,initials$a1,initials$b0,initials$b1,initials$b2)
        }       
        else{
          hats<-llsearch.LE3.CDS(x,y,n,lo,hi,start[1],start[2],start[3])
          est<-p.est.LE3.CDS(x,y,n,hats$jhat,start[1],start[2],start[3])
        }
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
        z$exp = NA
        ypred <- est$a0 + est$a1*est$xj
        z$exp[z$x >= est$xj] = est$b0 + est$b1*exp(est$b2 * (z$x[z$x >= est$xj] - est$xj))
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                    sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))     
      }
      
      else{ # c0: continuous
        if (start == "FALSE"){
          est<-llsearch.LE3.D.WITHOUT.I(x,y,n,lo,hi)
          
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
          z$exp = NA
          ypred <- est$a0 + est$a1*est$xj
          z$exp[z$x >= est$xj] = est$b0 + est$b1*exp(est$b2 * (z$x[z$x >= est$xj] - est$xj))
          
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          
          return(list(maxloglik=est$maxloglik,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))
        }
        else{
          hats<-llsearch.LE3.CD(x,y,n,lo,hi,start1,start2)
          est<-p.est.LE3.CD(x,y,n,hats$jhat,start1,start2)
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj]
          z$exp = NA
          ypred <- est$a0 + est$a1*est$xj
          z$exp[z$x >= est$xj] = ypred - est$b0 * (1 - exp(est$b1 * (z$x[z$x >= est$xj] - est$xj)) )
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))
        }
      }
    }   
  }
}

if (method=="E1L") {
  
    if(continuity == "Uncont"){ 
      
      if (variance == "Common"){
        hats<-llsearch.E1L.C(x,y,n,lo,hi)
        est<-p.est.E1L.C(x,y,n,hats$jhat)
        z$exp = NA
        z$exp[z$x <= est$xj] = est$a0 + est$a1 * exp(z$x[z$x <= est$xj])
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.E1L.D(x,y,n,lo,hi)
        est<-p.est.E1L.D(x,y,n,hats$jhat)
        z$exp = NA
        z$exp[z$x <= est$xj] = est$a0 + est$a1 * exp(z$x[z$x <= est$xj])
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }
    }
    
    if (continuity == "Cont"){
      
      if (variance =="Common") {
        
        if (smooth == "TRUE"){
          hats<-llsearch.E1L.CCS(x,y,n,lo,hi)
          est<-p.est.E1L.CCS(x,y,n,hats$jhat)
          z$exp = NA
          z$exp[z$x <= est$xj] = est$a0 + est$a1*exp(z$x[z$x <= est$xj])
          z$lm = NA
          z$lm[z$x >= est$xj] = est$b0 + est$b1*z$x[z$x >= est$xj]
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                      sigma2=est$sigma2,changepoint=est$xj,plot=plot))  
        }
        
        else{
        hats<-con.search.E1L.C(x,y,n,lo,hi)
        est<-con.vals.E1L.C(x,y,n,hats$jhat)
        z$exp = NA
        z$exp[z$x <= est$tau] = est$beta[1] + est$beta[2] * exp(z$x[z$x <= est$tau])
        z$lm = NA
        z$lm[z$x >= est$tau] = est$beta[3] + est$beta[4] * z$x[z$x >= est$tau]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
    }
      
      if (variance =="Diff") {
        
        if (smooth == "TRUE"){
          hats<-llsearch.E1L.CDS(x,y,n,lo,hi)
          est<-p.est.E1L.CDS(x,y,n,hats$jhat)
          z$exp = NA
          z$exp[z$x <= est$xj] = est$a0 + est$a1*exp(z$x[z$x <= est$xj])
          z$lm = NA
          z$lm[z$x >= est$xj] = est$b0 + est$b1*z$x[z$x >= est$xj]
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                      sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))  
        }
        
        else{
        hats<-con.search.E1L.D(x,y,n,lo,hi)
        est<-con.vals.E1L.D(x,y,n,hats$jhat)
        z$exp = NA
        z$exp[z$x <= est$tau] = est$beta[1] + est$beta[2] * exp(z$x[z$x <= est$tau])
        z$lm = NA
        z$lm[z$x >= est$tau] =  est$beta[3] + est$beta[4] * z$x[z$x >= est$tau]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
    }     
  }
}

if (method=="E2L") {
  
  if(continuity == "Uncont"){

      if (variance == "Common"){
        hats<-llsearch.E2L.C(x,y,n,lo,hi)
        est<-p.est.E2L.C(x,y,n,hats$jhat)
        z$log = NA
        z$log[z$x <= est$xj] = est$a0 + est$a1 * log(z$x[z$x <= est$xj])
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.E2L.D(x,y,n,lo,hi)
        est<-p.est.E2L.D(x,y,n,hats$jhat)
        z$log = NA
        z$log[z$x <= est$xj] = est$a0 + est$a1 * log(z$x[z$x <= est$xj])
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }  
    } 
  
  if (continuity == "Cont"){
    
      if (smooth == "TRUE"){
        
        if (variance == "Common"){
        hats<-llsearch.E2L.CCS(x,y,n,lo,hi)
        est<-p.est.E2L.CCS(x,y,n,hats$jhat)
        z$log = NA
        z$log[z$x <= est$xj] = est$a0 + est$a1 * log(z$x[z$x <= est$xj])
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.E2L.CDS(x,y,n,lo,hi)
        est<-p.est.E2L.CDS(x,y,n,hats$jhat)
        z$log = NA
        z$log[z$x <= est$xj] = est$a0 + est$a1 * log(z$x[z$x <= est$xj])
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))
      }
      
      }
      
      else{
        if (variance =="Common") {
          hats<-con.search.E2L.C(x,y,n,lo,hi)
          est<-con.vals.E2L.C(x,y,n,hats$jhat)
          z$log = NA
          z$log[z$x <= est$tau] = est$beta[1] + est$beta[2] * log(z$x[z$x <= est$tau])
          z$lm = NA
          z$lm[z$x >= est$tau] =  est$beta[3] + est$beta[4] * z$x[z$x >= est$tau]
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                      coe=est$beta,changepoint=est$tau,plot=plot))
        }
        
        if (variance =="Diff") {
          hats<-con.search.E2L.D(x,y,n,lo,hi)
          est<-con.vals.E2L.D(x,y,n,hats$jhat)
          z$log = NA
          z$log[z$x <= est$tau] = est$beta[1] + est$beta[2] * log(z$x[z$x <= est$tau])
          z$lm = NA
          z$lm[z$x >= est$tau] =  est$beta[3] + est$beta[4] * z$x[z$x >= est$tau]
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                      coe=est$beta,changepoint=est$tau,plot=plot))
        }
      } 
    }
  }
    
if (method=="E3L"){
  
  if (continuity == "Uncont"){
      
  if (variance == "Common"){
    hats<-llsearch.E3L.C(x,y,n,lo,hi,start[1],start[2],start[3])
    est<-p.est.E3L.C(x,y,n,hats$jhat,start[1],start[2],start[3])
    z$exp = NA
    z$exp[z$x <= est$xj] = est$a0 + est$a1 * exp(est$a2 * (z$x[z$x <= est$xj] - est$xj))
    z$lm = NA
    z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
    plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
      ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
      ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
      ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
    return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                sigma2=est$sigma2,changepoint=est$xj,plot=plot))
  }
  
  if (variance == "Diff"){
    hats<-llsearch.E3L.D(x,y,n,lo,hi)
    est<-p.est.E3L.D(x,y,n,hats$jhat)
    z$exp = NA
    z$exp[z$x <= est$xj] = est$a0 * (1 - exp(est$a1 * (z$x[z$x <= est$xj] - est$a2)))
    z$lm = NA
    z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
    plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
      ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
      ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
      ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
    return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
  }   
 }
 
 if (continuity == "Cont"){  
   if (smooth == "TRUE"){     
     if (variance =="Common") {
       
     if (start != "FALSE"){
       hats<-llsearch.E3L.CCS(x,y,n,lo,hi,start[1],start[2],start[3])
       est<- p.est.E3L.CCS(x,y,n,hats$jhat,start[1],start[2],start[3])
     }
       z$exp = NA
       z$exp[z$x <= est$xj] = est$a0 + est$a1 * exp(est$a2*(z$x[z$x <= est$xj]-est$xj))
       z$lm = NA
       z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
       plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
         ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
         ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
         ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
       return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                   sigma2=est$sigma2,changepoint=est$xj,plot=plot))
     }
     
     if (variance =="Diff") {
     if( start != "FALSE"){
       hats<-llsearch.E3L.CDS(x,y,n,lo,hi,start[1],start[2],start[3])
       est<- p.est.E3L.CDS(x,y,n,hats$jhat,start[1],start[2],start[3])
     }
       z$exp = NA
       z$exp[z$x <= est$xj] = est$a0 + est$a1 * exp(est$a2*(z$x[z$x <= est$xj]-est$xj))
       z$lm = NA
       z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
       plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
         ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
         ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
         ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
       return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                   sigma2=est$sigma2,tau2 = est$tau2, changepoint=est$xj,plot=plot))
     }
   }
     
   else {
       if (variance =="Common") {
         hats<-llsearch.E3L.CC(x,y,n,lo,hi,start[1],start[2]) # 0.4 2
         est<- p.est.E3L.CC(x,y,n,hats$jhat,start[1],start[2])
         z$exp = NA
         ypred <- est$b0 + est$b1*est$xj
         z$exp[z$x <= est$xj] = ypred - est$a1 + est$a1 * exp(est$a2*(z$x[z$x <= est$xj]-est$xj))
         z$lm = NA
         z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
         plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
           ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
           ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
           ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
         return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                     sigma2=est$sigma2,changepoint=est$xj,plot=plot))    
       }
       
       if (variance =="Diff") {
         hats<-llsearch.E3L.CD(x,y,n,lo,hi,start[1],start[2])
         est<- p.est.E3L.CD(x,y,n,hats$jhat,start[1],start[2])
         z$exp = NA
         ypred <- est$b0 + est$b1*est$xj
         z$exp[z$x <= est$xj] = ypred - est$a0 + est$a0 * exp(est$a1*(z$x[z$x <= est$xj]-est$xj))
         z$lm = NA
         z$lm[z$x >= est$xj] =  est$b0 + est$b1 * z$x[z$x >= est$xj]
         plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
           ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
           ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
           ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
         return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2, est$b0,est$b1),
                     sigma2=est$sigma2,tau2 = est$tau2, changepoint=est$xj,plot=plot))
       
     }
   }
   
 }
}

if (method=="QE1") {
  if (continuity == "Uncont"){
      if (variance == "Common"){
        hats<-llsearch.QE1.C(x,y,n,lo,hi)
        est<-p.est.QE1.C(x,y,n,hats$jhat)
        z$lm2 = NA
        z$lm2[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2 
        z$exp = NA
        z$exp[z$x >= est$xj] = est$b0 + est$b1 * exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.QE1.D(x,y,n,lo,hi)
        est<-p.est.QE1.D(x,y,n,hats$jhat)
        z$lm2 = NA
        z$lm2[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2 
        z$exp = NA
        z$exp[z$x >= est$xj] = est$b0 + est$b1 * exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }      
    }  
  if (continuity == "Cont"){
    
    if (smooth == "TRUE"){
      if (variance == "Common"){
        hats<-llsearch.QE1.CCS(x,y,n,lo,hi)
        est<- p.est.QE1.CCS(x,y,n,hats$jhat)
        z$exp = NA
        z$lm1 = NA
        z$lm1[z$x <= est$xj] = est$a0 + est$a1* z$x[z$x <= est$xj] + est$a2* (z$x[z$x <= est$xj])^2
        z$lm2 = NA
        z$lm2[z$x >= est$xj] = est$b0 + est$b1 * exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2, changepoint=est$xj,plot=plot))
      }
      if (variance == "Diff"){
        hats<-llsearch.QE1.CDS(x,y,n,lo,hi)
        est<- p.est.QE1.CDS(x,y,n,hats$jhat)
        z$exp = NA
        ypred <- est$b0 + est$b1*est$xj
        z$exp[z$x <= est$xj] = est$a0 + est$a1*z$x[z$x<=est$xj] + est$a2 * z$x[z$x<=est$xj]^2
        z$lm = NA
        z$lm[z$x >= est$xj] =  est$b0 + est$b1 * exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,tau2 = est$tau2, changepoint=est$xj,plot=plot))
      }
    }
    
    else{
    if (variance =="Common") {
      hats<-con.search.QE1.C(x,y,n,lo,hi)
      est<-con.vals.QE1.C(x,y,n,hats$jhat)
      z$lm2 = NA
      z$lm2[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau] +
        est$beta[3] * z$x[z$x <= est$tau]^2 
      z$exp = NA
      z$exp[z$x >= est$tau] = est$beta[4] + est$beta[5] * exp(z$x[z$x >= est$tau])
      plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
        ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
        ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
        ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
      return(list(maxloglik=hats$value,sigma2=1/est$eta[1],coe=est$beta,
                  changepoint=est$tau,plot=plot))
    }
    
    if (variance =="Diff") {
      hats<-con.search.QE1.D(x,y,n,lo,hi)
      est<-con.vals.QE1.D(x,y,n,hats$jhat)
      z$lm2 = NA
      z$lm2[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau] + 
        est$beta[3] * z$x[z$x <= est$tau]^2 
      z$exp = NA
      z$exp[z$x >= est$tau] = est$beta[4] + est$beta[5] * exp(z$x[z$x >= est$tau])
      plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
        ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
        ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
        ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
      return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                  coe=est$beta,changepoint=est$tau,plot=plot))
    }
  }
}
}

if (method=="QE2") {
    
  if (continuity == "Uncont"){
      if (variance == "Common"){
        hats<-llsearch.QE2.C(x,y,n,lo,hi)
        est<-p.est.QE2.C(x,y,n,hats$jhat)
        z$lm2 = NA
        z$lm2[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2 
        z$log = NA
        z$log[z$x >= est$xj] = est$b0 + est$b1* log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
    
      
      if (variance == "Diff"){
        hats<-llsearch.QE2.D(x,y,n,lo,hi)
        est<-p.est.QE2.D(x,y,n,hats$jhat)
        z$lm2 = NA
        z$lm2[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2 
        z$log = NA
        z$log[z$x >= est$xj] = est$b0 + est$b1* log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$a2,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }
  }
  
  if (continuity == "Cont"){
    
    if (smooth == "TRUE"){
      
      if (variance == "Common"){
        hats<-llsearch.QE2.CCS(x,y,n,lo,hi)
        est<- p.est.QE2.CCS(x,y,n,hats$jhat)
        z$lm = NA
        ypred <- est$b0 + est$b1*est$xj
        z$lm[z$x <= est$xj] = est$a0 + est$a1*z$x[z$x<= est$xj] + est$a2 * z$x[z$x<=est$xj]^2
        z$log = NA
        z$log[z$x >= est$xj] =  est$b0 + est$b1 * log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2, changepoint=est$xj,plot=plot))
      }
      
      if(variance == "Diff"){
        hats<-llsearch.QE2.CDS(x,y,n,lo,hi)
        est<- p.est.QE2.CDS(x,y,n,hats$jhat)
        z$lm = NA
        ypred <- est$b0 + est$b1*est$xj
        z$lm[z$x <= est$xj] = est$a0 + est$a1*z$x[z$x<=est$xj] + est$a2 * z$x[z$x<=est$xj]^2
        z$log = NA
        z$log[z$x >= est$xj] =  est$b0 + est$b1 * log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,tau2 = est$tau2, changepoint=est$xj, plot=plot))
        
      }
    }
    
    else{
      if (variance =="Common") {
        hats<-con.search.QE2.C(x,y,n,lo,hi)
        est<-con.vals.QE2.C(x,y,n,hats$jhat)
        z$lm2 = NA
        z$lm2[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau] + 
          est$beta[3] * z$x[z$x <= est$tau]^2 
        z$log = NA
        z$log[z$x >= est$tau] = est$beta[4] + est$beta[5] * log(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
      
      if (variance =="Diff") {
        hats<-con.search.QE2.D(x,y,n,lo,hi)
        est<-con.vals.QE2.D(x,y,n,hats$jhat)
        z$lm2 = NA
        z$lm2[z$x <= est$tau] = est$beta[1] + est$beta[2] * z$x[z$x <= est$tau] + 
          est$beta[3] * z$x[z$x <= est$tau]^2 
        z$log = NA
        z$log[z$x >= est$tau] = est$beta[4] + est$beta[5] * log(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm2), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
    }
  }    
}  

if (method=="QE3") {
  
  if (variance == "Common"){
    
    if (continuity == "Uncont"){
      
      if(start == "FALSE"){
        est<-llsearch.QE3.C.WITHOUT.I(x,y,n,lo,hi)
      }
      
      else{ # with initial values
        hats<-llsearch.QE3.C(x,y,n,lo,hi,start[1],start[2],start[3])
        est<-p.est.QE3.C(x,y,n,hats$jhat,start[1],start[2],start[3])
      }    
      z$lm = NA
      z$lm[z$x < est$xj] = est$a0 + est$a1 * z$x[z$x < est$xj] + est$a2 * z$x[z$x < est$xj]^2
      z$exp = NA
      z$exp[z$x > est$xj] = est$b0 + est$b1*(exp(est$b2*(z$x[z$x > est$xj] - est$xj)))
      
      plot <- ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) +
        ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
        ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
        ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
      
      return(list(maxloglik=est$maxloglik,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                  sigma2=est$sigma2,changepoint=est$xj,plot=plot))
    }
    
    if (continuity == "Cont"){
      if (smooth == "TRUE"){
        if (start == "FALSE"){
          initials <- llsearch.QE3.C.WITHOUT.I(x,y,n,lo,hi)
          hats<-llsearch.QE3.CCS(x,y,n,lo,hi,initials$a0,initials$a1,initials$a2,initials$b2)
          est<-p.est.QE3.CCS(x,y,n,hats$jhat,initials$a0,initials$a1,initials$a2,initials$b2)
        }       
        else{
          hats<-llsearch.QE3.CCS(x,y,n,lo,hi,start[1],start[2],start[3],start[4])
          est<-p.est.QE3.CCS(x,y,n,hats$jhat,start[1],start[2],start[3],start[4])
        }
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2
        z$exp = NA
        ypred <- est$a0 + est$a1*est$xj
        z$exp[z$x >= est$xj] = est$b0 + est$b1*exp(est$b2 * (z$x[z$x >= est$xj] - est$xj))
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))     
      }
      
      else{ # c0: continuous
        if (start == "FALSE"){
          est <- llsearch.QE3.CC.WITHOUT.I(x,y,n,lo,hi)
          
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2
          z$exp = NA
          z$exp[z$x >= est$xj] = est$b0 + est$b1*(exp(est$b2*(z$x[z$x >= est$xj] - est$xj)))
          
          plot <- ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) +
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          
          return(list(maxloglik=est$maxloglik,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,changepoint=est$xj,plot=plot))
        }
        else{
          hats<-llsearch.QE3.CC(x,y,n,lo,hi,start1,start2) 
          est<-p.est.QE3.CC(x,y,n,hats$jhat,start1,start2)
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2
          z$exp = NA
          ypred <- est$a0 + est$a1*est$xj + est$a2 * est$xj^2
          z$exp[z$x >= est$xj] = ypred - est$b0 * (1 - exp(est$b1 * (z$x[z$x >= est$xj] - est$xj)) )
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,changepoint=est$xj,plot=plot))
        }
      }
    }
  }
  
  if(variance == "Diff"){
    
    if (continuity == "Uncont"){
      
      if(start == "FALSE"){
        est<-llsearch.QE3.D.WITHOUT.I(x,y,n,lo,hi)
      }
      else{ # with initial values
        hats<-llsearch.QE3.D(x,y,n,lo,hi,start[1],start[2],start[3]) ## llsearch.QE3.D(x,y,101,36,42,4,6,0.3)
        est<-p.est.QE3.D(x,y,n,hats$jhat,start[1],start[2],start[3])
      }    
      z$lm = NA
      z$lm[z$x < est$xj] = est$a0 + est$a1 * z$x[z$x < est$xj] + est$a2 * z$x[z$x < est$xj]^2
      z$exp = NA
      ypred <- est$a0 + est$a1 * est$xj + est$a2 * est$xj^2
      z$exp[z$x > est$xj] = ypred - est$b1 + est$b1*(exp(est$b2*(z$x[z$x > est$xj] - est$xj))) 
      plot <- ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) +
        ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
        ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
        ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
      
      return(list(maxloglik=est$maxloglik,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                  sigma2=est$sigma2,eau2=est$tau2,changepoint=est$xj,plot=plot))
    }
    
    if (continuity == "Cont"){
      if (smooth == "TRUE"){
        if (start == "FALSE"){
          initials <- llsearch.QE3.C.WITHOUT.I(x,y,n,lo,hi)
          hats<-llsearch.QE3.CDS(x,y,n,lo,hi,initials$b0,initials$b1)
          est<-p.est.QE3.CDS(x,y,n,hats$jhat,initials$b0,initials$b1,initials$b2)
        }       
        else{
          hats<-llsearch.QE3.CDS(x,y,n,lo,hi,start[1],start[2])
          est<-p.est.QE3.CDS(x,y,n,hats$jhat,start[1],start[2])
        }
        z$lm = NA
        z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2
        z$exp = NA
        z$exp[z$x >= est$xj] = est$b0 + est$b1*exp(est$b2 * (z$x[z$x >= est$xj] - est$xj))
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                    sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))     
      }
      
      else{ # c0: continuous
        if (start == "FALSE"){
          est<-llsearch.QE3.D.WITHOUT.I(x,y,n,lo,hi) 
          
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2
          z$exp = NA
          z$exp[z$x >= est$xj] = est$b0 + est$b1*(exp(est$b2*(z$x[z$x >= est$xj] - est$xj)))
          
          plot <- ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) +
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) +
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          
          return(list(maxloglik=est$maxloglik,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))
        }
        else{ 
          hats<-llsearch.QE3.CD(x,y,n,lo,hi,start1,start2)
          est<-p.est.QE3.CD(x,y,n,hats$jhat,start1,start2)
          z$lm = NA
          z$lm[z$x <= est$xj] = est$a0 + est$a1 * z$x[z$x <= est$xj] + est$a2 * z$x[z$x <= est$xj]^2
          z$exp = NA
          ypred <- est$a0 + est$a1*est$xj + est$a2*est$xj^2
          z$exp[z$x >= est$xj] = ypred - est$b0 * (1 - exp(est$b1 * (z$x[z$x >= est$xj] - est$xj)) )
          plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
            ggplot2::geom_line(aes(x = x, y = lm), color = "deeppink",size=1) + 
            ggplot2::geom_line(aes(x = x, y = exp), color = "deeppink",size=1) +
            ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
          return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$b2),
                      sigma2=est$sigma2,tau2=est$tau2,changepoint=est$xj,plot=plot))
        }
      }
    }   
  }
}

if (method=="E1E1") {
  
  if (continuity == "Uncont"){ 
    
      if (variance == "Common"){
        hats<-llsearch.E1E1.C(x,y,n,lo,hi)
        est<-p.est.E1E1.C(x,y,n,hats$jhat)
        z$exp1 = NA
        z$exp1[z$x <= est$xj] = est$a0 + est$a1 * exp(z$x[z$x <= est$xj])
        z$exp2 = NA
        z$exp2[z$x >= est$xj] = est$b0 + est$b1 * exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.E1E1.D(x,y,n,lo,hi)
        est<-p.est.E1E1.D(x,y,n,hats$jhat)
        z$exp1 = NA
        z$exp1[z$x <= est$xj] = est$a0 + est$a1 * exp(z$x[z$x <= est$xj])
        z$exp2 = NA
        z$exp2[z$x >= est$xj] = est$b0 + est$b1 * exp(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }    
    }
  
  if (continuity == "Cont"){

      if (variance =="Common") {
        hats<-con.search.E1E1.C(x,y,n,lo,hi)
        est<-con.vals.E1E1.C(x,y,n,hats$jhat)
        z$exp1 = NA
        z$exp1[z$x <= est$tau] = est$beta[1] + est$beta[2] * exp(z$x[z$x <= est$tau])
        z$exp2 = NA
        z$exp2[z$x >= est$tau] = est$beta[3] + est$beta[4] * exp(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
      
      if (variance =="Diff") {
        hats<-con.search.E1E1.D(x,y,n,lo,hi)
        est<-con.vals.E1E1.D(x,y,n,hats$jhat)
        z$exp1 = NA
        z$exp1[z$x <= est$tau] = est$beta[1] + est$beta[2] * exp(z$x[z$x <= est$tau])
        z$exp2 = NA
        z$exp2[z$x >= est$tau] = est$beta[3] + est$beta[4] * exp(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = exp1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = exp2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
  }
  
 }

if (method=="E2E2") {
  if (continuity == "Uncont"){
      
      if (variance == "Common"){
        hats<-llsearch.E2E2.C(x,y,n,lo,hi)
        est<-p.est.E2E2.C(x,y,n,hats$jhat)
        z$log1 = NA
        z$log1[z$x <= est$xj] = est$a0 + est$a1 * log(z$x[z$x <= est$xj])
        z$log2 = NA
        z$log2[z$x >= est$xj] = est$b0 + est$b1 * log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = log1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=est$sigma2,changepoint=est$xj,plot=plot))
      }
      
      if (variance == "Diff"){
        hats<-llsearch.E2E2.D(x,y,n,lo,hi)
        est<-p.est.E2E2.D(x,y,n,hats$jhat)
        z$log1 = NA
        z$log1[z$x <= est$xj] = est$a0 + est$a1 * log(z$x[z$x <= est$xj])
        z$log2 = NA
        z$log2[z$x >= est$xj] = est$b0 + est$b1 * log(z$x[z$x >= est$xj])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = log1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$xj,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1),
                    sigma2=c(est$sigma2,est$tau2),changepoint=est$xj,plot=plot))
      }
    }
  
  if (continuity == "Cont"){
      
      if (variance =="Common") {
        hats<-con.search.E2E2.C(x,y,n,lo,hi)
        est<-con.vals.E2E2.C(x,y,n,hats$jhat)
        z$log1 = NA
        z$log1[z$x <= est$tau] = est$beta[1] + est$beta[2] * log(z$x[z$x <= est$tau])
        z$log2 = NA
        z$log2[z$x >= est$tau] = est$beta[3] + est$beta[4] * log(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = log1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=1/est$eta[1],
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
      
      if (variance =="Diff") {
        hats<-con.search.E2E2.D(x,y,n,lo,hi)
        est<-con.vals.E2E2.D(x,y,n,hats$jhat)
        z$log1 = NA
        z$log1[z$x <= est$tau] = est$beta[1] + est$beta[2] * log(z$x[z$x <= est$tau])
        z$log2 = NA
        z$log2[z$x >= est$tau] = est$beta[3] + est$beta[4] * log(z$x[z$x >= est$tau])
        plot<-ggplot2::ggplot(z) + ggplot2::geom_point(aes(x = x, y = y)) + 
          ggplot2::geom_line(aes(x = x, y = log1), color = "deeppink",size=1) + 
          ggplot2::geom_line(aes(x = x, y = log2), color = "deeppink",size=1) +
          ggplot2::geom_vline(xintercept =  est$tau,linetype="dotted",size=1)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),
                    coe=est$beta,changepoint=est$tau,plot=plot))
      }
      
    } 
  }

}
