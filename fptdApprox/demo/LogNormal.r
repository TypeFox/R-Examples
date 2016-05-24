#################################################################################################
#################################################################################################
###                                                                                           ###
###                                   LOGNORMAL PROCESS                                       ###
###                                                                                           ###
### Examples for:                                                                             ###
###               - Constant and exponential boundaries for which the f.p.t. density is known ###
###               - General boundaries                                                        ###
###                                                                                           ###
#################################################################################################
#################################################################################################

# Creating the diffproc object "Lognormal"
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/
                     (sigma*sqrt(t-s)*x)","plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))"))
Lognormal

# Theoretical expression of the f.p.t. density through an exponential boundary A*exp(B*t)
gL <- expression(abs(log(A)-log(x0)+B*t0)*exp(-(log(A)-log(x0)+B*t-(m-sigma^2/2)*(t-t0))^2/(2*sigma^2*(t-t0)))/
                 (sigma*sqrt(2*pi*(t-t0)^3)))

#####################################################################################################
# CONSTANT BOUNDARY                                                                                 #
#                                     S=2500                                                        #
# for a Lognormal process with infinitesimal moments A1(x,t) = 4*x  and A2(x,t) = (0.001^2)*(x^2),  #
# and x0 = 1.     (S > x0)                                                                          #
#                                                                                                   #
# Case of a very concentrated f.p.t. variable                                                       #
# Example 1 of the paper:                                                                           #
#                                                                                                   #
# P. Rom\u00E1n, J.J. Serrano, F. Torres (2008)                                                     #
# "First-passage-time location function: Application to determine first-pasage-times densities in   #
# diffusion processes"                                                                              #
# Computational Statistics and Data Analysis, 52, 4132-4146                                         #                                                                                              
#####################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y1L.FPTL <- FPTL(dp = Lognormal, t0 = 0, T = 2, x0 = 1, S = 2500, list(m = 4, sigma = 0.001))
# Displaying graphically the FPTL function
win.graph()
plot(y1L.FPTL)
# Displaying graphically the FPTL function  (with more adequate options)
win.graph()
plot(y1L.FPTL, from.t0 = FALSE, to.T = FALSE)
# Extracting the interesting information provided by the FPTL function
y1L.SFPTL <- summary(y1L.FPTL)

# Approximating the f.p.t density with default options (variableStep = TRUE, from.t0 = FALSE, 
# to.T = FALSE, skip = TRUE)
y1L.g <- Approx.cfpt.density(y1L.SFPTL)
# Reporting information of the approximation process of the f.p.t. density
report(y1L.g)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(y1L.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y1L.g$x, eval(gL, list(t0=0,x0=1,m=4,sigma=0.001,A=2500,B=0,t=y1L.g$x)), type="l", col=2)
# Displaying graphically the approximation of the f.p.t density together with the information provided by 
# the FPTL function
win.graph()
plot(y1L.g, growth.points=T, instants=T)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y1L.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 2, id = 1, S = 2500, list(m = 4, sigma = 0.001))
# Reporting information of the approximation process of the f.p.t. density
report(y1L.g1, report.sfptl = TRUE)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(y1L.g1)

# Approximating the f.p.t density from t0 and variable integration step by using p and alpha arguments to reduce
# significantly the computational cost.

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################
 
y2L.g <- Approx.cfpt.density(y1L.SFPTL, from.t0=TRUE, n=200, p=0.1, alpha=0.9)
y2L.g
plot(y2L.g, from.t0=FALSE)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y2L.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 2, id = 1, S = 2500, list(m = 4, sigma = 0.001), from.t0=TRUE, n=200, p=0.1, alpha=0.9)
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y2L.g1, report.sfptl = TRUE)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(y2L.g1, from.t0=FALSE)

# Approximating the f.p.t density from t0 and fixed integration step
# Total iterations = 193021 (it is recommended not run)
# Approx.cfpt.density(yyCSDA, variableStep = FALSE, from.t0=TRUE)

#######################################################################################################
# EXPONENTIAL BOUNDARY                                                                                #
#                                     S(t) = 0.5*exp(0.4*t)                                           #
# for a Lognormal process with infinitesimal moments A1(x,t) = 0.2*x  and  A2(x,t) = (0.25^2)*(x^2),  #
# t0 = 0 and x0 = 1.     (S(t0) < x0)                                                                 #
#######################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y3L.FPTL <- FPTL(dp = Lognormal, t0 = 0, T = 25, x0 = 1, S = "0.5*exp(0.4*t)", list(m = 0.2, sigma = 0.25))
# Displaying graphically the FPTL function
win.graph()
plot(y3L.FPTL)
# Extracting the interesting information provided by the FPTL function
y3L.SFPTL <- summary(y3L.FPTL)
# Approximating the f.p.t density 
y3L.g <- Approx.cfpt.density(y3L.SFPTL)
y3L.g
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y3L.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y3L.g$x, eval(gL, list(t0=0,x0=1,m=0.2,sigma=0.25,A=0.5,B=0.4,t=y3L.g$x)), type="l", col=2)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y3L.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 25, id = 1, S = "0.5*exp(0.4*t)", list(m = 0.2, sigma = 0.25))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y3L.g1, report.sfptl = TRUE)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(y3L.g1)

#######################################################################################################
# GENERAL BOUNDARY                                                                                    #
#                           S(t) = 7 + 3.2*t + 1.4*t*sin(1.75*t)                                      #
# for a Lognormal process with infinitesimal moments A1(x,t) = 0.48*x  and  A2(x,t) = (0.07^2)*(x^2), #
# t0 = 0 and x0 = 1.     (S(t0) > x0)                                                                 #
#######################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y4L.FPTL <- FPTL(dp = Lognormal, t0 = 0, T = 10, x0 = 1, S = "7+3.2*t+1.4*t*sin(1.75*t)", list(m = 0.48, sigma = 0.07))
# Displaying graphically the FPTL function
win.graph()
plot(y4L.FPTL)
# Extracting the interesting information provided by the FPTL function
y4L.SFPTL <- summary(y4L.FPTL)
# Reporting the interesting information provided by the FPTL function
report(y4L.SFPTL)
# Reporting the interesting information provided by the FPTL function (in Latex format)
report(y4L.SFPTL, tex=TRUE)
# Approximating the f.p.t density            
y4L.g <- Approx.cfpt.density(y4L.SFPTL)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(y4L.g)
# Reporting information of the approximation process of the f.p.t. density
report(y4L.g)
# Reporting information of the approximation process of the f.p.t. density (in Latex format)
report(y4L.g, tex=TRUE)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y4L.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 10, id = 1, S = "7+3.2*t+1.4*t*sin(1.75*t)", list(m = 0.48, sigma = 0.07))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y4L.g1, report.sfptl = TRUE)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(y4L.g1)

#######################################################################################################
# GENERAL BOUNDARY                                                                                    #
#                           S(t) = 4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))                           #
# for a Lognormal process with infinitesimal moments A1(x,t) = 0.48*x  and  A2(x,t) = (0.07^2)*(x^2), #
# t0 = 0 and x0 = 1.     (S(t0) > x0)                                                                 #
#######################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y5L.FPTL <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = "4.5+4*t^2+7*t*sqrt(t)*sin(6*sqrt(t))", list(m = 0.48, sigma = 0.07))
# Displaying graphically the FPTL function
win.graph()
plot(y5L.FPTL)
# Extracting the interesting information provided by the FPTL function
y5L.SFPTL <- summary(y5L.FPTL)
# Reporting the interesting information provided by the FPTL function
report(y5L.SFPTL)
# Approximating the f.p.t density 
y5L.g <- Approx.cfpt.density(y5L.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y5L.g)
# Reporting information of the approximation process of the f.p.t. density
report(y5L.g)

# Extracting the interesting information provided by the FPTL function with a value for the zeroSlope 
# lower than the default considered
y5La.SFPTL <- summary(y5L.FPTL, zeroSlope=0.001)
# Reporting the interesting information provided by the FPTL function
report(y5La.SFPTL)
# Displaying graphically the FPTL function
# Note that in this case the FPTL exhibits 3 growth subintervals
win.graph()
plot(y5L.FPTL, y5La.SFPTL)
# Approximating the f.p.t density 
y5La.g <- Approx.cfpt.density(y5La.SFPTL)
y5La.g
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y5La.g) 

# Approximating the f.p.t density with variable integration step and not avoiding subintervals
y5Lb.g <- Approx.cfpt.density(y5L.SFPTL, skip=FALSE)
y5Lb.g  

# Approximating the f.p.t density with variable integration step, not avoiding subintervals, from t0 and to T
y5Lc.g <- Approx.cfpt.density(y5L.SFPTL, from.t0=TRUE, to.T=TRUE, skip=FALSE)
y5Lc.g
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y5Lc.g) 

# Approximating the f.p.t density with fixed integration step, from t0 and not avoiding subintervals
y5Ld.g <- Approx.cfpt.density(y5L.SFPTL, variableStep=FALSE, from.t0=TRUE, skip=FALSE, n=200)
y5Ld.g
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y5Ld.g) 

#############################################
# Approximating directly the f.p.t density. #
#############################################

y5L.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5+4*t^2+7*t*sqrt(t)*sin(6*sqrt(t))", list(m = 0.48, sigma = 0.07))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y5L.g1, report.sfptl = TRUE)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(y5L.g1)

# Approximating the f.p.t density with a value for the zeroSlope lower than the default considered
y5La.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5+4*t^2+7*t*sqrt(t)*sin(6*sqrt(t))", list(m = 0.48, sigma = 0.07), zeroSlope=0.001)
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y5La.g1, report.sfptl = TRUE)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(y5La.g1)

# Approximating the f.p.t density with variable integration step and not avoiding subintervals
y5Lb.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5+4*t^2+7*t*sqrt(t)*sin(6*sqrt(t))", list(m = 0.48, sigma = 0.07), skip=FALSE)
y5Lb.g1  

# Approximating the f.p.t density with variable integration step, not avoiding subintervals, from t0 and to T
y5Lc.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5+4*t^2+7*t*sqrt(t)*sin(6*sqrt(t))", list(m = 0.48, sigma = 0.07), from.t0=TRUE, to.T=TRUE, skip=FALSE)
y5Lc.g1

# Approximating the f.p.t density with fixed integration step, from t0 and not avoiding subintervals
y5Ld.g1 <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5+4*t^2+7*t*sqrt(t)*sin(6*sqrt(t))", list(m = 0.48, sigma = 0.07), variableStep=FALSE, from.t0=TRUE, skip=FALSE, n=200)
y5Ld.g1
