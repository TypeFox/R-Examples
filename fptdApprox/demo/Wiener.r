##########################################################################################
##########################################################################################
###                                                                                    ###
###                                     WIENER PROCESS                                 ###
###                                                                                    ###
### Examples for constant and linear boundaries for which the f.p.t. density is known  ###
###                                                                                    ###
##########################################################################################
##########################################################################################

# Creating the diffproc object "Wiener"
Wiener <- diffproc(c("m","sigma^2","dnorm((x-(y+m*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s))",
                     "pnorm(x,y+m*(t-s),sigma*sqrt(t-s))"))
Wiener

# Theoretical expression of the f.p.t. density through a linear boundary S(t) = a * t + b
gW <- expression(abs(a * t0 + b - x0) * dnorm(a * t + b, x0 + m * (t - t0), sigma * sqrt(t - t0))/(t - t0))

###################################################################################
# CONSTANT BOUNDARY                                                               #
#                                    S = 4                                        #
# for a Wienner process with infinitesimal moments A1(x,t) = 1  and  A2(x,t) = 1, #
# and x0 = 0   (S > x0)                                                           #
###################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y1W.FPTL <- FPTL(dp = Wiener, t0 = 0, T = 20, x0 = 0, S = 4, list(m = 1, sigma = 1))
# Displaying graphically the FPTL function
plot(y1W.FPTL)
# Extracting and showing the interesting information provided by the FPTL function
y1W.SFPTL <- summary(y1W.FPTL)
y1W.SFPTL
# Reporting the interesting information provided by the FPTL function
report(y1W.SFPTL)
# Approximating the f.p.t density
y1W.g <- Approx.cfpt.density(y1W.SFPTL)
# Reporting information of the approximation process of the f.p.t. density
report(y1W.g)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y1W.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y1W.g$x, eval(gW, list(t0 = 0, x0 = 0, m = 1, sigma = 1, a = 0, b = 4, t = y1W.g$x)), type="l", col=2)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y1W.g1 <- Approx.fpt.density(dp = Wiener, t0 = 0, T = 20, id = 0, S = 4, list(m = 1, sigma = 1))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y1W.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y1W.g1)

###################################################################################
# LINEAR BOUNDARY                                                                 #
#                             S(t) = 10 - 0.5*t                                   #
# for a Wienner process with infinitesimal moments A1(x,t) = 1  and  A2(x,t) = 1, #
# t0 = 0 and x0 = 0   (S(t0) > x0)                                                #
###################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y2W.FPTL <- FPTL(dp = Wiener, t0 = 0, T = 20, x0 = 0, S = "10-0.5*t", list(m = 1, sigma = 1))
# Displaying graphically the FPTL function
win.graph()
plot(y2W.FPTL)
# Extracting and showing the interesting information provided by the FPTL function
y2W.SFPTL <- summary(y2W.FPTL)
y2W.SFPTL
# Reporting the interesting information provided by the FPTL function
report(y2W.SFPTL)
# Approximating the f.p.t density
y2W.g <- Approx.cfpt.density(y2W.SFPTL)
# Reporting information of the approximation process of the f.p.t. density
report(y2W.g)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y2W.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y2W.g$x,eval(gW, list(t0 = 0, x0 = 0, m = 1, sigma = 1, a = -0.5, b = 10, t = y2W.g$x)), type="l", col=2)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y2W.g1 <- Approx.fpt.density(dp = Wiener, t0 = 0, T = 20, id = 0, S = "10-0.5*t", list(m = 1, sigma = 1))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y2W.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y2W.g1)

###################################################################################
# LINEAR BOUNDARY                                                                 #
#                              S(t) = -1 + t/2                                    #
# for a Wienner process with infinitesimal moments A1(x,t) = 0  and  A2(x,t) = 1, #
# t0 = 1 and x0 = 1   (S(t0) < x0)                                                #
###################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y3W.FPTL <- FPTL(dp = Wiener, t0 = 1, T = 40, x0 = 1, S = "-1+t/2", list(m = 0, sigma = 1))
# Displaying graphically the FPTL function
win.graph()
plot(y3W.FPTL)
# Extracting the interesting information provided by the FPTL function
y3W.SFPTL <- summary(y3W.FPTL)
# Reporting the interesting information provided by the FPTL function
report(y3W.SFPTL)
# Approximating the f.p.t density
y3W.g <- Approx.cfpt.density(y3W.SFPTL)
# Reporting information of the approximation process of the f.p.t. density
report(y3W.g, report.sfptl=TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y3W.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y3W.g$x, eval(gW, list(t0 = 1, x0 = 1, m = 0, sigma = 1, a = 0.5, b = -1, t = y3W.g$x)), type="l",col=2)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y3W.g1 <- Approx.fpt.density(dp = Wiener, t0 = 1, T = 40, id = 1, S = "-1+t/2", list(m = 0, sigma = 1))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y3W.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y3W.g1)  
