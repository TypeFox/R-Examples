########################################################
########################################################
###                                                  ###
###            ORNSTEIN UHLENBECK PROCESS            ###
###                                                  ###
########################################################
########################################################

# Creating the diffproc object "OU"
OU <- diffproc(c("alpha*x + beta","sigma^2","dnorm((x-(y*exp(alpha*(t-s)) - beta*(1 - exp(alpha*(t-s)))/alpha))/
              (sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha))),0,1)/(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))",
              "pnorm(x, y*exp(alpha*(t-s)) - beta*(1 - exp(alpha*(t-s)))/alpha,sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))"))
OU
            
######################################################################################################
#                                      CONSTANT BOUNDARY S = 2                                       #
#  for an Ornstein Uhlenbeck process with infinitesimal moments A1(x,t) = - 2x + 2 and A2(x,t) = 1,  #
#  t0 = 0 and x0 = 1.    (x0 < S)                                                                    #
######################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y1OU.FPTL <- FPTL(OU, 0, 70, 1, 2, list(alpha=-2,beta=2,sigma=1))
# Displaying graphically the FPTL function
win.graph()
plot(y1OU.FPTL)
# Extracting the interesting information provided by the FPTL function
y1OU.SFPTL <- summary(y1OU.FPTL)
# Reporting the interesting information provided by the FPTL function
report(y1OU.SFPTL)
# Approximating the f.p.t density
y1OU.g <- Approx.cfpt.density(y1OU.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y1OU.g)
# Reporting information of the approximation process of the f.p.t. density
report(y1OU.g)

# The value of the cumulative integral of the approximation is less than 1 - tol.
# The maximum value of the FPTL function is very small and this function increases up to tmax 
# and then remains constant. Thus, it may be appropriate to increase the value of the final  
# stopping instant or approximate the density again to T (with to.T = TRUE). Since T is much 
# larger than the final stopping instant, we try to increase that value using k argument to 
# summary the fptl class object
y1OUa.SFPTL <- summary(y1OU.FPTL, k=4)
y1OUa.SFPTL
y1OUa.g <- Approx.cfpt.density(y1OUa.SFPTL)
y1OUa.g
win.graph()
plot(y1OUa.g)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y1OU.g1 <- Approx.fpt.density(OU, 0, 70, 1, 2, list(alpha=-2,beta=2,sigma=1))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y1OU.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y1OU.g1)

# Approximating the f.p.t density with k=4 to increase the value of the final stopping instant
y1OUa.g1 <- Approx.fpt.density(OU, 0, 70, 1, 2, list(alpha=-2,beta=2,sigma=1), k=4)
y1OUa.g1

#####################################################################################################
#                              OSCILLATING BOUNDARY S(t) = 2+cos(2*pi*t/5)                          #
#  for an Ornstein Uhlenbeck process with infinitesimal moments A1(x,t) = - x + 1 and A2(x,t) = 2,  #
#  t0 = 0 and x0 = 0.     (x0 < S(t0))                                                              #
#####################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y2OU.FPTL <- FPTL(OU, 0, 70, 0, "2+cos(2*pi*t/5)", list(alpha=-1,beta=1,sigma=sqrt(2)))
# Displaying graphically the FPTL function
win.graph()
plot(y2OU.FPTL)
# Extracting the interesting information provided by the FPTL function
y2OU.SFPTL <- summary(y2OU.FPTL)
# Approximating the f.p.t density
y2OU.g <- Approx.cfpt.density(y2OU.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y2OU.g)
# Reporting information of the approximation process of the f.p.t. density
report(y2OU.g, sfptl=TRUE)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y2OU.g1 <- Approx.fpt.density(OU, 0, 70, 0, "2+cos(2*pi*t/5)", list(alpha=-1,beta=1,sigma=sqrt(2)))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y2OU.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y2OU.g1)

####################################################################################################################
#  Boundary                                                                                                        #
#                             S(t) = - (beta/alpha) + A*exp(alpha*t) + B*exp(-alpha*t)                             #
#  for an Ornstein Uhlenbeck process with infinitesimal moments A1(x,t) = alpha*x + beta  and  A2(x,t) = sigma^2.  #  
#  In this case the f.p.t. density is known                                                                        #
####################################################################################################################

# Theoretical expression of the f.p.t. density through the boundary 
gOU <- parse(text="2*alpha*abs(A*exp(alpha*t0)+B*exp(-alpha*t0)-x0-(beta/alpha))*dnorm(-(beta/alpha)+A*exp(alpha*t)+B*exp(-alpha*t), 
            x0*exp(alpha*(t-t0))-beta*(1-exp(alpha*(t-t0)))/alpha, sigma*sqrt((exp(2*alpha*(t-t0))-1)/(2*alpha)))/(exp(alpha*(t-t0))-exp(-alpha*(t-t0)))")

###############################################################
#  A1(x,t) = -2*x + 2, A2(x,t) = 1, x0 = -1, S = 1  (x0 < S)  #
###############################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y3OU.FPTL <- FPTL(OU, 0, 70, -1, 1, list(alpha=-2,beta=2,sigma=1))
# Displaying graphically the FPTL function
win.graph()
plot(y3OU.FPTL)
# Extracting the interesting information provided by the FPTL function
y3OU.SFPTL <- summary(y3OU.FPTL)
y3OU.SFPTL
# Approximating the f.p.t density
y3OU.g <- Approx.cfpt.density(y3OU.SFPTL)
# Reporting information of the approximation process of the f.p.t. density
report(y3OU.g)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y3OU.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y3OU.g$x, eval(gOU, list(alpha=-2, beta=2, sigma=1, A=0, B=0, t0=0, x0=-1, t=y3OU.g$x)), type="l", col=2)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y3OU.g1 <- Approx.fpt.density(OU, 0, 70, -1, 1, list(alpha=-2,beta=2,sigma=1))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y3OU.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y3OU.g1)

##############################################################
#  A1(x,t) = -2*x + 2, A2(x,t) = 1, x0 = 3, S = 1  (x0 > S)  #
##############################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y4OU.FPTL <- FPTL(OU, 0, 70, 3, 1, list(alpha=-2,beta=2,sigma=1))
# Displaying graphically the FPTL function
win.graph()
plot(y4OU.FPTL)
# Extracting the interesting information provided by the FPTL function
y4OU.SFPTL <- summary(y4OU.FPTL)
# Approximating the f.p.t density
y4OU.g <- Approx.cfpt.density(y4OU.SFPTL)
# Reporting information of the approximation process of the f.p.t. density
report(y4OU.g)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y4OU.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y4OU.g$x, eval(gOU, list(alpha=-2, beta=2, sigma=1, A=0, B=0, t0=0, x0=3, t=y4OU.g$x)), type="l", col=2)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y4OU.g1 <- Approx.fpt.density(OU, 0, 70, 3, 1, list(alpha=-2,beta=2,sigma=1))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y4OU.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y4OU.g1)

###################################################################################################
#  A1(x,t) = -0.2*x + 1, A2(x,t) = 1, x0 = 1, t0 = 0, S(t) = 5 + 10*exp(-0.2*t)     (x0 < S(t0))  #
###################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function         
y5OU.FPTL <- FPTL(OU, 0, 70, 1, "- beta/alpha + A*exp(alpha*t)", list(alpha=-0.2,beta=1,sigma=1,A=10))
# Displaying graphically the FPTL function
win.graph()
plot(y5OU.FPTL)
# Extracting the interesting information provided by the FPTL function
y5OU.SFPTL <- summary(y5OU.FPTL)
# Approximating the f.p.t density
y5OU.g <- Approx.cfpt.density(y5OU.SFPTL)
# Reporting information of the approximation process of the f.p.t. density
report(y5OU.g)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y5OU.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y5OU.g$x, eval(gOU, list(alpha=-0.2, beta=1, sigma=1, A=10, B=0, t0=0, x0=1, t=y5OU.g$x)), type="l", col=2)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y5OU.g1 <- Approx.fpt.density(OU, 0, 70, 1, "- beta/alpha + A*exp(alpha*t)", list(alpha=-0.2,beta=1,sigma=1,A=10))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y5OU.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y5OU.g1)

##############################################################################################################
#  A1(x,t) = -0.1*x + 2, A2(x,t) = 1, x0 = 1, t0 = 0, S(t) = 20 + exp(-0.1*t) - exp(0.1*t)     (x0 < S(t0))  #
##############################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y6OU.FPTL <- FPTL(OU, 0, 70,1, "- beta/alpha + A*exp(alpha*t) + B*exp(-alpha*t)", list(alpha=-0.1,beta=2,sigma=1,A=1,B=-1))
# Displaying graphically the FPTL function
win.graph()
plot(y6OU.FPTL)
# Extracting the interesting information provided by the FPTL function
y6OU.SFPTL <- summary(y6OU.FPTL)
# Approximating the f.p.t density
y6OU.g <- Approx.cfpt.density(y6OU.SFPTL)
# Reporting information of the approximation process of the f.p.t. density
report(y6OU.g)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y6OU.g)
# Superimposing the plot of the theoretical f.p.t. density
points(y6OU.g$x, eval(gOU, list(alpha=-0.1, beta=2, sigma=1, A=1, B=-1, t0=0, x0=1, t=y6OU.g$x)), type="l", col=2)

#############################################
# Approximating directly the f.p.t density. #
#############################################

y6OU.g1 <- Approx.fpt.density(OU, 0, 70,1, "- beta/alpha + A*exp(alpha*t) + B*exp(-alpha*t)", list(alpha=-0.1,beta=2,sigma=1,A=1,B=-1))
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y6OU.g1, report.sfptl = TRUE)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y6OU.g1)
                  