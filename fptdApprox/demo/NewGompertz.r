#############################################################################################
#############################################################################################
###                                                                                       ###
###                                       GOMPERTZ-TYPE PROCESS                           ###
###                                                                                       ###
### Examples for constant and another boundaries associated to some time random variables ###
### related to the growth of a Gompertz-type diffusion process with infinitesimal moments ###
### A1(x,t) = m*x*exp(-beta*t)  and  A2(x,t) = (sigma^2)*(x^2)                            ###
#############################################################################################
#############################################################################################
                                                              
# Creating the diffproc object "Gompertz"
NewGompertz <- diffproc(c("m*x*exp(-beta*t)","sigma^2*x^2","dnorm((log(x)-log(y) + (m/beta)*(exp(-beta*t) - exp(-beta*s))
                     + (t-s)*sigma^2/2)/(sigma*sqrt(t-s)),0,1)/(x*sigma*sqrt(t-s))","plnorm(x,log(y)
                     - (m/beta)*(exp(-beta*t) - exp(-beta*s)) - (t-s)*sigma^2/2,sigma*sqrt(t-s))"))
NewGompertz

################################################################################################
# CONSTANT BOUNDARY  S                                                                         #
# Application to real data of the paper:                                                       #
#                                                                                              #
# R. Guti\u00E9rrez, P. Rom\u00E1n, D. Romero, J.J. Serrano, and F. Torres (2008)              #
# "Some time random variables related to a Gompertz-type diffusion process"                    #
# Cybernetics and Systems: An International Journal, 39, 1-13                                  #
################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y1NG.FPTL <- FPTL(dp = NewGompertz, t0 = 1, T = 30, x0 = 144, S = 2500, list(m = 0.755152, beta=0.183128, sigma = 0.0708605))
# Displaying graphically the FPTL function
plot(y1NG.FPTL)
# Extracting the interesting information provided by the FPTL function
y1NG.SFPTL <- summary(y1NG.FPTL)
# Approximating the f.p.t density
y1NG.g <- Approx.cfpt.density(y1NG.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y1NG.g)
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y1NG.g, report.sfptl=TRUE)

#############################################
# Approximating directly the f.p.t density. #
#############################################

# Approximating the f.p.t density
y1NG.g1 <- Approx.fpt.density(dp = NewGompertz, t0 = 1, T = 30, id = 144, S = 2500, list(m = 0.755152, beta=0.183128, sigma = 0.0708605))
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y1NG.g1)
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y1NG.g1, report.sfptl=TRUE)

################################################################################################
# CONSTANT BOUNDARY  S = x0*exp((m/beta)*exp(-beta*t0) - 1)                                    #
# Inflection time in Application to real data of the paper:                                    #
#                                                                                              #
# R. Guti\u00E9rrez, P. Rom\u00E1n, D. Romero, J.J. Serrano, and F. Torres (2008)              #
# "Some time random variables related to a Gompertz-type diffusion process"                    #
# Cybernetics and Systems: An International Journal, 39, 1-13                                  #
################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y2NG.FPTL <- FPTL(dp = NewGompertz, t0 = 1, T = 30, x0 = 144, S = "x0*exp((m/beta)*exp(-beta*t0) - 1)", list(m = 0.755152, beta=0.183128, sigma = 0.0708605))
# Displaying graphically the FPTL function
win.graph()
plot(y2NG.FPTL)
# Extracting the interesting information provided by the FPTL function
y2NG.SFPTL <- summary(y2NG.FPTL)
# Approximating the f.p.t density
y2NG.g <- Approx.cfpt.density(y2NG.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y2NG.g)
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y2NG.g, report.sfptl=TRUE)

#############################################
# Approximating directly the f.p.t density. #
#############################################

# Approximating the f.p.t density
y2NG.g1 <- Approx.fpt.density(dp = NewGompertz, t0 = 1, T = 30, id = 144, S = "x0*exp((m/beta)*exp(-beta*t0) - 1)", list(m = 0.755152, beta=0.183128, sigma = 0.0708605))
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y2NG.g1)
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y2NG.g1, report.sfptl=TRUE)

################################################################################################
# CONSTANT BOUNDARY  S = p*x0*exp((m/beta)*exp(-beta*t0))                                      #
# Time at which the process achieves a percentage 100*p% of the total growth in                # 
# Application to real data of the paper:                                                       #
#                                                                                              #
# R. Guti\u00E9rrez, P. Rom\u00E1n, D. Romero, J.J. Serrano, and F. Torres (2008)              #
# "Some time random variables related to a Gompertz-type diffusion process"                    #
# Cybernetics and Systems: An International Journal, 39, 1-13                                  #
################################################################################################

#################################################
# Approximating step-by-step the f.p.t density. #
#################################################

# Evaluating the FPTL function
y3NG.FPTL <- FPTL(dp = NewGompertz, t0 = 1, T = 30, x0 = 144, S = "p*x0*exp((m/beta)*exp(-beta*t0))", list(m = 0.755152, beta=0.183128, sigma = 0.0708605, p=0.5))
# Displaying graphically the FPTL function
win.graph()
plot(y3NG.FPTL)
# Extracting the interesting information provided by the FPTL function
y3NG.SFPTL <- summary(y3NG.FPTL)
# Approximating the f.p.t density
y3NG.g <- Approx.cfpt.density(y3NG.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y3NG.g)
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y3NG.g, report.sfptl=TRUE)

#############################################
# Approximating directly the f.p.t density. #
#############################################

# Approximating the f.p.t density
y3NG.g1 <- Approx.fpt.density(dp = NewGompertz, t0 = 1, T = 30, id = 144, S = "p*x0*exp((m/beta)*exp(-beta*t0))", list(m = 0.755152, beta=0.183128, sigma = 0.0708605, p=0.5))
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y3NG.g1)
# Reporting information of the approximation process of the f.p.t. density and the interesting 
# information provided by the FPTL function
report(y3NG.g1, report.sfptl=TRUE)
