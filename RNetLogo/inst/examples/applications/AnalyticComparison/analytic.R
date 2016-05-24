###############################################
# Preparation of Ryacas package 
###############################################

# install Ryacas, if currently not installed
#install.packages("Ryacas")

# load Ryacas
library(Ryacas)

# install Yacas, if currently not installed 
# (just for Windows - see Ryacas documentation for other systems)
#yacasInstall()

# tells us, where Yacas is installed
# (just for Windows - see RYacas documentation for other systems)
#yacasFile()



###############################################
# A. ABM Simulation/Empirical part 
###############################################

#-------------------
# A.I. Initialize simulation
#-------------------
library(RNetLogo)
#----
# TODO: adapt this path to your NetLogo installation folder
#----
nl.path <- "C:/Program Files/NetLogo 5.3/app"
# initialize NetLogo
NLStart(nl.path, gui=FALSE)
# load Gas Lab model from model library
model.path <- "models/Sample Models/Chemistry & Physics/GasLab/GasLab Free Gas.nlogo"
NLLoadModel(paste(nl.path,model.path,sep="/"))
# initialize simulation
NLCommand("set number-of-particles 500", "no-display", "setup")
# run simulation for 40 times of 50 steps (= 2000 simulation steps)
# save speed of particles after every 50 simulation step interval
particles.speed <- NLDoReport(40, "repeat 50 [go]", "[speed] of particles")
# flatten the list of lists (one list for each of the 40 runs) to one big vector
particles.speed.vector <- unlist(particles.speed)
# plot histogram
hist(particles.speed.vector, breaks=max(particles.speed.vector), freq=FALSE)


###############################################
# B. Theoretical part 
###############################################

#-------------------
# B.I. Equation B
#-------------------
# get mean energy from NetLogo simulation
energy.mean <- NLReport("mean [energy] of particles")
# definition of function B
B <- function(v, m=1, k=1) v * exp((-m*v^2)/(2*k*energy.mean))
# register function B in Yacas
yacas(B)
# test the function B with a value
#testB1 <- expression(B(1))
#t <- yacas(N(testB1))
#Eval(t)

#-------------------
# B.II. Integration
#-------------------
# integrate function B from 0 to endless
B.integr <- expression(integrate(B,0,Infinity))
# register intergration expression in Yacas
yacas(B.integr)
# calculate a numerical approximation using Yacas function N()
normalizer.yacas <- yacas(N(B.integr))
# get result from Yacas in R
normalizer <- Eval(normalizer.yacas)
# the numeric result is in column value
#normalizer$value

#-------------------
# B.III. Calculation of theoretical value
#-------------------
# get max. speed from NetLogo simulation
maxspeed <- max(particles.speed.vector) 
# create a sequence vector from 0 to maxspeed + stepsize, by stepsize
stepsize <- 0.25
v.vec <- seq(0, maxspeed, stepsize)
# calculate the theoretical value at the points of the sequence vector
theoretical <- B(v.vec) / normalizer$value
# plot the theoretical values
plot(v.vec ,theoretical)


###############################################
# C. Plot simulation and theoretical part together 
###############################################

# merge both plots
hist(particles.speed.vector, breaks=max(particles.speed.vector)*5, freq=FALSE,
     xlim=c(0,as.integer(maxspeed)+5), 
     main="Histogram of empirical speed together \nwith theoretical Maxwell-Boltzmann distribution",
     xlab="speed of particles")
lines(v.vec, theoretical, lwd=2, col="blue")

NLQuit()