### R code from vignette source 'RSGHB_HowTo.rnw'

###################################################
### code chunk number 1: RSGHB_HowTo.rnw:44-72
###################################################

library(RSGHB)

# Load example data
data(choicedata)

# This data set is organized as one row per choice observation. This isn't necessary 
# but it does make writing the likelihood function straightforward.
head(choicedata)

# Any variables or transformations of variables from the choicedata data.frame that are
# needed for evaluating the likelihood function can be extracted as a series of vectors
# for convenience (e.g., TT1, TT2, TOLL2). Alternatively, the likelihood function could
# refer to the choicedata data.frame directly. In this example, each alternative is
# defined by travel times and toll costs.
TT1   <- choicedata$tt1
TT2   <- choicedata$tt2
TOLL2 <- choicedata$toll2

# Similarly for the vectors of choices. Note that in this example, there are only two 
# alternatives. Dummying coding the choice vector is not necessary but again makes
# writing the likelihood function straightforward.
choice1 <- (choicedata$Choice==1)
choice2 <- (choicedata$Choice==2)

# Frequency of choices
table(choicedata$Choice)



###################################################
### code chunk number 2: RSGHB_HowTo.rnw:83-102
###################################################

# The likelihood function
likelihood <- function(fc, b) {

  # Assign fixed parameters to named variables for convenience    
  cc <- 1
  ASC1  <- fc[cc]; cc <- cc + 1
  Btime <- fc[cc]; cc <- cc + 1
  Btoll <- fc[cc]; cc <- cc + 1

  # Utility functions
  v1 <- ASC1 + Btime * TT1                   
  v2 <-        Btime * TT2 + Btoll * TOLL2   
 
  # MNL probability statement
  p  <- (exp(v1) * choice1 + exp(v2) * choice2) / (exp(v1) + exp(v2))
	
  return(p)
}


###################################################
### code chunk number 3: RSGHB_HowTo.rnw:109-146
###################################################

### Setting control list for estimation (see ?doHB for more estimation options)

# The model name/description
modelname <- "MNL"          

# gVarNamesFixed contains the names for the fixed (non-random) variables in the model
# These will be used in the output and also when displaying iteration detail to 
# the screen
gVarNamesFixed <- c("ASC1", "BTime", "BCost")

# FC contains the starting values for the fixed coefficients
FC <- c(0, 0, 0)                 

# gNCREP contains the number of iterations to use prior to convergence
gNCREP <- 2500
# gNEREP contains the number of iterations to keep for averaging after convergence 
# has been reached
gNEREP <- 2500
# gNSKIP contains the number of iterations between retaining draws for averaging
gNSKIP <- 1		       
# gINFOSKIP controls how frequently to print info about the iteration process
gINFOSKIP <- 10
# gSeed ensures reproducible results
gSeed <- 1987

# To simplify the doHB function call, all of the control arguments are placed in 
# a single list that can be passed directly to doHB
control <- list(modelname = modelname,
                gVarNamesFixed = gVarNamesFixed,
                FC = FC,
                gNCREP = gNCREP,
                gNEREP = gNEREP,
                gNSKIP = gNSKIP,
                gINFOSKIP = gINFOSKIP,
                gSeed = gSeed)



###################################################
### code chunk number 4: RSGHB_HowTo.rnw:157-163
###################################################

# Estimate model
control$verbose <- FALSE
control$nodiagnostics <- TRUE
model <- doHB(likelihood, choicedata, control)



###################################################
### code chunk number 5: RSGHB_HowTo.rnw:185-189
###################################################

# Model iteration details
head(model[["iter.detail"]])



###################################################
### code chunk number 6: RSGHB_HowTo.rnw:194-198
###################################################

# Plot model statistics
plot(model)



###################################################
### code chunk number 7: RSGHB_HowTo.rnw:203-207
###################################################

# Plot parameter estimates (see ?plot.RSGHB for more uses)
plot(model, type = "F")



###################################################
### code chunk number 8: RSGHB_HowTo.rnw:222-243
###################################################

likelihood <- function(fc, b) {
     
     # Note that the fc argument is still supplied, but is unused
     
     # Using b instead of fc is the only change
     cc     <- 1
     ASC1   <- b[, cc]; cc <- cc + 1
     Btime  <- b[, cc]; cc <- cc + 1
     Btoll  <- b[, cc]; cc <- cc + 1  
  
     # Utility functions
     v1 <- ASC1 + Btime * TT1                   
     v2 <-        Btime * TT2 + Btoll * TOLL2   
 
     # MNL probability statement
     p  <- (exp(v1) * choice1 + exp(v2) * choice2) / (exp(v1) + exp(v2))
     
     return(p)
}



###################################################
### code chunk number 9: RSGHB_HowTo.rnw:250-301
###################################################
### Setting control list for estimation (see ?doHB for more estimation options)

# The model name/description
modelname <- "MMNL"

# gVarNamesNormal provides names for the random parameters in the same way
# gVarNamesFixed does for the fixed parameters
gVarNamesNormal <- c("ASC1","BTime","BCost")

# svN contains the starting values for the means of the normal distributions for each 
# of the random parameters
svN <- c(0, 0, 0)

# gDIST specifies the type of continuous distribution to use for the random parameters
# gDIST must have an entry for each value in gVarNamesNormal
# The options are:
# 1. normal
# 2. log-nomal
# 3. negative log-normal
# 4. normal with all values below zero massed at zero
# 5. normal with all values greater than zero massed at zero
# 6. Johnson SB with a specified min and max

# In this example, normal distributions are used for all 3 parameters
gDIST <- c(1, 1, 1)

# gNCREP contains the number of iterations to use prior to convergence
gNCREP <- 2500
# gNEREP contains the number of iterations to keep for averaging after convergence 
# has been reached
gNEREP <- 2500
# gNSKIP contains the number of iterations between retaining draws for averaging
gNSKIP <- 1		       
# gINFOSKIP controls how frequently to print info about the iteration process
gINFOSKIP <- 10
# gSeed ensures reproducible results
gSeed <- 1987

# To simplify the doHB function call, all of the control parameters are placed in
# a single list that can be passed directly to doHB
control <- list(modelname = modelname,
                gVarNamesNormal = gVarNamesNormal,
                gDIST = gDIST,
                svN = svN,
                gNCREP = gNCREP,
                gNEREP = gNEREP,
                gNSKIP = gNSKIP,
                gINFOSKIP = gINFOSKIP,
                gSeed = gSeed)




###################################################
### code chunk number 10: RSGHB_HowTo.rnw:312-318
###################################################

# Estimate model
control$verbose <- FALSE
control$nodiagnostics <- TRUE
model <- doHB(likelihood, choicedata, control)



###################################################
### code chunk number 11: RSGHB_HowTo.rnw:336-343
###################################################

# Model iteration details
head(model[["iter.detail"]])

# Plot model statistics
plot(model)



###################################################
### code chunk number 12: RSGHB_HowTo.rnw:348-355
###################################################

# Sample-level means
head(model[["A"]])

# Plot sample-level means
plot(model, type = "A")



###################################################
### code chunk number 13: RSGHB_HowTo.rnw:360-367
###################################################

# Average individual-level draws
head(model[["B"]])

# Standard deviations of individual-level draws
head(model[["Bsd"]])



###################################################
### code chunk number 14: RSGHB_HowTo.rnw:374-381
###################################################

# Average individual-level draws (transformed; if applicable)
head(model[["C"]])

# Standard deviations of individual-level draws (transformed; if applicable)
head(model[["Csd"]])



