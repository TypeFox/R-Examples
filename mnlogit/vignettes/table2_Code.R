library(mnlogit)
library(mlogit)
source("simChoiceModel.R") #  Has makeModel() to generate simulated data

############################################################################
### Code for producing results in Table 2 of the paper
############################################################################

## Parameter Setting

# List of number of alternatives to keep in test multinomial logit models
numChoicesVec <- c(3, 5, 7)

## Following setting used to generate Table 2 in paper 
## May take about 4 hours to run - depending on the machine 
## numChoicesVec <- c(10, 20, 30)

if (length(numChoicesVec) < 1) stop("Can't be empty!")
if (any(numChoicesVec <= 2)) 
  stop("Can't specify any model with fewer than 2 choices!")

############################################################################
## Generate simulated data and model formulas
############################################################################
# Type 'X' problems
# Make formula 
vars <- paste("X", 1:50, sep="", collapse=" + ")
fmX <- formula(paste("response ~ 1|", vars, " - 1 | 1"))

# Make formula for type 'Y' problems 
vars <- paste("X", 1:50, sep="", collapse=" + ")
fmY <- formula(paste("response ~ 1| - 1 | ", vars))

# Formula for type 'Z' problems
vars <- paste("X", 1:50, sep="", collapse=" + ")
fmZ <- formula(paste("response ~ ", vars, "| - 1 | 1"))

# Formula for type 'YZ' problems
# 5 variables of type 'Z' and 45 variables of type 'Y'
vars <- paste("X", 1:45, sep="", collapse=" + ")
fmYZ <- formula(paste("response ~ X46 + X47 + X48 + X49 + X50| - 1 | ", vars))

############################################################################
## Run code and generate results
############################################################################

# Tables for recording runnign times ratio wrt to mnlogit
# See Table 2 of the paper
# Newton-Raphson ratios
nr.ratio.table <- matrix(rep(0, 4 * length(numChoicesVec)), nrow = 4)
colnames(nr.ratio.table) <- numChoicesVec
rownames(nr.ratio.table)<-c("Problem X", "Problem Y", "Problem YZ","Problem Z")
# BFGS ratios
bfgs.ratio.table <- matrix(rep(0, 4 * length(numChoicesVec)), nrow = 4)
colnames(bfgs.ratio.table) <- numChoicesVec
rownames(bfgs.ratio.table)<-c("Problem X","Problem Y","Problem YZ","Problem Z")

for (i in 1:length(numChoicesVec)) {
  numChoices <- numChoicesVec[i]

  cat(paste0("\nRunning speed test with ", numChoices, " alternatives.\n"))
  ## Problem X
  dataX <- makeModel('X', K=numChoices) # generate data 
  # Default args set: p = 50 variables, N = K * p * 20 observations
  # Data for mlogit
  mdatX <- mlogit.data(dataX[order(dataX$indivID), ], "response", shape="long",
    alt.var="choices")
  # Run mnlogit on 1 proc
  mnlogit.time <- system.time(fit.mnlogit <- mnlogit(fmX, dataX, "choices"))[3] 
  # Run mlogit
  nr.time <- system.time(fit.mlogit <- mlogit(fmX, mdatX))[3] # Newton-Raphson
  bfgs.time <- system.time(fit.mlogit <- mlogit(fmX, mdatX, method='bfgs'))[3]
  # Compute & Store ratios
  nr.ratio.table[1, i] <- nr.time / mnlogit.time
  bfgs.ratio.table[1, i] <- bfgs.time / mnlogit.time
  cat(paste0("\n\tDone problem X.\n"))

  # Data for Type 'Y, 'Z' & 'YZ' problems
  data <- makeModel('Y', K=numChoices)  # generate data 
  # Default args set: p = 50 variables, N = K * p * 20 observations
  mdat <- mlogit.data(data[order(data$indivID), ], "response", shape="long",
    alt.var="choices")

  ## Problem Y
  # Run mnlogit on 1 proc
  mnlogit.time <- system.time(fit.mnlogit <- mnlogit(fmY, data, "choices"))[3] 
  # Run mlogit
  nr.time <- system.time(fit.mlogit <- mlogit(fmY, mdat))[3] # Newton-Raphson
  bfgs.time <- system.time(fit.mlogit <- mlogit(fmY, mdat, method='bfgs'))[3]
  # Compute & Store ratios
  nr.ratio.table[2, i] <- nr.time / mnlogit.time
  bfgs.ratio.table[2, i] <- bfgs.time / mnlogit.time
  cat(paste0("\n\tDone problem Y.\n"))

  ## Problem YZ
  # Run mnlogit on 1 proc
  mnlogit.time <- system.time(fit.mnlogit <- mnlogit(fmYZ, data, "choices"))[3] 
  # Run mlogit
  nr.time <- system.time(fit.mlogit <- mlogit(fmYZ, mdat))[3] # Newton-Raphson
  bfgs.time <- system.time(fit.mlogit <- mlogit(fmYZ, mdat, method='bfgs'))[3]
  # Compute & Store ratios
  nr.ratio.table[3, i] <- nr.time / mnlogit.time
  bfgs.ratio.table[3, i] <- bfgs.time / mnlogit.time
  cat(paste0("\n\tDone problem YZ.\n"))

  ## Problem Z
  # Run mnlogit on 1 proc
  mnlogit.time <- system.time(fit.mnlogit <- mnlogit(fmZ, data, "choices"))[3] 
  # Run mlogit
  nr.time <- system.time(fit.mlogit <- mlogit(fmZ, mdat))[3] # Newton-Raphson
  bfgs.time <- system.time(fit.mlogit <- mlogit(fmZ, mdat, method='bfgs'))[3]
  # Compute & Store ratios
  nr.ratio.table[4, i] <- nr.time / mnlogit.time
  bfgs.ratio.table[4, i] <- bfgs.time / mnlogit.time
  cat(paste0("\n\tDone problem Z.\n"))
}

cat("\nFinished running tests, printing runtime ratios (See Table 2 of paper).")
cat("\nNewton-Raphson (mlogit / mnlogit)\n")
print(nr.ratio.table)
cat("\nBFGS(mlogit / mnlogit)\n")
print(bfgs.ratio.table)
