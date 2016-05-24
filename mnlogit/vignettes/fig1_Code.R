##############################################################################
## Code to generate vignette's Fig. 1 (Parallel Hessian Calculation Speedup) # 
##                                                                           #
## Speed tests single & multi threaded mnlogit using simulated data.         #
##                                                                           #
## Running this script will produce single and multicore runtimes on your    #
## LOCAL machine. Runtimes WILL fluctuate if other programs are running.     #
##############################################################################

#############################################################################
## 			Main user parameters 				   ##
#############################################################################
# Vector with number of cores to run tests on
cores.run <- c(1, 2, 4, 8, 16) # Seprate timing runs for each 
# Number of choices
num.choices <- 100             # Controls problem size 
# With these setting it takes about 4 hours to run this script on dev machine.
# 50 covariates is default choice set below in 'nvars' parameter & formulas.

#############################################################################
## File to plot Hessian calculation speedup data (NULL disables plotting)   #
#############################################################################
hessSpeedupFile <- paste0("speedupHess_K", num.choices, ".pdf") # NULL
#hessSpeedupFile <- "HessianSpeedups.pdf" # NULL
#############################################################################

library(mnlogit)
source("simChoiceModel.R")  # for data generation

# Runs mnlogit and returns runtimes (all in seconds)
time.mnlogit <- function(formula, data, num.cores=1)
{
    cat("\n==============\nTiming Run Starts\n==============\n")
    cat(paste0("Using ", num.cores, " processors.\n"))
    elapsed <- system.time(
        # Increase print.level to 1 to enable in-iteration printing
        fit <- mnlogit(formula, data, choiceVar="choices", ncores=num.cores,
                       maxiter=50, ftol=1e-12, gtol=1e-3, print.level=0)
        )
    return(list(tot.time  = elapsed[3],
                hess.time = fit$est.stat$hessMins * 60))
    cat("\n==============\nTiming Run Ends\n==============\n")
}

# Creates a table of runtimes (all in sec), after running over multiple cores
run.tests <- function(formula, data, num.cores.vec)
{
    if (!is.vector(num.cores.vec))
        stop("Invalid num.cores.vec argument.")

    timings <- matrix(rep(0, 3*length(num.cores.vec)), 
                      nrow=length(num.cores.vec), ncol=3)
    colnames(timings) <- c("ncores", "hess.time", "tot.time")

    for (i in c(1:length(num.cores.vec))) {
        ncores <- num.cores.vec[i]
        stats <- time.mnlogit(formula, data, num.cores=ncores)
        timings[i, 1] <- ncores
        timings[i, 2] <- stats$hess.time
        timings[i, 3] <- stats$tot.time
    }
    return(timings)
}

#############################################################################
## Formulas of the 4 problem classes. See section 3 of the vignette.
#############################################################################

# 50 variables of type 'X'
fmX <- formula(response ~ 1| X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26 + X27 + X28 + X29 + X30 + X31 + X32 + X33 + X34 + X35 + X36 + X37 + X38 + X39 + X40 + X41 + X42 + X43 + X44 + X45 + X46 + X47 + X48 + X49 + X50 - 1| 1)

# 50 variables of type 'Y'
fmY <- formula(response ~ 1| -1| X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26 + X27 + X28 + X29 + X30 + X31 + X32 + X33 + X34 + X35 + X36 + X37 + X38 + X39 + X40 + X41 + X42 + X43 + X44 + X45 + X46 + X47 + X48 + X49 + X50)

# 50 variables of type 'Z'
fmZ <- formula(response ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26 + X27 + X28 + X29 + X30 + X31 + X32 + X33 + X34 + X35 + X36 + X37 + X38 + X39 + X40 + X41 + X42 + X43 + X44 + X45 + X46 + X47 + X48 + X49 + X50| -1 | 1)

# 45 variables of type 'Y', 5 variables of type 'Z'
fmYZ <- formula(response ~ X1 + X2 + X3 + X4 + X5 | -1 | X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26 + X27 + X28 + X29 + X30 + X31 + X32 + X33 + X34 + X35 + X36 + X37 + X38 + X39 + X40 + X41 + X42 + X43 + X44 + X45 + X46 + X47 + X48 + X49 + X50)

nvars <- 50 # number of variables in each of the 4 problem classes

#############################################################################
# Run problem type 'X' timing tests
#############################################################################
cat("\n--------------\nCASE X\n--------------\n")
df <- makeModel('X', K=num.choices, numCovars=nvars)
timeX <- run.tests(fmX, df, cores.run) 
cat("\n----------------------------\n")

#############################################################################
# Run problem type 'Y' timing tests
#############################################################################
cat("\n--------------\nCASE Y\n--------------\n")
df <- makeModel('Y', K=num.choices, numCovars=nvars)
timeY <- run.tests(fmY, df, cores.run) 
cat("\n----------------------------\n")

#############################################################################
# Run problem type 'Z' timing tests
#############################################################################
cat("\n--------------\nCASE Z\n--------------\n")
df <- makeModel('Y', K=num.choices, numCovars=nvars)
timeZ <- run.tests(fmZ, df, cores.run) 
cat("\n----------------------------\n")

#############################################################################
# Run problem type 'YZ' timing tests
#############################################################################
cat("\n--------------\nCASE YZ\n--------------\n")
df <- makeModel('Y', K=num.choices, numCovars=nvars)
timeYZ <- run.tests(fmYZ, df, cores.run)
cat("\n----------------------------\n")


#############################################################################
cat("\n***** FINAL STATS *****\n")
#############################################################################
cat("\n--------------\nCASE X\n--------------\n")
print(timeX) 
cat("\n--------------\nCASE Y\n--------------\n")
print(timeY) 
cat("\n--------------\nCASE Z\n--------------\n")
print(timeZ) 
cat("\n--------------\nCASE YZ\n--------------\n")
print(timeYZ)

#############################################################################
# Hessian calculation speedups
#############################################################################
hessTimes <- cbind(timeX[ ,2], timeY[ ,2], timeZ[ ,2], timeYZ[ ,2])
hessSpeedup <- apply(hessTimes, 2, function(vec) 1.0/(vec/vec[1]))
hessSpeedup <- cbind(timeX[ ,1], hessSpeedup)
colnames(hessSpeedup) <- c("ncores", "caseX", "caseY", "caseZ", "caseYZ")
cat("\n--------------\nHessian parallel speedup\n--------------\n")
print(hessSpeedup)

#############################################################################
## Plotting hessian speedup. See Fig 1 of the vignette.
#############################################################################
if (!is.null(hessSpeedupFile)) {
  pdf(hessSpeedupFile)
  matplot(hessSpeedup[ , 1], hessSpeedup[ , 2:5], type="b", pch=1:4, col=1:4,
          xlab="procs", ylab="Speedup")
  legend("topleft", c("case X", "case Y", "case YZ", "case Z"), pch=1:4, col=1:4)
  dev.off()
  cat(paste("\nPlotted speedups to file:", hessSpeedupFile, "\n"))
}
