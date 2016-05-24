# 2014-09-01 CJS Removed prompts; showed how to call openbugs
# 2009-12-05
# This is a demonstration of the error messages in calling Time Stratified Petersen with Diagonal Entries (TSPDE) program

# load library
library("BTSPAS")

# Get the data. In many cases, this is stored in a *.csv file and read into the program
# using a read.csv() call. In this demo, the raw data is assigned directly as a vector.
#

# Indicator for the week.
demo.jweek <- c(9,   10,   11,   12,   13,   14,   15,   16,   17,   18,
          19,   20,   21,   22,   23,   24,   25,   26,   27,   28, 
          29,   30,   31,   32,   33,   34,   35,   36,   37,   38,
          39,   40,   41,   42,   43,   44,   45,   46)

# Number of marked fish released in each week.
demo.n1 <- c(   0, 1465, 1106,  229,   20,  177,  702,  633, 1370,  283,
         647,  276,  277,  333, 3981, 3988, 2889, 3119, 2478, 1292,
        2326, 2528, 2338, 1012,  729,  333,  269,   77,   62,   26,
          20, 4757, 2876, 3989, 1755, 1527,  485,  115)

# Number of marked fish recaptured in the week of release. No marked fish
# are recaptured outside the week of release.
demo.m2 <- c(   0,   51,  121,   25,    0,   17,   74,   94,   62,   10,
          32,   11,   13,   15,  242,   55,  115,  198,   80,   71, 
         153,  156,  275,  101,   66,   44,   33,    7,    9,    3,
           1,  188,    8,   81,   27,   30,   14,    4)

# Number of unmarked fish captured at the trap in each week.
demo.u2 <- c(4135,10452, 2199,  655,  308,  719,  973,  972, 2386,  469,
         897,  426,  407,  526,39969,17580, 7928, 6918, 3578, 1713, 
        4212, 5037, 3315, 1300,  989,  444,  339,  107,   79,   41,
          23,35118,34534,14960, 3643, 1811,  679,  154)

# What fraction of the week was sampled?
demo.sampfrac<-c(3,   8,    6,    7,    7,    7,    7,    7,    7,    7,
            7,   7,    7,    7,    7,    7,    7,    7,    7,    7,
            6,   7,    7,    7,    7,    7,    7,    7,    7,    7,
            7,   7,    7,    7,    7,    7,    7,    5)/7

# After which weeks is the spline allowed to jump?
demo.jump.after <- c(22,39)  # julian weeks after which jump occurs

# Which julian weeks have "bad" recapture values. These will be set to missing and estimated.
demo.bad.m2     <- c(41)   # list julian week with bad m2 values

# The prefix for the output files:
demo.prefix <- "demo-JC-2003-CH-TSPDE" 

# Title for the analysis
demo.title <- "Junction City 2003 Chinook "




cat("********** CHECKING ERROR MESSAGES IN TSPDE *********\n\n")

# Generate the possible error messages.
# Lengths of the n1, m2, u2, sampfrac, time not all the same
demo.jweek_error <- demo.jweek[-1]
demo.jc.2003.ch.tspde <- TimeStratPetersenDiagError_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek_error,
                  n1=demo.n1, 
                  m2=demo.m2, 
                  u2=demo.u2,
                  sampfrac=demo.sampfrac,
                  jump.after=demo.jump.after,
                  bad.m2=demo.bad.m2,
		  #engine="openbugs",
                  debug=TRUE  # this generates only 10,000 iterations of the MCMC chain for checking.
                  )
demo.logitP.cov_error <- matrix(0,nrow=length(demo.n1)+2,ncol=2)
demo.jc.2003.ch.tspde <- TimeStratPetersenDiagError_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek,
                  n1=demo.n1, 
                  m2=demo.m2, 
                  u2=demo.u2,
                  logitP.cov=demo.logitP.cov_error,
                  sampfrac=demo.sampfrac,
                  jump.after=demo.jump.after,
                  bad.m2=demo.bad.m2,
		  #engine="openbugs",
                  debug=TRUE  # this generates only 10,000 iterations of the MCMC chain for checking.
                  )

# Make sure that all m2 <= n1
demo.m2_error <- demo.m2
demo.m2_error[1] <- demo.n1[1]+5
demo.jc.2003.ch.tspde <- TimeStratPetersenDiagError_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek,
                  n1=demo.n1, 
                  m2=demo.m2_error, 
                  u2=demo.u2,
                  sampfrac=demo.sampfrac,
                  jump.after=demo.jump.after,
                  bad.m2=demo.bad.m2,
		  #engine="openbugs",
                  debug=TRUE  # this generates only 10,000 iterations of the MCMC chain for checking.
                  )


# check that elements of bad.m2 belong to time;
demo.bad.m2_error <- demo.bad.m2
demo.bad.m2_error <- c(demo.bad.m2_error, max(demo.jweek)+1)
demo.jc.2003.ch.tspde <- TimeStratPetersenDiagError_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek,
                  n1=demo.n1, 
                  m2=demo.m2,
                  u2=demo.u2,
                  sampfrac=demo.sampfrac,
                  jump.after=demo.jump.after,
                  bad.m2=demo.bad.m2_error,
		  #engine="openbugs",
                  debug=TRUE  # this generates only 10,000 iterations of the MCMC chain for checking.
                  )


# check that elements of jump.after belong to time;
demo.jump.after_error <- demo.jump.after
demo.jump.after_error <- c(demo.jump.after_error, max(demo.jweek)+1)
demo.jc.2003.ch.tspde <- TimeStratPetersenDiagError_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek,
                  n1=demo.n1, 
                  m2=demo.m2,
                  u2=demo.u2,
                  sampfrac=demo.sampfrac,
                  jump.after=demo.jump.after_error,
                  bad.m2=demo.bad.m2,
		  #engine="openbugs",
                  debug=TRUE  # this generates only 10,000 iterations of the MCMC chain for checking.
                  )

# check to see if specifying wrong path to openbugs/winbugs works;
demo.jc.2003.ch.tspde <- TimeStratPetersenDiagError_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek,
                  n1=demo.n1, 
                  m2=demo.m2,
                  u2=demo.u2,
                  sampfrac=demo.sampfrac,
                  jump.after=demo.jump.after,
                  bad.m2=demo.bad.m2,
		  #engine="openbugs",
                  debug=TRUE  # this generates only 10,000 iterations of the MCMC chain for checking.
                  )

# check to see if specifying wrong path to openbugs/winbugs works;
demo.jc.2003.ch.tspde <- TimeStratPetersenDiagError_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek,
                  n1=demo.n1, 
                  m2=demo.m2,
                  u2=demo.u2,
                  sampfrac=demo.sampfrac,
                  jump.after=demo.jump.after,
                  bad.m2=demo.bad.m2,
		  #engine="openbugs",
                  debug=TRUE  # this generates only 10,000 iterations of the MCMC chain for checking.
                  )


# move up the directory
setwd("..")


cat("\n\n\n ***** End of Demonstration *****\n\n\n")


