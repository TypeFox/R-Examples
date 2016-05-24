# 2011-03-15
# This is a demonstration of the error messages in calling the
#    TimeStratPetersenNonDiagErrorNP program
# 2014-09-01 CJS jags converstion; remove prompts

# library("BTSPAS")

# Get the data.
demo.data <- textConnection(
"jweek,marked,0,1,2,3,4,5,6,7,unmarked
29,1,0,0,0,0,0,0,0,0,2
30,35,0,5,7,2,0,0,0,0,65
31,186,1,35,11,4,0,0,0,0,325
32,292,9,33,16,6,0,0,0,0,873
33,460,6,41,16,9,3,0,2,1,976
34,397,4,44,7,5,1,1,0,1,761
35,492,7,31,12,1,4,1,1,0,869
36,151,3,6,2,1,1,0,0,0,473
37,130,3,2,2,0,0,1,0,0,332
38,557,8,27,11,2,5,0,0,0,197
39,46,0,7,0,0,0,0,0,0,177
40,143,14,6,3,0,0,0,0,0,282
41,26,2,1,0,0,0,0,0,0,82
42,0,0,0,0,0,0,0,0,0,100")

demo.fish <- read.csv(demo.data, header=TRUE)

demo.prefix <- "MT-2010"
demo.title  <- "Moricetown 2010 data with errors in input data"
cat("*** Starting ",demo.title, "\n\n")

# Extract the various data values
# The weekly data
demo.jweek <- demo.fish$jweek

# The number of strata
demo.st.mark  <- nrow(demo.fish)
demo.st.recov <- nrow(demo.fish)

demo.n1 <- demo.fish[1:demo.st.mark,"marked"]
demo.n1[which(demo.n1==0)] <- 1          # Correction for strata with 0 marks released

## Matrix of recoveries
demo.m2 <- as.matrix(demo.fish[1:demo.st.mark, paste("X",0:(ncol(demo.fish)-4),sep="")])

## Number of unmarked fish captured 
demo.u2 <- demo.fish[,"unmarked"]

## Sampling fraction for each stratum 
demo.sampfrac <- rep(1,demo.st.recov)

## Identify any expected jumps in abundance
demo.jump.after <- NULL 

## Identify spurious values in n1, m2, and u2 that should be removed
demo.bad.n1     <- c()     # list sample times of bad n1 values
demo.bad.m2     <- c()     # list sample times of bad m2 values
demo.bad.u2     <- c()     # list sample times of bad u2 values

## Fix capture probabilities for strata when traps not operated
demo.logitP.fixed <- NULL
demo.logitP.fixed.values <- rep(-10,length(demo.logitP.fixed))

##### Run Model with valid data #####
#demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
#                  title=      demo.title,
#                  prefix=     demo.prefix,
#                  time=       demo.jweek,
#                  n1=         demo.n1, 
#                  m2=         demo.m2, 
#                  u2=         demo.u2,
#                  sampfrac=   demo.sampfrac,
#                  jump.after= demo.jump.after,
#                  bad.n1=     demo.bad.n1,
#                  bad.m2=     demo.bad.m2,
#                  bad.u2=     demo.bad.u2,
#                  logitP.fixed=demo.logitP.fixed,
#                  logitP.fixed.values=demo.logitP.fixed.values,
#                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
#                  debug=FALSE,
#                  #engine="openbugs"
#)


cat("********** CHECKING ERROR MESSAGES IN TimeStratPetersenNonDiagErrorNP_fit *********\n\n")

# Generate the possible error messages.
# Lengths of the n1, m2, u2, sampfrac, time not all the same
demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1[-1], 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"
)

demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2[-1],
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"  # show how to call openbugs
)

demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  logitP.cov= rep(1,20),
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"    # show how to call openbugs
)


## 2. Check that rowsum of m2 <= n1
demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2*100, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"   # show how to call openbugs
)


## 3. Check that elements of bad.m2 etc belong to time

demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     max(demo.jweek)+1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"     # show how to call openbugs
)

demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     max(demo.jweek)+1,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"   # show how to call openbugs
)

demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     max(demo.jweek)+1,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"   # show how to call openbugs
)

demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= max(demo.jweek)+1,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"  # show how to call openbugs
)

## Check for existence of openbugs/winbugs directory
demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  engine="openbugs",
                  OPENBUGS.directory="SLFDJS",
) 


## 5. Check that index of logitP.fixed is ok

demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=max(demo.jweek)+1,
                  logitP.fixed.values=c(-10),
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"  # show how to call openbugs
)


## 6. Check that length of prior for muTT is correct length
demo.results <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jweek,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed,
                  logitP.fixed.values=demo.logitP.fixed.values,
                  prior.muTT = c(23), 
                  n.chains=3, n.iter=200, n.burnin=100, n.sims=20,
                  debug=FALSE,
                  #engine="openbugs"  # show how to call openbugs
)


cat("\n\n\n ***** End of Demonstration *****\n\n\n")


