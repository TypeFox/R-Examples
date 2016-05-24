############################################################
#                                                          #
# Analysis of the Steelhead data allowing for all back
# 60 radio tagged fish were released but only 40 passed through
# the tagging site.
#
# For an example of the analysis of this data, read the report at
# http://people.stat.sfu.ca/~cschwarz/Consulting/Moricetown/Report-2011-06-01.pdf
#                                                          #
############################################################
# 2014-09-01 CJS fix up for jags

if(file.access("demo-TSPNDENP-fall-back")!=0){ 
  dir.create("demo-TSPNDENP-fall-back", showWarnings=TRUE)}  # Test and then create the directory
setwd("demo-TSPNDENP-fall-back")

## Load the BTSPAS library
library("BTSPAS")

## Settings

## Output file prefix and title
FB.prefix <- "FB-"
FB.title  <- "FB-demo"

FB.csv <- textConnection("
jweek,marked,0 , 1 , 2 ,3 ,4 ,5 ,6 ,7     ,unmarked
29 ,  1 ,    0 , 0 , 0 ,0 ,0 ,0 ,0 ,0     ,  2
30 , 35 ,    0 , 5 , 7 ,2 ,0 ,0 ,0 ,0     ,  65
31 ,186 ,    1 ,35 ,11 ,4 ,0 ,0 ,0 ,0     ,325
32 ,292 ,    9 ,33 ,16 ,6 ,0 ,0 ,0 ,0     ,873
33 ,460 ,    6 ,41 ,16 ,9 ,3 ,0 ,2 ,1     ,976
34 ,397 ,    4 ,44 , 7 ,5 ,1 ,1 ,0 ,1     ,761
35 ,492 ,    7 ,31 ,12 ,1 ,4 ,1 ,1 ,0     ,869
36 ,151 ,    3 , 6 , 2 ,1 ,1 ,0 ,0 ,0     ,473
37 ,130 ,    3 , 2 , 2 ,0 ,0 ,1 ,0 ,0     ,332
38 ,557 ,    8 ,27 ,11 ,2 ,5 ,0 ,0 ,0     ,197
39 , 46 ,    0 , 7 , 0 ,0 ,0 ,0 ,0 ,0     ,177
40 ,143 ,   14 , 6 , 3 ,0 ,0 ,0 ,0 ,0     ,282
41 , 26 ,    2 , 1 , 0 ,0 ,0 ,0 ,0 ,0     , 82
42 ,  0 ,    0 , 0 , 0 ,0 ,0 ,0 ,0 ,0     ,100")

# Read data
FB.data <- read.csv(FB.csv,header=TRUE, as.is=TRUE, strip.white=TRUE)

## Extract components from data
## Total number of marking and recovery and recovery strata
FB.st.mark  <- nrow(FB.data)
FB.st.recov <- nrow(FB.data)

## Julian week of release/recapture
FB.jweek <- FB.data$jweek

## Number of fish marked
FB.n1 <- FB.data[1:FB.st.mark,"marked"]
FB.n1[which(FB.n1==0)] <- 1          # Correction for strata with 0 marks released

## Matrix of recoveries
FB.m2 <- as.matrix(FB.data[1:FB.st.mark, paste("X",0:(ncol(FB.data)-4),sep="")])

## Number of unmarked fish captured 
FB.u2 <- FB.data[,"unmarked"]

## Sampling fraction for each stratum   ****** THIS NEEDES TO BE REVISED
FB.sampfrac <- rep(1,FB.st.recov)

## Identify any expected jumps in abundance
FB.jump.after <- NULL 

## Identify spurious values in n1, m2, and u2 that should be removed
FB.bad.n1     <- c()     # list sample times of bad n1 values
FB.bad.m2     <- c()     # list sample times of bad m2 values
FB.bad.u2     <- c()     # list sample times of bad u2 values

## Fix capture probabilities for strata when traps not operated
FB.logitP.fixed <- NULL
FB.logitP.fixed.values <- rep(-10,length(FB.logitP.fixed))

##### Run Model #####
FB.results <- TimeStratPetersenNonDiagErrorNPMarkAvail_fit(
                  title=      FB.title,
                  prefix=     FB.prefix,
                  time=       FB.jweek,
                  n1=         FB.n1, 
                  m2=         FB.m2, 
                  u2=         FB.u2,
                  sampfrac=   FB.sampfrac,
                  jump.after= FB.jump.after,
                  bad.n1=     FB.bad.n1,
                  bad.m2=     FB.bad.m2,
                  bad.u2=     FB.bad.u2,
                  logitP.fixed=FB.logitP.fixed,
                  logitP.fixed.values=FB.logitP.fixed.values,
                  marked_available_n=66, marked_available_x=40,  # 40/66 fish did NOT fall back
	                #engine="openbugs",
                  debug=TRUE)
 
# estimate the time to target runsize and produce the plot
FB.results$TimeToTargetRunSize <- TimeToTargetRunSize( 
        U=FB.results$sims.list$U,
        time=FB.results$data$time,
        targetU=10000,  # when do 10000 unmarked fish pass the canyon?
        file_prefix = FB.prefix)

# Rename files that were created.

file.copy("data.txt",       paste(FB.prefix,".data.txt",sep=""),      overwrite=TRUE)
file.copy("CODAindex.txt",  paste(FB.prefix,".CODAindex.txt",sep=""), overwrite=TRUE)
file.copy("CODAchain1.txt", paste(FB.prefix,".CODAchain1.txt",sep=""),overwrite=TRUE)
file.copy("CODAchain2.txt", paste(FB.prefix,".CODAchain2.txt",sep=""),overwrite=TRUE)
file.copy("CODAchain3.txt", paste(FB.prefix,".CODAchain3.txt",sep=""),overwrite=TRUE)
file.copy("inits1.txt",     paste(FB.prefix,".inits1.txt",sep=""),    overwrite=TRUE)
file.copy("inits2.txt",     paste(FB.prefix,".inits2.txt",sep=""),    overwrite=TRUE)
file.copy("inits3.txt",     paste(FB.prefix,".inits3.txt",sep=""),    overwrite=TRUE)

file.remove("data.txt"       )       
file.remove("CODAindex.txt"  )
file.remove("CODAchain1.txt" )
file.remove("CODAchain2.txt" )
file.remove("CODAchain3.txt" )
file.remove("inits1.txt"     )
file.remove("inits2.txt"     )
file.remove("inits3.txt"     )
 

save(list=c("FB.results"), file="FB-results-saved.Rdata")  # save the results from this run

cat("\n\n\n ***** FILES and GRAPHS saved in \n    ", getwd(), "\n\n\n")
print(dir())

# move up the directory
setwd("..")

cat("***** END OF RUN *****", date(), "\n")

setwd("..")

