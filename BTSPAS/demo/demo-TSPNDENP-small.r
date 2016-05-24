# 2014-09-01 CJS updated for jages; removed prompts
# 2010-03-21 Conne River 2009 analysis.
# 2010-11-29 Updated demo to show how to estimate time to get target value

# Demonstration of the Time-Stratified Petersen experiment with
# NON-Diagonal entries, travel times for the marked fish modelled
# non-parametrically, and some capture probabilities fixed to 0
# because no sampling occurred.


# Create a directory to store the results Test and then create the
# directory
if(file.access("demo-TSPNDENP-small")!=0){
  dir.create("demo-TSPNDENP-small", showWarnings=TRUE)
}
setwd("demo-TSPNDENP-small")

## Load BTSPAS library
library(BTSPAS)

## Define data
demo.data <- textConnection(
"Date,tagged,   0, 1, 2, 3,   Tot-recoveries,Untagged,Tot-caught,WaterTemp
2009-04-29, 25, 1, 2, 1, 0,     4,   0,   0,   7.8
2009-04-30, 75, 3, 8, 0, 0,    19, 133, 133,   7.0
2009-05-01, 97, 4,16, 0, 0,    20, 158, 161,   7.5
2009-05-02,127, 6, 0, 0, 3,    9,  128, 130,   7.7
2009-05-03,  0, 0, 0, 0, 0,     0,   0,   0,   9.3
2009-05-04,  0, 0, 0, 0, 0,     0,   0,   0,  11.7
2009-05-05,216, 8,12, 4, 1,    25,  859, 893,  11.6
2009-05-06,215, 1,26, 2, 0,    29,  427, 445,  11.7
2009-05-07,205,11,13,11, 0,    35,  849, 892,  11.1
2009-05-08,202, 3,23, 2, 0,    28,  488, 507,   9.8
2009-05-09  ,0, 0, 0, 0, 0,     0,   22,  22,  10.9
2009-05-10,  0, 0, 0, 0, 0,     0,   22,  22,  12.2
2009-05-11,  0, 0, 0, 0, 0,     0,    1,   1, 14.3")

demo.Fish <- read.csv(demo.data, header=TRUE, as.is=TRUE, strip.white=TRUE)

# Now to extract the subset of data, do any fancy adjustments, and fit the data.

demo.prefix <- "Small-TSPNDENP"
demo.title  <- "Small NP"
cat("*** Starting ",demo.title, "\n\n")

# We do some fixups because the data entry forces you to enter extra rows of zeros
# at the end because the u2 vector can be longer than the n1 vector.
demo.n1 <- demo.Fish$tagged
demo.n1 <- demo.n1[ 1:(length(demo.n1)-3)]  # last rows have no fish released, but have fish recaptured

# Notice that OpenBugs/WinBugs has problems with binomial/multinomial distributions with an index of 0.
# Consequently, we set the number of fish released to 1. This "approximation" will work fine in most
# real situations.
demo.n1[c(5,6)] <- 1  # no fish released, or a blow out on recoveries


demo.m2 <- demo.Fish[, paste("X",0:3,sep="")]
demo.m2 <- as.matrix(demo.m2)
demo.m2 <- demo.m2[ 1:length(demo.n1),]     # last rows have no fish released

demo.u2 <- demo.Fish$Untagged

# what fraction of the day was sampled?
demo.sampfrac <- rep(1,length(demo.u2)) # values are on a daily basis

# what is the strata identification number (julian day since start of year)?
demo.jday <- 119+ 1:length(demo.u2)

# are there any jumps in the abundance?
demo.jump.after <- NULL    # list sample times after which jump in number occurs

# are there any bad values that need to be adjusted?
demo.bad.n1     <- c()     # list sample times of bad n1 values
demo.bad.m2     <- c()     # list sample times of bad m2 values
demo.bad.u2     <- c()     # list sample times of bad u2 values

# are there any days where the capture probability is fixed in advance?
# On days 13 and 14 the second trap could not operate because of high water levels.
demo.logitP.fixed <- 119+ c(5,6)
demo.logitP.fixed.values <- rep(-10,length(demo.logitP.fixed))  # logitP=-10 === p[i]=0

## Run TSPNDE model
demo.cr.2009.as.tspndenp <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.jday,
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
		              #engine="openbugs",
                  debug=TRUE
                  )
# Rename files that were created.

file.rename("data.txt",       paste(demo.prefix,".data.txt",sep=""))
file.rename("CODAindex.txt",  paste(demo.prefix,".CODAindex.txt",sep=""))
file.rename("CODAchain1.txt", paste(demo.prefix,".CODAchain1.txt",sep=""))
file.rename("CODAchain2.txt", paste(demo.prefix,".CODAchain2.txt",sep=""))
file.rename("CODAchain3.txt", paste(demo.prefix,".CODAchain3.txt",sep=""))
file.rename("inits1.txt",     paste(demo.prefix,".inits1.txt",sep=""))
file.rename("inits2.txt",     paste(demo.prefix,".inits2.txt",sep=""))
file.rename("inits3.txt",     paste(demo.prefix,".inits3.txt",sep=""))

# Extract the total abundance over time, and extract the time needed to reach a target value

demo.cr.2009.as.tspndenp$TimeToTargetRunSize <- TimeToTargetRunSize(
        U=demo.cr.2009.as.tspndenp$sims.list$U,
        time=demo.cr.2009.as.tspndenp$data$time,
        targetU=30000,
        file_prefix = demo.prefix)


# Save the information for later retreival if needed

save(list=c("demo.cr.2009.as.tspndenp"), file="Small-tspndenp-saved.Rdata")  # save the results from this run

# move up the directory
setwd("..")

cat("\n\n\n ***** End of Demonstration *****\n\n\n")
