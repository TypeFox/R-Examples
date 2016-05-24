# 2014-09-01 CJS modified for JAGS; removed prompts
# 2010-03-03 Conne River 2009 analysis.

# Demonstration of the Time-Stratified Petersen experiment with NON-Diagonal entries
# with some capture probabilities fixed to 0 because no sampling occurred.

# Uses the Conne River 2009 data as provided by Brian Dempson (DFO, St. John's) on 2010-02-08.

# There are two traps on the Conne River.
# At the first trap, fish are captured and marked with individually numbered tags on a
# daily basis. They are released, swim downstream, and  may be recaptured at a secondary 
# trap. The number of recaptured fish and the number of unmarked fish at the second trap
# are also recorded on a daily basis.
#
# In days 13 and 14, the water levels were too high so no fish could be released
# nor unmarked fish captured. Consequently, the logitP is set to -10 (corresponding to p[i]=0) for these days

#
# The vector n1[i] is the number of fish marked and released on day [i] at the first trap.
#
# The matrix m2[i,j] is the number of fish released on day [i] and recaptured j-1 days later
# at the second trap, i.e. the first column are the number of fish recaptured 
# the same day, the second column the number of fish
# recaptured the day after release etc. 
# All marked fish are assumed to have passed the second trap by day i+(ncol(m2)-1) days later.
#
# The vector u2[j] is the number of unmarked fish captured at the second trap on day [j].
# 
# As in the diagonal case, "bad" values of n1, m2, and u2 are allowed to be removed.
# This is tricker for the m2 matrix as only a single entry may be in error.
# 
# For more details refer to:
#   Schwarz, C.J., and Dempson, B.D. (1994). 
#   Mark-recapture estimation of a salmon smolt population. 
#   Biometrics 50: 98-108.
#
# In this implementation, the parameters of the log-normal movement distribution are
# smoothed over time.

# Create a directory to store the results
if(file.access("demo-TSPNDE")!=0){ dir.create("demo-TSPNDE", showWarnings=TRUE)}  # Test and then create the directory
setwd("demo-TSPNDE")

library("BTSPAS")

# Get the raw data and read it in
demo.data <- textConnection(
"Date ,   tagged , 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 ,  Tot-recoveries , Untagged, Dummy1, Dummy2, WaterTemp
2009-04-29 , 25  , 0 , 0 , 1 , 0 , 0 , 1 , 0 , 0 ,   2              , 0 , 0 ,  , 7.8
2009-04-30 , 75  , 0 , 2 , 2 , 2 , 2 , 1 , 0 , 0 ,   9              , 133 , 133 ,  , 7.0
2009-05-01 , 97  , 0 , 0 , 1 , 4 , 1 , 0 , 0 , 1 ,   7              , 158 , 161 ,  , 7.5
2009-05-02 , 127 , 0 , 2 , 5 , 0 , 0 , 1 , 0 , 0 ,   8              , 128 , 130 ,  , 7.7
2009-05-03 , 194 , 0 , 6 , 2 , 1 , 0 , 0 , 0 , 0 ,   9              , 197 , 202 ,  , 9.3
2009-05-04 , 200 , 2 , 22 , 4 , 1 , 0 , 0 , 0 , 0 ,  29             , 522 , 542 ,  , 11.7
2009-05-05 , 216 , 8 , 12 , 4 , 1 , 1 , 0 , 0 , 0 ,  26             , 859 , 893 ,  , 11.6
2009-05-06 , 215 , 1 , 26 , 2 , 0 , 1 , 0 , 0 , 1 ,  31             , 427 , 445 ,  , 11.7
2009-05-07 , 205 , 11 , 13 , 11 , 0 , 0 , 0 , 0 , 0, 35             , 849 , 892 ,  , 11.1
2009-05-08 , 202 , 3 , 23 , 2 , 0 , 0 , 0 , 0 , 0  , 28             , 488 , 507 ,  , 9.8
2009-05-09 , 222 , 7 , 1 , 1 , 0 , 0 , 0 , 0 , 0 ,    9             , 1080 , 1123 ,  , 9.5
2009-05-10 , 158 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 350 , 354 ,  , 9.2
2009-05-11 , 0   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 0 , 0 ,  , 6.6
2009-05-12 , 0   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 0 , 0 ,  , 7.4
2009-05-13 , 161 , 0 , 12 , 6 , 2 , 0 , 0 , 0 , 0 ,  20             , 58 , 59 ,  , 9.6
2009-05-14 , 198 , 3 , 12 , 1 , 0 , 0 , 0 , 1 , 0 ,  17             , 336 , 351 ,  , 11.4
2009-05-15 , 61  , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,    1             , 118 , 137 ,  , 11.6
2009-05-16 , 66  , 0 , 1 , 4 , 1 , 0 , 0 , 0 , 0 ,    6             , 72 , 75 ,  , 12.3
2009-05-17 , 31  , 0 , 0 , 3 , 0 , 0 , 0 , 0 , 0 ,    3             , 34 , 36 ,  , 12.9
2009-05-18 , 11  , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 26 , 30 ,  , 13.9
2009-05-19 , 12  , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 21 , 25 ,  , 14.4
2009-05-20 , 5   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 16 , 16 ,  , 13.6
2009-05-21 , 0   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 22 , 22 ,  , 10.9
2009-05-22 , 0   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 22 , 22 ,  , 12.2
2009-05-23 , 0   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 1 , 1 ,  , 14.3
2009-05-24 , 0   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 1 , 1 ,  , 13.0
2009-05-25 , 0   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 5 , 6 ,  , 13.9
2009-05-26 , 0   , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,    0             , 3 , 3 ,  , 13.2     ")

demo.Fish <- read.csv(demo.data, header=TRUE)

demo.Fish[1:5,]

# Now to extract the subset of data, do any fancy adjustments, and fit the data.

demo.prefix <- "CR-2009-AS-TSPNDE"
demo.title  <- "Conne River 2009 Atlantic Salmon "
cat("*** Starting ",demo.title, "\n\n")

# We do some fixups because the data entry forces you to enter extra rows of zeros
# at the end because the u2 vector can be longer than the n1 vector.
demo.n1 <- demo.Fish$tagged
demo.n1 <- demo.n1[ 1:(length(demo.n1)-6)]  # last rows have no fish released, but have fish recaptured

# Notice that OpenBugs has problems with binomial/multinomial distributions with an index of 0.
# Consequently, we set the number of fish released to 1. This "approximation" will work fine in most
# real situations.
demo.n1[c(13,14)] <- 1  # no fish released, or a blow out on recoveries

demo.m2 <- demo.Fish[, paste("X",0:7,sep="")]
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
demo.logitP.fixed <- 119+ c(13,14)
demo.logitP.fixed.values <- rep(-10,length(demo.logitP.fixed))  # logitP=-10 === p[i]=0


demo.cr.2009.as.tspnde <- TimeStratPetersenNonDiagError_fit(
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
		  #engine="openbugs",   # use the openbugs engine
                  debug=TRUE)

# Rename files that were created.

file.copy("data.txt",       paste(demo.prefix,".data.txt",sep=""),      overwrite=TRUE)
file.copy("CODAindex.txt",  paste(demo.prefix,".CODAindex.txt",sep=""), overwrite=TRUE)
file.copy("CODAchain1.txt", paste(demo.prefix,".CODAchain1.txt",sep=""),overwrite=TRUE)
file.copy("CODAchain2.txt", paste(demo.prefix,".CODAchain2.txt",sep=""),overwrite=TRUE)
file.copy("CODAchain3.txt", paste(demo.prefix,".CODAchain3.txt",sep=""),overwrite=TRUE)
file.copy("inits1.txt",     paste(demo.prefix,".inits1.txt",sep=""),    overwrite=TRUE)
file.copy("inits2.txt",     paste(demo.prefix,".inits2.txt",sep=""),    overwrite=TRUE)
file.copy("inits3.txt",     paste(demo.prefix,".inits3.txt",sep=""),    overwrite=TRUE)

file.remove("data.txt"       )       
file.remove("CODAindex.txt"  )
file.remove("CODAchain1.txt" )
file.remove("CODAchain2.txt" )
file.remove("CODAchain3.txt" )
file.remove("inits1.txt"     )
file.remove("inits2.txt"     )
file.remove("inits3.txt"     )
 

save(list=c("demo.cr.2009.as.tspnde"), file="demo.cr-2009-as-tspnde-saved.Rdata")  # save the results from this run

cat("\n\n\n ***** FILES and GRAPHS saved in \n    ", getwd(), "\n\n\n")
print(dir())

# move up the directory
setwd("..")

cat("\n\n\n ***** End of Demonstration *****\n\n\n")

