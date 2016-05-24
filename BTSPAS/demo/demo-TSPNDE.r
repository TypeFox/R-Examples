# 2014-09-01 CJS fixup for JAGS, removal of prompts
# 2009-12-07 

# Demonstration of the Time-Stratified Petersen experiment with NON-Diagonal entries
# This uses a simulated subset of the the original Conne River data.
#
# There are two traps on the Conne River.
# At the first trap, fish are captured and marked with individually numbered tags on a
# daily basis. They are released, swim downstream, and  may be recaptured at a secondary 
# trap. The number of recaptured fish and the number of unmarked fish at the second trap
# are also recorded on a daily basis.
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
"Date , tagged , 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , Tot-recoveries , Untagged
1987-04-26 , 8 , 0 , 0 , 0 , 0 , 2 , 0 , 0 , 0 , 0 , 0 , 2 , 0
1987-04-27 , 5 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 2
1987-04-28 , 6 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1
1987-04-29 , 17 , 0 , 0 , 2 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 4 , 2
1987-04-30 , 66 , 0 , 1 , 0 , 2 , 3 , 2 , 0 , 0 , 0 , 0 , 8 , 39
1987-05-01 , 193 , 0 , 1 , 7 , 7 , 2 , 2 , 0 , 0 , 0 , 0 , 19 , 226
1987-05-02 , 90 , 0 , 2 , 0 , 0 , 0 , 2 , 0 , 0 , 0 , 0 , 4 , 75
1987-05-03 , 260 , 0 , 0 , 14 , 6 , 1 , 1 , 1 , 0 , 1 , 0 , 24 , 129
1987-05-04 , 368 , 0 , 9 , 46 , 4 , 2 , 1 , 0 , 3 , 0 , 1 , 66 , 120
1987-05-05 , 506 , 0 , 38 , 33 , 11 , 0 , 1 , 3 , 0 , 0 , 0 , 86 , 380
1987-05-06 , 317 , 1 , 27 , 26 , 3 , 1 , 4 , 0 , 0 , 0 , 0 , 62 , 921
1987-05-07 , 43 , 0 , 4 , 3 , 0 , 2 , 0 , 0 , 0 , 0 , 0 , 9 , 1005
1987-05-08 , 259 , 1 , 42 , 5 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 50 , 1181
1987-05-09 , 259 , 1 , 32 , 27 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 61 , 1087
1987-05-10 , 249 , 1 , 85 , 3 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 90 , 1108
1987-05-11 , 250 , 3 , 21 , 19 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 45 , 1685
1987-05-12 , 298 , 42 , 16 , 11 , 9 , 1 , 0 , 0 , 0 , 0 , 0 , 79 , 671
1987-05-13 , 250 , 1 , 7 , 25 , 6 , 4 , 0 , 0 , 0 , 0 , 0 , 43 , 1766
1987-05-14 , 193 , 0 , 9 , 18 , 8 , 0 , 0 , 0 , 0 , 0 , 0 , 35 , 636
1987-05-15 , 207 , 0 , 17 , 21 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 40 , 483
1987-05-16 , 175 , 0 , 18 , 10 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 29 , 170
1987-05-17 , 141 , 0 , 12 , 14 , 7 , 1 , 1 , 0 , 0 , 0 , 0 , 35 , 269
1987-05-18 , 155 , 0 , 1 , 19 , 13 , 6 , 2 , 0 , 0 , 0 , 0 , 41 , 212
1987-05-19 , 123 , 0 , 5 , 22 , 5 , 0 , 0 , 0 , 1 , 0 , 0 , 33 , 260
1987-05-20 , 128 , 0 , 6 , 17 , 2 , 1 , 0 , 0 , 0 , 0 , 0 , 26 , 154
1987-05-21 , 72 , 0 , 11 , 9 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 22 , 145
1987-05-22 , 57 , 0 , 6 , 8 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 15 , 99
1987-05-23 , 49 , 0 , 4 , 2 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 7 , 58
1987-05-24 , 57 , 14 , 2 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 17 , 74
1987-05-25 , 18 , 0 , 3 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 3 , 40
1987-05-26 , 20 , 0 , 3 , 4 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 7 , 50
1987-05-27 , 16 , 0 , 3 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 3 , 59
1987-05-28 , 15 , 0 , 0 , 2 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 3 , 40
1987-05-29 , 10 , 0 , 1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 9
1987-05-30 , 13 , 0 , 0 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 14
1987-05-31 , 8 , 0 , 3 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 4 , 13
1987-06-01 , 2 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 22
1987-06-02 , 23 , 0 , 6 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 6 , 24
1987-06-03 , 20 , 0 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 33
1987-06-04 , 10 , 0 , 4 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 5 , 19
1987-06-05 , 10 , 3 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 4 , 12
1987-06-06 , 5 , 0 , 2 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 3 , 7
1987-06-07 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 4
1987-06-08 , 2 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0
1987-06-09 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0
1987-06-10 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 59")
demo.Fish <- read.csv(demo.data, header=TRUE)


demo.Fish[1:5,]

# Now to extract the subset of data, do any fancy adjustments, and fit the data.

demo.prefix <- "demo-CR-1987-AS-TSPNDE"
demo.title  <- "Conne River 1987 Simulated Atlantic Salmon "
cat("*** Starting ",demo.title, "\n\n")

# We do some fixups because the data entry forces you to enter extra rows of zeros
# at the end because the u2 vector can be longer than the n1 vector.
demo.n1 <- demo.Fish$tagged
demo.n1 <- demo.n1[ 1:(length(demo.n1)-2)]  # last two entries are ficticious

demo.m2 <- demo.Fish[, paste("X",0:9,sep="")]
demo.m2 <- as.matrix(demo.m2)
demo.m2 <- demo.m2[ 1:length(demo.n1),]     # last two rows are ficticious

demo.u2 <- demo.Fish$Untagged

# what fractio of the day was sampled?
demo.sampfrac <- rep(1,length(demo.u2)) # values are on a daily basis

# what is the strata identification number?
demo.day <- 1:length(demo.u2)

# are there any jumps in the abundance?
demo.jump.after <- NULL    # list sample times after which jump in number occurs

# are there any bad values that need to be adjusted?
demo.bad.n1     <- c()     # list sample times of bad n1 values
demo.bad.m2     <- c()     # list sample times of bad m2 values
demo.bad.u2     <- c()     # list sample times of bad u2 values


demo.cr.1987.as.tspnde <- TimeStratPetersenNonDiagError_fit(
                  title=      demo.title,
                  prefix=     demo.prefix,
                  time=       demo.day,
                  n1=         demo.n1, 
                  m2=         demo.m2, 
                  u2=         demo.u2,
                  sampfrac=   demo.sampfrac,
                  jump.after= demo.jump.after,
                  bad.n1=     demo.bad.n1,
                  bad.m2=     demo.bad.m2,
                  bad.u2=     demo.bad.u2,
                  #engine="openbugs", # how to specify openbugs engine
                  debug=TRUE
                  )

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
 

save(list=c("demo.cr.1987.as.tspnde"), file="demo-cr-1987-as-tspnde-saved.Rdata")  # save the results from this run

cat("\n\n\n ***** FILES and GRAPHS saved in \n    ", getwd(), "\n\n\n")
print(dir())

# move up the directory
setwd("..")

cat("\n\n\n ***** End of Demonstration *****\n\n\n")

