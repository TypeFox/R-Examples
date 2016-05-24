# 2009-12-06 CJS First edition
# 2014-09-01 CJS remove prompts; jags; engine

# This is a demonstration of how to call the Time Stratified Petersen with Diagonal Entries (TSPDE) program
# SEPARATING WILD VS HATCHERY STEELHEAD SALMON.
#
# This differs from the Wild vs Hatchery Chinook salmon in that all hatchery raised steelhead are marked,
# so there is complete separation by age and (wild/hatchery).
# There are 3 population of interest, Wild.YoY, Hatchery.Age1+, and  Wild.Age1+.

# It is based on the analysis of California Junction City 2003 Steelhead data and is the example used
# in the Trinity River Project.
#
# In each julian week j, n1[j] are marked and released above the rotary screw trap.
# Of these, m2[j] are recaptured. All recaptures take place in the week of release, i.e. the matrix of
# releases and recoveries is diagonal.
# The n1[j] and m2[j] establish the capture efficiency of the trap.
#
# At the same time, u2[j] unmarked fish are captured at the screw trap.
# These fish are a mixture of wild and hatchery raised steelhead salmon. 
# The u2[j] are separated into u2.W.YoY[j] (wild, YoY steelhead),
#                              u2.W.1  [j] (wild, age 1+ steelhead), and
#                              u2.H.1  [j] (hatchery, age 1+ steelhead)
#
# The program assumes that the trap was operating all days of the week. The sampfrac[j] variable
# gives the proportion of days the trap was operating. For example, if the trap was operating for 3 of the
# 7 days in a week, then sampfrac[j]<- 3/7
#
#  The vectors bad.m2 etc indicates which julian weeks something went wrong. 
# 
#  The prefix is used to identify the output files for this run.
#  The title  is used to title the output.

#  Notes:
#    Marking and releasing of steelhead took place in only a few weeks! The hierarchical model
#    will extrapolate outside these weeks to estimate the capture rate.

# Create a directory to store the results Test and then create the
# directory
if(file.access("demo-demo-TSPDE-WHsteel")!=0){
  dir.create("demo-demo-TSPDE-WHsteel", showWarnings=TRUE)
}
setwd("demo-demo-TSPDE-WHsteel")

## Load BTSPAS library
library(BTSPAS)

# Get the data. In many cases, this is stored in a *.csv file and read into the program
# using a read.csv() call. In this demo, the raw data is assigned directly as a vector.
#

# Indicator for the week.
demo.jweek <- c(9,   10,   11,   12,   13,   14,   15,   16,   17,   18,
          19,   20,   21,   22,   23,   24,   25,   26,   27,   28, 
          29,   30,   31,   32,   33,   34,   35,   36,   37,   38,
          39,   40,   41,   42,   43,   44,   45,   46)

# Number of marked fish released in each week.
demo.n1 <- c(   0,    0,    0,  999, 1707, 1947, 2109,  972,  687,    0,
           0,    0,    0,    0,    0,    0,    0,    3,    0,    0,
           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
           0,    0,    0,    0,    0,    0,    0,    0)

# Number of marked fish recaptured in the week of release. No marked fish
# are recaptured outside the week of release.
demo.m2 <- c(   0,    0,    0,    5,   13,   39,    7,    1,    0,    0,
           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
           0,    0,    0,    0,    0,    0,    0,    0)

# Number of unmarked fish captured at the trap in each week separated by wild/hatchery and age
# Wild YoY
demo.u2.W.YoY <-c(0,  0,    0,    0,   11,    0,    0,    0,    1,   33,
            31, 11,   78,   46,   35,   30,  309,  278,  207,  196,
           613,764,  556,  250,  106,  413,  995,  357,  181,   53,
            29,  3,    5,   14,    8,   19,   46,  229)

# Wild age 1+
demo.u2.W.1 <- c(58,357,  720,  850,  585,  532,  873,  303,  291,   12,
           101, 47,   49,   44,   50,   38,   58,   36,   13,    5,
            12, 15,   11,   12,   13,   12,   28,   10,    8,    3,
             2,  0,    0,    4,   10,    7,    4,    7)

# Hatchery age 1+
demo.u2.H.1 <- c( 0,  2,    0, 4643, 5758, 4220, 2328, 1474,  875,   39,
            15, 13,   26,   22,   59,   15,    8,    4,    2,    0,
             0,  0,    0,    0,    0,    0,    1,    0,   27,    2,
             0,  0,    0,    0,    0,    0,    0,    0)

# What fraction of the week was sampled?
demo.sampfrac<-c(3,   8,    6,    7,    7,    7,    7,    7,    7,    7,
            7,   7,    7,    7,    7,    7,    7,    7,    7,    7,
            6,   7,    7,    7,    7,    7,    7,    7,    7,    7,
            7,   7,    7,    7,    7,    7,    7,    5)/7

# After which weeks do the hatchery fish start to arrive. Prior to this point, all fish are wild and it is not
# necessary to separate out the wild vs hatchery
demo.hatch.after <- c(11)  # julian weeks after which hatchery fish arrive.

# Which julian weeks have "bad"  values. These will be set to missing and estimated.
demo.bad.m2       <- c()     # list of julian weeks with bad m2 values
demo.bad.u2.W.YoY <- c()     # list of julian weeks with bad u2.W.YoY values
demo.bad.u2.W.1   <- c()     # list of julian weeks with bad u2.W.1   values
demo.bad.u2.H.1   <- c()     # list of julian weeks with bad u2.H.1   values



# The prefix for the output files:
demo.prefix <- "demo-JC-2003-ST-TSPDE-WH" 

# Title for the analysis
demo.title <- "Junction City 2003 Steelhead - Separation of Wild and Hatchery YoY and Age 1+ Steelhead"



cat("*** Starting ",demo.title, "\n\n")

# Make the call to fit the model and generate the output files
demo.jc.2003.st.tspde <- TimeStratPetersenDiagErrorWHSteel_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek,
                  n1=demo.n1,
                  m2=demo.m2,
                  u2.W.YoY=demo.u2.W.YoY, u2.W.1=demo.u2.W.1, u2.H.1=demo.u2.H.1,
                  sampfrac=demo.sampfrac,
                  hatch.after=demo.hatch.after,
                  bad.m2=demo.bad.m2,
                  bad.u2.W.YoY=demo.bad.u2.W.YoY,
                  bad.u2.W.1  =demo.bad.u2.W.1,
                  bad.u2.H.1  =demo.bad.u2.H.1,
		  #engine="openbugs", # show how to call openbugs
                  debug=TRUE  # this generates only 10,000 iterations of the MCMC chain for checking.
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
 
# save the results in a data dump that can be read in later using the load() command.
# Contact Carl Schwarz (cschwarz@stat.sfu.ca) for details.
save(list=c("demo.jc.2003.st.tspde"), file="demo-jc-2003-st-tspde-saved.Rdata")  # save the results from this run


cat("\n\n\n ***** FILES and GRAPHS saved in \n    ", getwd(), "\n\n\n")
print(dir())

# move up the directory
setwd("..")


cat("\n\n\n ***** End of Demonstration *****\n\n\n")

