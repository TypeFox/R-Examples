#
# 2010-03-30 CJS first edition
# 2014-09-01 CJS remove prompts; jags; engine

# This is a demonstration of how to call the Time Stratified Petersen with Diagonal Entries (TSPDE) program
# SEPARATING WILD VS HATCHERY CHINOOK SALMON when Age1 chinook are present (residualized) from last year.

# It is based on the analysis of California North Fork Chinook data.
#
# In each julian week j, n1[j] are marked and released above the rotary screw trap.
# Of these, m2[j] are recaptured. All recaptures take place in the week of release, i.e. the matrix of
# releases and recoveries is diagonal.
# The n1[j] and m2[j] establish the capture efficiency of the trap.
#
# At the same time, u2[j] unmarked fish are captured at the screw trap.
# These fish are a mixture of YoY and Age1 wild and hatchery raised chinook salmon. 
# A portion (clip.rate.H.YoY, clip.rate.H.1) of the YoY and Age1 hatchery raised fish 
#    are adipose fin clipped and can be recognized as hatchery raised.
# The unclipped fish are a mixture of wild and hatchery fish which must be separated.
# Hence the u2[j] are separated into u2.A.YoY[j] (YoY adipose clipped fish known to be hatchery), and
#                                    u2.N.YoY[j] (YoY unclipped fish) which are mixture of hatchery and wild fish
#                                    u2.A.1  [j] (Age1 adipose clipped fish known to be hatchery), and
#                                    u2.N.1  [j] (Age1 unclipped fish) which are mixture of hatchery and wild fish
#
# The program assumes that the trap was operating all days of the week. The sampfrac[j] variable
# gives the proportion of days the trap was operating. For example, if the trap was operating for 3 of the
# 7 days in a week, then sampfrac[j]<- 3/7
#
#  The program tries to fit a single spline to the entire dataset. However, in julian week 23
#  hatchery released YoY fish started to arrive at the trap resulting in sudden jump
#  in abundance. The hatch.after vector gives the julian weeks just BEFORE the sudden jump in YoY,
#  i.e. the spline is allowed to jump AFTER the julian weeks in jump.after.
#
#  The Age1 fish are fit with a single spline and NO jumps are allowed
#
#  The vector bad.m2 indicates which julian weeks something went wrong. For example, the
#  number of recoveries in julian week 41 is far below expectations and leads to impossible
#  Petersen estimate for julian week 41.
# 
#  The prefix is used to identify the output files for this run.
#  The title  is used to title the output.

# Create a directory to store the results Test and then create the
# directory
if(file.access("demo-TSPDE-WHchinook2")!=0){
  dir.create("demo-TSPDE-WHchinook2", showWarnings=TRUE)
}
setwd("demo-TSPDE-WHchinook2")

## Load BTSPAS library
library(BTSPAS)

# Get the data. In many cases, this is stored in a *.csv file and read into the program
# using a read.csv() call. In this demo, the csv data is stored as a data stream rather than a 
# separate file
#

# Get the raw data and read it in
demo.data <- textConnection(
"Site , Trap_id , Year, Species, Spp, jweek, TotalCatch, TotalMarks,TotalRecaps,Week0_,Week_1_,Week_2_,Week_3_,Week_4_,Week_5_,DaysOperating,CH_YOY_NC,CH_YOY_AD,CH_1_NC,CH_1_AD
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,2,25,0,0,0,0,0,0,0,0,2,15,0,10,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,3,185,0,0,0,0,0,0,0,0,7,94,0,90,1
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,4,499,833,52,52,0,0,0,0,0,7,385,0,112,2
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,5,1314,852,67,66,0,0,1,0,0,7,1162,0,140,12
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,6,705,1495,77,76,0,1,0,0,0,7,592,0,103,10
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,7,1249,1356,182,180,1,1,0,0,0,5,1151,0,94,4
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,8,2386,1889,145,142,0,3,0,0,0,6,2258,0,121,7
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,9,1205,2934,89,62,26,0,0,1,0,7,1123,0,80,2
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,10,2339,1546,53,53,0,0,0,0,0,7,2277,0,57,5
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,11,2523,4001,232,231,1,0,0,0,0,7,2492,0,27,4
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,12,1681,2955,158,158,0,0,0,0,0,7,1579,0,88,14
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,13,1096,529,14,12,2,0,0,0,0,7,1046,0,45,5
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,14,782,1172,49,49,0,0,0,0,0,7,766,0,13,3
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,15,2712,3204,232,121,111,0,0,0,0,7,2702,0,9,1
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,16,10428,1328,57,57,0,0,0,0,0,7,10408,0,18,2
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,17,12163,3540,114,114,0,0,0,0,0,7,12145,0,15,3
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,18,187,4791,45,40,4,1,0,0,0,7,186,0,1,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,19,409,4808,11,10,1,0,0,0,0,7,407,0,2,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,20,862,5952,44,15,27,2,0,0,0,7,862,0,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,21,465,3852,55,52,2,1,0,0,0,5,465,0,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,22,751,2621,17,15,2,0,0,0,0,7,724,0,27,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,23,5714,2131,37,35,2,0,0,0,0,7,4860,854,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,24,4334,5002,152,147,5,0,0,0,0,7,3539,794,1,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,25,5501,3706,120,117,3,0,0,0,0,7,4597,904,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,26,4527,1225,44,43,1,0,0,0,0,7,3819,708,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,27,4062,723,45,42,3,0,0,0,0,5,3300,762,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,28,6816,2895,167,166,1,0,0,0,0,7,5460,1356,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,29,3532,1395,117,117,0,0,0,0,0,7,2918,614,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,30,2672,479,77,77,0,0,0,0,0,7,2252,420,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,31,1529,964,74,0,74,0,0,0,0,7,1240,289,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,32,515,2803,288,288,0,0,0,0,0,7,428,87,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,33,578,952,51,49,2,0,0,0,0,4,464,114,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,34,568,880,126,126,0,0,0,0,0,4,515,53,0,0
Pear Tree,TRN-5; TRN-8,2009,Chinook,1,35,114,0,0,0,0,0,0,0,0,2,93,21,0,0 ")

demo.Fish <- read.csv(demo.data, header=TRUE)


# Extract the relevant data
demo.jweek <- demo.Fish$jweek

# Number of marked fish and recaptures of same (number recovered outside of week of release is negligibl)
demo.n1  <- demo.Fish$TotalMarks
demo.m2  <- demo.Fish$TotalRecaps

# Number of unmarked fish captured at the trap in each week separated by adipose and wild clips and YoY and Age1.
demo.u2.A.YoY <- demo.Fish$CH_YOY_AD
demo.u2.N.YoY <- demo.Fish$CH_YOY_NC

demo.u2.A.1   <- demo.Fish$CH_1_AD
demo.u2.N.1   <- demo.Fish$CH_1_NC

# What fraction of the week was sampled?
demo.sampfrac <- demo.Fish$DaysOperating/7

# After which weeks do the YoY hatchery fish start to arrive. 
# Prior to this point, all YoY fish are wild and it is not
# necessary to separate out the YoY wild vs hatchery
demo.hatch.after.YoY <- c(22)  # julian weeks after which YoY hatchery fish arrive.

# Which julian weeks have "bad"  values. These will be set to missing and estimated.
demo.bad.m2       <- c()     # list of julian weeks with bad m2 values
demo.bad.u2.A.YoY <- c()     # list of julian weeks with bad u2.A.YoY values
demo.bad.u2.N.YoY <- c()     # list of julian weeks with bad u2.N.YoY values
demo.bad.u2.A.1   <- c()     # list of julian weeks with bad u2.A.YoY values
demo.bad.u2.N.1   <- c()     # list of julian weeks with bad u2.N.YoY values

# The clipping fraction for the current YoY and last year's YoY (which are now Age 1 fish)
demo.clip.frac.H.YoY <- .25    # what fraction of the YoY  hatchery fish are adipose fin clipped
demo.clip.frac.H.1   <- .25    # what fraction of the Age1 hatchery fish are adipose fin clipped

# The prefix for the output files:
demo.prefix <- "demo-NF-2009-CH-TSPDE-WHC2" 

# Title for the analysis
demo.title <- "North Fork 2009 Chinook - Separation of YoY and Age 1 Wild and Hatchery Chinook"


cat("*** Starting ",demo.title, "\n\n")

# Make the call to fit the model and generate the output files
demo.nf.2009.ch2.tspde <- TimeStratPetersenDiagErrorWHChinook2_fit(
                  title=demo.title,
                  prefix=demo.prefix,
                  time=demo.jweek,
                  n1=demo.n1,          
                  m2=demo.m2,         
                  u2.A.YoY=demo.u2.A.YoY,
                  u2.N.YoY=demo.u2.N.YoY,
                  u2.A.1  =demo.u2.A.1,
                  u2.N.1  =demo.u2.N.1,
                  sampfrac=demo.sampfrac,
                  clip.frac.H.YoY= demo.clip.frac.H.YoY,
                  clip.frac.H.1  = demo.clip.frac.H.1,
                  hatch.after.YoY=demo.hatch.after.YoY,
                  bad.m2=demo.bad.m2, 
                  bad.u2.A.YoY=demo.bad.u2.A.YoY, bad.u2.N.YoY=demo.bad.u2.N.YoY,
                  bad.u2.A.1  =demo.bad.u2.A.1  , bad.u2.N.1  =demo.bad.u2.N.1,
		  #engine="openbugs",   # show how to call openbugs
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
save(list=c("demo.nf.2009.ch2.tspde"), file="demo-nf-2009-ch-tspde-whc2-saved.Rdata")  # save the results from this run


cat("\n\n\n ***** FILES and GRAPHS saved in \n    ", getwd(), "\n\n\n")
print(dir())

# move up the directory
setwd("..")


cat("\n\n\n ***** End of Demonstration *****\n\n\n")


