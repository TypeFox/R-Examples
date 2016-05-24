#  Set up the data objects used in the examples

examples <- "c:\\Program Files\\Insightful\\Splus70\\fdaS\\examples\\"

#  ------------  gait  --------------------------

hip  <- matrix(scan(paste(examples,"gait\\hip.dat",sep=""),  0), 20, 39)
knee <- matrix(scan(paste(examples,"gait\\knee.dat",sep=""), 0), 20, 39)

#  set up the argument values

gaittime  <- (1:20)/21
gaitrange <- c(0,1)

#  set up a three-dimensional array of function values

gaitarray <- array(0, c(20, 39, 2))
dimnames(gaitarray) <- list(NULL, NULL, c("Hip Angle", "Knee Angle"))
gaitarray[,,1] <- hip
gaitarray[,,2] <- knee

#  ------------  goods index  --------------------------

temp    <- matrix(scan(paste(examples,"goodsindex\\nondurprod.dat",sep=""),0), 18, 81)
tempmat <- temp[2:13,]
tempmat[12,81] <- 0
nondurables <- matrix(tempmat, 12*81, 1)
nondurables <- nondurables[1:971]
ndur <- 971

#  for completeness, make dec 99 equal to dec 98, jan 00 equal to jan 99

nondurables <- c(nondurables,nondurables[961])
nondurables <- c(nondurables,nondurables[962])
ndur <- 973

#  set up time values

durtime <- (0:(ndur-1))/12 + 1919
goodsrange <- c(1919,2000)

#  compute log nondurables

lognondur <- log10(nondurables)

#  ------------  Berkeley growth  --------------------------

ncasem <- 39
ncasef <- 54
nage   <- 31

hgtm <- t(matrix(scan(paste(examples,"growth\\hgtm.dat",sep="") ,0), ncasem, nage, byrow=T))
hgtf <- t(matrix(scan(paste(examples,"growth\\hgtf.dat",sep="") ,0), ncasem, nage, byrow=T))

age <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))

#  ------------  handwriting  --------------------------

temp <- array(scan(paste(examples,"handwrit\\fdareg.dat",sep=""),0), c(20,2,1401))

#  set up a three-dimensional array

fdaarray <- array(0, c(1401, 20, 2))
fdaarray[,,1] <- t(temp[,1,])/1000
fdaarray[,,2] <- t(temp[,2,])/1000
dimnames(fdaarray) <- list(NULL, NULL, c("X", "Y") )

#  Set up time values and range.
#  It is best to choose milliseconds as a time scale
#  in order to make the ratio of the time
#  unit to the inter-knot interval not too
#  far from one.  Otherwise, smoothing parameter values
#  may be extremely small or extremely large.

fdatime  <- seq(0, 2300, len=1401)

#  ------------  lip  --------------------------

lipmat <- matrix(scan(paste(examples,"lip\\lip.dat",sep=""), 0), 51, 20)

liptime  <- seq(0,1,.02)

#  ------------  melanoma  --------------------------

tempmat <- t(matrix(scan(paste(examples,"melanoma\\melanoma.dat",sep=""), 0), 3, 37))

year  <- tempmat[,2]
mela  <- tempmat[,3]
nyear <- length(year)

#  ------------  pinch  --------------------------

pinchmat   <- matrix(scan(paste(examples,"pinch\\pinch.dat",sep=""),0), 151, 20, byrow=T)

pinchtime  <- seq(0,150,len=151)/600

#  ------------  refinery  --------------------------

refinery <- t(matrix(scan(paste(examples,"refinery\\refinery.dat",sep=""), 0), 3, 193))

tval <- refinery[,1]  #  observation time
uval <- refinery[,2]  #  reflux flow
yval <- refinery[,3]  #  tray 47 level

#  center the data on mean values prior to change

uval <- uval - mean(uval[1:60])
yval <- yval - mean(yval[1:60])

#  ------------  daily weather  --------------------------

tempav <- matrix(scan(paste(examples,"weather\\dailtemp.dat",sep=""),0), 365, 35)
precav <- matrix(scan(paste(examples,"weather\\dailprec.dat",sep=""),0), 365, 35)

#  set up the times of observation at noon

daytime   <- (1:365)-0.5
dayrange  <- c(0,365)
dayperiod <- 365

#  day values roughly in weeks

weeks <- seq(0,365,length=53)   

#  define 11-character names for stations

place <- c(
"Arvida     ", "Bagottville", "Calgary    ", "Charlottvl ", "Churchill  ", "Dawson     ",
"Edmonton   ", "Fredericton", "Halifax    ", "Inuvik     ", "Iqaluit    ", "Kamloops   ",
"London     ", "Montreal   ", "Ottawa     ", "Pr. Albert ", "Pr. George ", "Pr. Rupert ",
"Quebec     ", "Regina     ", "Resolute   ", "Scheffervll", "Sherbrooke ", "St. Johns  ",
"Sydney     ", "The Pas    ", "Thunderbay ", "Toronto    ", "Uranium Cty", "Vancouver  ",
"Victoria   ", "Whitehorse ", "Winnipeg   ", "Yarmouth   ", "Yellowknife")

dimnames(tempav) <- list(NULL,place)
dimnames(precav) <- list(NULL,place)

#  set up indices that order the stations from east to west to north

geogindex <- c(24,  9, 25, 34,  4,  8, 22,  1,  2, 19, 23, 14, 15, 28, 13, 
               27, 33, 26,  5, 20, 16, 29,  7,  3, 12, 30, 31, 17, 18, 32, 
                6, 35, 11, 10, 21)

#  put the stations in geographical order, from east to west to north
#  rather in the original alphatical order.

tempav <- tempav[,geogindex]
precav <- precav[,geogindex]
place  <- place[geogindex]



