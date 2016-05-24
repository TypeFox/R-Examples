opar <- par(ask = dev.interactive(orNone = TRUE))


cat("******************************************************\n",
    "We give here an example of management of an object of \n",
    "ltraj",
    "\n******************************************************\n")


###############################################
###
### Example 1: preparing the "ltraj" from raw data


data(puechabon)

locs <- puechabon$locs
locs[1:4,]
xy <- locs[,c("X","Y")]

### Conversion of the date to the format POSIX
da <- as.character(locs$Date)
da <- as.POSIXct(strptime(as.character(locs$Date),"%y%m%d"))

## object of class "ltraj"
(ltr <- as.ltraj(xy, da, id = locs$Name))


## look at the data:
ltr

## the data are not regular: see the distribution of dt (in hours)
## according to the date
is.regular(ltr)

plotltr(ltr, "dt/3600/24")

## First, note that Chou is monitored two successive summers
## Now, cut the trajectory of Chou into two summers

foo <- function(dt) {
    return(dt> (100*3600*24))
}

## The function foo returns TRUE if dt is longer than 100 days
## We use it to cut ltr:

l2 <- cutltraj(ltr, "foo(dt)", nextr = TRUE)
l2
plot(l2, perani = FALSE)

## Now, look again at the time lag:

plotltr(l2, "dt/3600/24")


## The relocations have been collected daily, and there are many
## missing values
## We set the missing values in this trajectory. We first
## set the reference date: the hour should be exact (i.e. minutes=0):
refda <- strptime("00:00", "%H:%M")
refda

## Set the missing values
l3 <- setNA(l2, refda, 1, units = "day")

## look at the data:
l3


## The trajectory is now regular, but there is a lot of
## missing values!!
summaryNAltraj(l3)


## Are the missing values randomly distributed in the trajectory?
runsNAltraj(l3)

## yes, it seems so...
## Trajectory is ready for the analysis

cat("*******************************************************\n",
    "The deeply commented source for this demo can be found in the file:\n",
    file.path(system.file(package = "adehabitat"), "demo", "managltraj.r\n"),
    "Examples of analysis are given in demo(analysisltraj)\n",
    "******************************************************\n")
