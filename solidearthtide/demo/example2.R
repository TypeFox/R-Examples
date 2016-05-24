library(Rgnuplot)
# Example 2
# calculate a global snapshot of solid
# earth tide displacements at an array of grid points, for a single snapshot in time.
# original MATLAB code:
# Chad Greene, 2015
# Solid Earth Tide Predictions
# http://www.mathworks.com/matlabcentral/fileexchange/50327-sold-earth-tide-predictions/content/earthtide

# define latitude and longitude
lonIn<-matrix(rep(seq(-170,170,20),18),18,byrow=TRUE)
latIn<-matrix(rep(seq(85,-85,-10),18),18,byrow=FALSE)
# define start and end dates
timeIn <- Datenum(as.Date("may 16, 2007", "%b %d, %Y"))+(9*3600+23*60)/86400 # may 16, 2007 9:23	733178.390972222
upOut <- matrix(0,18,18)
northOut <- matrix(0,18,18)
eastOut <- matrix(0,18,18)
# loop to fill the data matrices
for (n in 1:18^2){
tmp <- Earthtide(latIn[n], lonIn[n], timeIn,TRUE)
upOut[n] <- tmp$upOut
northOut[n] <- tmp$northOut
eastOut[n] <- tmp$eastOut
}
# use Rgnuplot for the graphics
h1<-Gpinit()
Gpsetwd(h1)
# create data file with all the output
a<-cbind(rev(c(lonIn)),rev(c(latIn)),c(upOut),c(eastOut),c(northOut))
write.table(a, file = "allOut.tab", sep = "\t",row.names =FALSE, col.names =FALSE, quote =FALSE)
# plot heatmap + arrows
Gprun ( 'set view map
set cblabel "Vertical displacement (cm)"
set xlabel "Longitude (deg)"
set ylabel "Latitude (deg)"
splot "allOut.tab" using 1:2:($3*100) with image notitle, "allOut.tab"  u 1:2:(0):(-$4*150):(-$5*150):(0)  with vectors arrowstyle 1 notit' , TRUE )

