library(RFOC)



print("Demonstration of Stereonet Usage")

print("Draw an Equal-Angle Stereonet (Schmidt Net)")

print("Create a stereo net structure and save, then plot")

MN = makenet()

pnet(MN)

readline("To Continue Hit Enter Key\n")
pnet(MN)
mtext(side=3,"Left Click 3 times on net to record points")
print("Left Click 3 times on net to record points")

rpoints = AXpoint(UP=TRUE, col=2, n=3)

mtext(side=1, "Left Click 3 times on net to record more points")


print("Left Click 3 times on net to record more points")

bpoints = AXpoint(UP=TRUE, col=4, n=3)

readline("To Continue Hit Enter Key\n")
mtext(side=3,"POints are saved, and replotted:")
pnet(MN)

qpoint(rpoints$az, rpoints$dip, UP=TRUE,  col=2)
qpoint(bpoints$az, bpoints$dip, UP=TRUE, col=4)


readline("To Continue Hit Enter Key\n")

pnet(MN)
qpoint(rpoints$az[1], rpoints$dip[1], UP=TRUE,  col=2)
dip = rpoints$dip[1]
rplane = lowplane(rpoints$az[1]-90, rpoints$dip[1], UP=TRUE, col='red')
##########  must subtract 90 from the azimuth to get correct orientation
mtext(side=3,"First point is chosen and plotted as a pole:")

readline("To Continue Hit Enter Key\n")

pnet(MN)

addsmallcirc(34, 52, 10, lty=2)
mtext(side=3,"Click three times in small black circle")

rpoints = AXpoint(UP=TRUE, col=2, n=3)

###########  plot great circles associated with poles
for(i in 1:length(rpoints$az)) {
  rplane = lowplane(rpoints$az[i]-90, rpoints$dip[i], UP=TRUE, col='red')
}

########## 
readline("To Continue Hit Enter Key\n")
########## ##########   create a bunch of poles and calcualte statistics

az = rnorm(25, m=-134, sd=20)
dip = rnorm(25, m=52, sd=10)
pnet(MN)
qpoint(az, dip, UP=FALSE,  col=2)

A95 = alpha95(az, dip)

qpoint(A95$Dr, 90-A95$Ir, UP=FALSE,  col='black', pch=6)

addsmallcirc(A95$Dr, 90-A95$Ir, A95$Alph95, lty=1, lwd=2)
mtext(side=3,"Random Points and 95% Confidence of Mean Vector")

readline("To Continue Hit Enter Key\n")
########## ##########

library(MASS)

#########  create structure of points:
   PZZ = focpoint(az, 90-dip, col = "red", pch = 3, lab = "",
        UP = FALSE, PLOT = FALSE)

KP = kde2d(PZZ$x, PZZ$y, n = 100, lims = c(-1, 1, -1, 1))

##########  use meshgrid to set up image:
######                blank out points outside of sphere
  M = meshgrid(KP$x, KP$y)
        flag = sqrt(M$x^2 + M$y^2) > 1
        KP$z[flag] = NA

pnet(MN)

###########  location of center of focal sphere:
fx = 0
fy = 0

##   size of focal sphere:
siz = 1

##########  do image:
        image(fx + siz * KP$x, fy + siz * KP$y, KP$z, add = TRUE,
            col = terrain.colors(100))
pnet(MN, add=TRUE)

mtext(side=3,"Smoothed Random Points")


readline("To Continue Hit Enter Key\n")

##########  add contours

cont.col ='blue'
 cline.list = contourLines(fx + siz * KP$x, fy + siz * KP$y,
            KP$z)
        templines <- function(clines) {
            cx = clines[[2]]
            cy = clines[[3]]
            flag = sqrt((cx - fx)^2 + (cy - fy)^2) > siz
            cx[flag] = NA
            cy[flag] = NA
            lines(cx, cy, col = cont.col)
        }
        invisible(lapply(cline.list, templines))

points(PZZ$x, PZZ$y)


points(PZZ$x, PZZ$y)


