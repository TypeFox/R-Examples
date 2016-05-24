
data(USAmap)


####  set up limits for USA map
USALL=list()
USALL$lat=c(24.72853,49.62741)
USALL$lon=c(229.29389,296.41803)

## set UTM projection
PROJ = setPROJ(type = 2, LAT0 =mean(USALL$lat), LON0 = mean(USALL$lon) )


####

####  plot using no projection:
#### plotGEOmap(USAmap, LIM= c(USALL$lon[1], USALL$lat[1], USALL$lon[2], USALL$lat[2]    )  ,  add=FALSE, shiftlon=0)

#### readline("Hit Return>")

####   get(getOption("device"))(width=13, height=10)  
#### dev.new()

####  plot with UTM  projection:
plotGEOmapXY(USAmap, LIM= c(USALL$lon[1], USALL$lat[1], USALL$lon[2], USALL$lat[2]    )  , PROJ=PROJ, add=FALSE, shiftlon=0)


##########  after clicking on the screen for the jet stream
js1=list()
js1$x=c(-2264.91,-1959.38,-1719.32,-1102.81, -808.19, -480.83, -104.38,  272.08,  757.66, 1030.45, 1275.97, 1772.46, 2285.31, 2634.49)
js1$y=c(1803.45,1594.03,1357.04,1125.57, 883.07, 602.00, 249.27, 199.67, 276.83, 453.19, 684.66, 965.74,1197.21,1312.95)


jetstream  =getsplineG(js1$x, js1$y, kdiv=20)


###  after clicking on the screen for the cold front
c1=list()
c1$x=c(-508.113,  26.566, 326.641, 490.318)
c1$y=c( 673.64,1037.39,1384.60,1753.85)

coldfront =   getsplineG(c1$x, c1$y, kdiv=20)

###  after clicking on the screen for the warm front1

w1=list()
w1$x=c( 499.668, 614.244, 824.299, 996.162,1101.190,1153.704,1210.992,1249.184)
w1$y=c(1758.352,1657.081,1406.314,1150.725, 986.761, 885.490, 750.462, 678.125)

warmfront =   getsplineG(w1$x, w1$y, kdiv=20)


#####  make plot of usa, projected with UTM parameters provided above:
plotusa()
par(xpd=TRUE)

###  add jet stream:
lines(jetstream, col='green', lwd=5)
g = PointsAlong(jetstream$x, jetstream$y, N=8)        
  pin = par("pin")
  u = par("usr")
   umm =   (u[4]-u[3])/pin[2]/25.4
  lenarr =  umm*5
rot=list(cs = g$rot$sn , sn = -g$rot$cs )
teeth(g$x, g$y, lenarr, rot, col='green', border='green')

###  add cold front:
lines(coldfront, col='blue', lwd=2)

##########  put in teeth
thrust(coldfront$x, coldfront$y,  h=3, N=10, REV=TRUE, col='blue', lty=2)

### add warm front:
lines(warmfront, col='red', lwd=2)

###  get points along the warm front:
g = PointsAlong(warmfront$x, warmfront$y, N=6)


##  put in half moons:
horseshoe(g$x  , g$y , r1=3, r2=3*.8, h2=0, h1=0, rot=g$rot , col='red', fill=TRUE)

##    zeb=locator(1)
###   set location for horizontal kilometer scale
zeb=list()
zeb$x=c(197.727896066)
zeb$y=c(-1155.81158234)

zebra(zeb$x[1],zeb$y[1], 1000, 100, 60, lab="Km", cex=.6)


LOW=list()
LOW$x=c(554.03313,618.79783,595.66758,387.49533,211.70544,234.83568,373.61718, 554.03313)
LOW$y=c(1841.8045,1916.5726,2024.0517,2052.0898,1967.9757,1832.4585,1771.7095, 1841.8045)


lowp =   getsplineG(LOW$x, LOW$y, kdiv=20)


sk = 1.3
z = PointsAlong(lowp$x, lowp$y, N=20, endtol = 0)

polygon(lowp, border='blue', col=NULL )

 perpen(z$x,z$y, h=sk, rot=z$rot, lwd=1, col='blue' )

text(mean(lowp$x), mean(lowp$y), labels="Low")

          

#####  tags = locator(3)


tags=list()
tags$x=c(1025.890219, -38.101260,1798.440553)
tags$y=c(1304.40892,1117.48871, 883.83844)

angs=list()
angs$x=c( 928.743171,1118.411217,-149.126457,  68.297888,1664.285106,1937.222050)
angs$y=c(1435.25307,1173.56477,1024.02860,1206.27581, 837.10839, 953.93352)


rots =  180*atan2( -(angs$y[c(1,3,5)]-angs$y[c(2,4,6)]), -(angs$x[c(1,3,5)]-angs$x[c(2,4,6)]) )/pi

############  set up angles for roation based on screen locator()
######  angs = locator(6) ; segments(angs$x[c(1,3,5)],angs$y[c(1,3,5)] , angs$x[c(2,4,6)],angs$y[c(2,4,6)])

LABS = c("Warm Front", "Cold Front", "Jet Stream")
COLS = c('red', 'blue', 'green')
########    
###  i = 1;  plotusa()


for(i in 1:length(LABS))
text(tags$x[i], tags$y[i], labels=LABS[i], font=2, cex=1.4, col=COLS[i], srt=rots[i])


