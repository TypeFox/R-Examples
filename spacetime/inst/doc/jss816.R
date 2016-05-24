### R code from vignette source 'jss816.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: jss816.Rnw:86-88
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
set.seed(1331)


###################################################
### code chunk number 2: jss816.Rnw:185-187
###################################################
library(foreign)
read.dbf(system.file("shapes/sids.dbf", package="maptools"))[1:5,c(5,9:14)]


###################################################
### code chunk number 3: jss816.Rnw:198-200
###################################################
data("wind", package = "gstat")
wind[1:6,1:12]


###################################################
### code chunk number 4: jss816.Rnw:214-216
###################################################
data("Produc", package = "plm")
Produc[1:5,1:9]


###################################################
### code chunk number 5: jss816.Rnw:248-321
###################################################
opar = par()
par(mfrow=c(2,2))
# 1:
s = 1:3
t = c(1, 1.75, 3, 4.5)
g = data.frame(rep(t, each=3), rep(s,4))
col = 'blue'
pch = 16
plot(g, xaxt = 'n', yaxt = 'n', xlab = "Time points",
    ylab = "Spatial features", xlim = c(.5,5.5), ylim = c(.5,3.5),
	pch = pch, col = col)
abline(h=s, col = grey(.8))
abline(v=t, col = grey(.8))
points(g)
axis(1, at = t, labels = c("1st", "2nd", "3rd", "4th"))
axis(2, at = s, labels = c("1st", "2nd", "3rd"))
text(g, labels = 1:12, pos=4)
title("STF: full grid layout")
# 2:
s = 1:3
t = c(1, 2.2, 3, 4.5)
g = data.frame(rep(t, each=3), rep(s,4))
sel = c(1,2,3,5,6,7,11)
plot(g[sel,], xaxt = 'n', yaxt = 'n', xlab = "Time points",
    ylab = "Spatial features", xlim = c(.5,5.5), ylim = c(.5,3.5),
	pch = pch, col = col)
abline(h=s, col = grey(.8))
abline(v=t, col = grey(.8))
points(g[sel,])
axis(1, at = t, labels = c("1st", "2nd", "3rd", "4th"))
axis(2, at = s, labels = c("1st", "2nd", "3rd"))
text(g[sel,], labels = paste(1:length(sel), "[",c(1,2,3,2,3,1,2),",",c(1,1,1,2, 2,3,4),"]", sep=""), pos=4)
title("STS: sparse grid layout")
# 3:
s = c(1,2,3,1,4)
t = c(1, 2.2, 2.5, 4, 4.5)
g = data.frame(t,s)
plot(g, xaxt = 'n', yaxt = 'n', xlab = "Time points",                
    ylab = "Spatial features", xlim = c(.5,5.5), ylim = c(.5,4.5),
	pch = pch, col = col)
#abline(h=s, col = grey(.8))
#abline(v=t, col = grey(.8))
arrows(t,s,0.5,s,.1,col='red')
arrows(t,s,t,0.5,.1,col='red')
points(g)
axis(1, at = sort(unique(t)), labels = c("1st", "2nd", "3rd", "4th", "5th"))
axis(2, at = sort(unique(s)), labels = c("1st,4th", "2nd", "3rd", "5th"))
text(g, labels = 1:5, pos=4)
title("STI: irregular layout")
# 4: traj
ns = 400
nt = 100
s = sort(runif(ns))
t = sort(runif(nt))
g = data.frame(t[1:30],s[1:30])
plot(g, xaxt = 'n', yaxt = 'n', xlab = "Time points",                
    ylab = "Spatial features", 
	type='l', col = 'blue', xlim = c(0,1), ylim = c(0,s[136]))
lines(data.frame(t[41:60],s[31:50]), col = 'blue')
lines(data.frame(t[91:100],s[51:60]), col = 'blue')
lines(data.frame(t[21:40],s[61:80]), col = 'red')
lines(data.frame(t[51:90],s[81:120]), col = 'red')
lines(data.frame(t[11:25],s[121:135]), col = 'green')
#abline(h=s, col = grey(.8))
#abline(v=t, col = grey(.8))
#arrows(t,s,0.5,s,.1,col='red')
#arrows(t,s,t,0.5,.1,col='red')
axis(1, at = sort(unique(t)), labels = rep("", length(t)))
axis(2, at = sort(unique(s)), labels = rep("", length(s)))
#text(g, labels = 1:5, pos=4)
title("STT: trajectory")
opar$cin = opar$cra = opar$csi = opar$cxy = opar$din = opar$page = NULL
par(opar)


###################################################
### code chunk number 6: jss816.Rnw:582-586
###################################################
sp = cbind(x = c(0,0,1), y = c(0,1,1))
row.names(sp) = paste("point", 1:nrow(sp), sep="")
library(sp)
sp = SpatialPoints(sp)


###################################################
### code chunk number 7: jss816.Rnw:590-591
###################################################
time = as.POSIXct("2010-08-05", tz = "GMT")+3600*(10:13)


###################################################
### code chunk number 8: jss816.Rnw:594-598
###################################################
m = c(10,20,30) # means for each of the 3 point locations
values = rnorm(length(sp)*length(time), mean = rep(m, 4))
IDs = paste("ID",1:length(values), sep = "_")
mydata = data.frame(values = signif(values, 3), ID=IDs)


###################################################
### code chunk number 9: jss816.Rnw:601-602
###################################################
library(spacetime)


###################################################
### code chunk number 10: jss816.Rnw:604-606 (eval = FALSE)
###################################################
## library(spacetime)
## stfdf = STFDF(sp, time, data = mydata)


###################################################
### code chunk number 11: jss816.Rnw:611-612
###################################################
stfdf = STFDF(sp, time, mydata, time+60)


###################################################
### code chunk number 12: jss816.Rnw:663-664 (eval = FALSE)
###################################################
## air_quality[2:3, 1:10, "PM10"]


###################################################
### code chunk number 13: jss816.Rnw:668-669 (eval = FALSE)
###################################################
## air_quality[Germany, "2008::2009", "PM10"]


###################################################
### code chunk number 14: jss816.Rnw:688-691
###################################################
xs1 = as(stfdf, "Spatial")
class(xs1)
xs1


###################################################
### code chunk number 15: jss816.Rnw:695-696
###################################################
attr(xs1, "time")


###################################################
### code chunk number 16: jss816.Rnw:700-704
###################################################
x = as(stfdf, "STIDF")
xs2 = as(x, "Spatial")
class(xs2)
xs2[1:4,]


###################################################
### code chunk number 17: jss816.Rnw:751-753 (eval = FALSE)
###################################################
## scales=list(x=list(rot = 45))
## stplot(wind.data, mode = "xt", scales = scales, xlab = NULL)


###################################################
### code chunk number 18: jss816.Rnw:801-811 (eval = FALSE)
###################################################
## # code to create figure 5.
## library(lattice)
## library(RColorBrewer)
## b = brewer.pal(12, "Set3")
## par.settings = list(superpose.symbol = list(col = b, fill = b), 
## 	superpose.line = list(col = b),
## 	fontsize = list(text=9)) 
## stplot(wind.data, mode = "ts",  auto.key=list(space="right"), 
## 	xlab = "1961", ylab = expression(sqrt(speed)),
## 	par.settings = par.settings)


###################################################
### code chunk number 19: jss816.Rnw:853-855
###################################################
library(xts)
.parseISO8601('2010-05')


###################################################
### code chunk number 20: jss816.Rnw:860-861
###################################################
.parseISO8601('2010-05-01T13:30/2010-05-01T13:39')


###################################################
### code chunk number 21: jss816.Rnw:932-936
###################################################
library(maptools)
fname = system.file("shapes/sids.shp", package="maptools")[1]
nc = readShapePoly(fname, 
	proj4string = CRS("+proj=longlat +datum=NAD27 +ellps=clrk66"))


###################################################
### code chunk number 22: jss816.Rnw:939-941
###################################################
time = as.POSIXct(c("1974-07-01", "1979-07-01"), tz = "GMT")
endTime = as.POSIXct(c("1978-06-30", "1984-06-30"), tz = "GMT")


###################################################
### code chunk number 23: jss816.Rnw:944-948
###################################################
data = data.frame(
	BIR = c(nc$BIR74, nc$BIR79),
	NWBIR = c(nc$NWBIR74, nc$NWBIR79),
	SID = c(nc$SID74, nc$SID79))


###################################################
### code chunk number 24: jss816.Rnw:951-952
###################################################
nct = STFDF(sp = as(nc, "SpatialPolygons"), time, data, endTime)


###################################################
### code chunk number 25: jss816.Rnw:961-966
###################################################
library(maps)
states.m = map('state', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(states.m$names, ":"), function(x) x[1])
library(maptools)
states = map2SpatialPolygons(states.m, IDs=IDs)


###################################################
### code chunk number 26: jss816.Rnw:969-971
###################################################
yrs = 1970:1986
time = as.POSIXct(paste(yrs, "-01-01", sep=""), tz = "GMT")


###################################################
### code chunk number 27: jss816.Rnw:974-976
###################################################
library(plm)
data("Produc")


###################################################
### code chunk number 28: jss816.Rnw:982-984
###################################################
# deselect District of Columbia, polygon 8, which is not present in Produc:
Produc.st = STFDF(states[-8], time, Produc[order(Produc[2], Produc[1]),])


###################################################
### code chunk number 29: jss816.Rnw:986-988
###################################################
library(RColorBrewer)
stplot(Produc.st[,,"unemp"], yrs, col.regions = brewer.pal(9, "YlOrRd"),cuts=9)


###################################################
### code chunk number 30: jss816.Rnw:1011-1013
###################################################
zz <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
	data = as.data.frame(Produc.st), index = c("state", "year"))


###################################################
### code chunk number 31: jss816.Rnw:1034-1040
###################################################
library(gstat)
data("wind")
wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates(wind.loc) = ~x+y
proj4string(wind.loc) = "+proj=longlat +datum=WGS84"


###################################################
### code chunk number 32: jss816.Rnw:1045-1050
###################################################
library(mapdata)
plot(wind.loc, xlim = c(-11,-5.4), ylim = c(51,55.5), axes=T, col="red",
	cex.axis =.7)
map("worldHires", add=TRUE, col = grey(.5))
text(coordinates(wind.loc), pos=1, label=wind.loc$Station, cex=.7)


###################################################
### code chunk number 33: jss816.Rnw:1062-1063
###################################################
wind[1:3,]


###################################################
### code chunk number 34: jss816.Rnw:1067-1069
###################################################
wind$time = ISOdate(wind$year+1900, wind$month, wind$day)
wind$jday = as.numeric(format(wind$time, '%j'))


###################################################
### code chunk number 35: jss816.Rnw:1074-1081
###################################################
stations = 4:15
windsqrt = sqrt(0.5148 * as.matrix(wind[stations])) # knots -> m/s
Jday = 1:366
windsqrt = windsqrt - mean(windsqrt)
daymeans = sapply(split(windsqrt, wind$jday), mean)
meanwind = lowess(daymeans ~ Jday, f = 0.1)$y[wind$jday]
velocities = apply(windsqrt, 2, function(x) { x - meanwind })


###################################################
### code chunk number 36: jss816.Rnw:1086-1090
###################################################
wind.loc = wind.loc[match(names(wind[4:15]), wind.loc$Code),]
pts = coordinates(wind.loc[match(names(wind[4:15]), wind.loc$Code),])
rownames(pts) = wind.loc$Station
pts = SpatialPoints(pts, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))


###################################################
### code chunk number 37: jss816.Rnw:1096-1099
###################################################
library(rgdal)
utm29 = CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84")
pts = spTransform(pts, utm29)


###################################################
### code chunk number 38: jss816.Rnw:1103-1106
###################################################
wind.data = stConstruct(velocities, space = list(values = 1:ncol(velocities)), 
	time = wind$time, SpatialObj = pts, interval = TRUE)
class(wind.data)


###################################################
### code chunk number 39: jss816.Rnw:1109-1115
###################################################
library(maptools)
m = map2SpatialLines(
	#map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
	map("worldHires", xlim = c(-11.5,-6.0), ylim = c(51.3,55.0), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84"
m = spTransform(m, utm29)


###################################################
### code chunk number 40: jss816.Rnw:1118-1120
###################################################
grd = SpatialPixels(SpatialPoints(makegrid(m, n = 300)),
	proj4string = proj4string(m))


###################################################
### code chunk number 41: jss816.Rnw:1123-1124
###################################################
wind.data = wind.data[, "1961-04"]


###################################################
### code chunk number 42: jss816.Rnw:1128-1132
###################################################
n = 10
library(xts)
tgrd = seq(min(index(wind.data)), max(index(wind.data)), length=n)
pred.grd = STF(grd, tgrd)


###################################################
### code chunk number 43: jss816.Rnw:1136-1140
###################################################
v = vgmST("separable", space = vgm(1, "Exp", 750000), time = vgm(1, "Exp", 1.5 * 3600 * 24),
         sill=0.6)
wind.ST = krigeST(values ~ 1, wind.data, pred.grd, v)
colnames(wind.ST@data) <- "sqrt_speed"


###################################################
### code chunk number 44: jss816.Rnw:1146-1151 (eval = FALSE)
###################################################
## layout = list(list("sp.lines", m, col='grey'),
## 	list("sp.points", pts, first=F, cex=.5))
## stplot(wind.ST, col.regions=brewer.pal(11, "RdBu")[-c(10,11)],
## 	at=seq(-1.375,1,by=.25),
## 	par.strip.text = list(cex=.7), sp.layout = layout)


###################################################
### code chunk number 45: jss816.Rnw:1154-1161
###################################################
pdf("wind.pdf", height=4.5)
layout = list(list("sp.lines", m, col='grey'),
	list("sp.points", pts, first=F, cex=.5))
print(stplot(wind.ST, col.regions=brewer.pal(11, "RdBu")[-c(10,11)],
	at=seq(-1.375,1,by=.25),
	par.strip.text = list(cex=.7), sp.layout = layout))
dev.off()


###################################################
### code chunk number 46: jss816.Rnw:1164-1175
###################################################
pdf("windts.pdf", height = 4)
library(lattice)
library(RColorBrewer)
b = brewer.pal(12,"Set3")
par.settings = list(superpose.symbol = list(col = b, fill = b), 
	superpose.line = list(col = b),
	fontsize = list(text=9)) 
print(stplot(wind.data, mode = "ts",  auto.key=list(space="right"), 
	xlab = "1961", ylab = expression(sqrt(speed)),
	par.settings = par.settings))
dev.off()


###################################################
### code chunk number 47: jss816.Rnw:1178-1183 (eval = FALSE)
###################################################
## pdf("hov.pdf")
## scales=list(x=list(rot=45))
## stplot(wind.data, mode = "xt", scales = scales, xlab = NULL, 
## 	col.regions=brewer.pal(11, "RdBu"),at = seq(-1.625,1.125,by=.25))
## dev.off()


###################################################
### code chunk number 48: jss816.Rnw:1189-1190 (eval = FALSE)
###################################################
## eof.data = eof(wind.data)


###################################################
### code chunk number 49: jss816.Rnw:1193-1194
###################################################
eof.int = eof(wind.ST)


###################################################
### code chunk number 50: jss816.Rnw:1198-1199 (eval = FALSE)
###################################################
## eof.xts = eof(wind.ST, "temporal")


###################################################
### code chunk number 51: jss816.Rnw:1208-1210
###################################################
print(spplot(eof.int[1:4], col.regions=bpy.colors(),
	par.strip.text = list(cex=.5), as.table = TRUE, sp.layout = layout))


###################################################
### code chunk number 52: jss816.Rnw:1227-1236
###################################################
library(adehabitatLT)
data("puechabonsp")
locs = puechabonsp$relocs
xy = coordinates(locs)
da = as.character(locs$Date)
da = as.POSIXct(strptime(as.character(locs$Date),"%y%m%d", tz = "GMT"))
ltr = as.ltraj(xy, da, id = locs$Name)
foo = function(dt) dt > 100*3600*24
l2 = cutltraj(ltr, "foo(dt)", nextr = TRUE)


###################################################
### code chunk number 53: jss816.Rnw:1240-1242 (eval = FALSE)
###################################################
## sttdf = as(l2, "STTDF")
## stplot(sttdf, by="time*id")


###################################################
### code chunk number 54: jss816.Rnw:1248-1250
###################################################
sttdf = as(l2, "STTDF")
print(stplot(sttdf, by="time*id"))


###################################################
### code chunk number 55: jss816.Rnw:1263-1267
###################################################
library(cshapes)
cs = cshp()
names(cs)
row.names(cs) = paste(as.character(cs$CNTRY_NAME), 1:244)


###################################################
### code chunk number 56: jss816.Rnw:1277-1282
###################################################
cs = cs[cs$COWSYEAR > -1,]
begin = as.POSIXct(paste(cs$COWSYEAR, cs$COWSMONTH, cs$COWSDAY, sep="-"), 
    tz = "GMT")
end = as.POSIXct(paste(cs$COWEYEAR, cs$COWEMONTH,cs$COWEDAY, sep="-"), 
    tz = "GMT")


###################################################
### code chunk number 57: jss816.Rnw:1285-1286
###################################################
st = STIDF(geometry(cs), begin, as.data.frame(cs), end)


###################################################
### code chunk number 58: jss816.Rnw:1290-1292
###################################################
pt = SpatialPoints(cbind(7, 52), CRS(proj4string(cs)))
as.data.frame(st[pt,,1:5])


###################################################
### code chunk number 59: jss816.Rnw:1316-1318 (eval = FALSE)
###################################################
## library(spacetime)
## demo(CressieWikle)


###################################################
### code chunk number 60: jss816.Rnw:1329-1331 (eval = FALSE)
###################################################
## library(gstat)
## vignette("st")


###################################################
### code chunk number 61: jss816.Rnw:1337-1339 (eval = FALSE)
###################################################
## library(spacetime)
## vignette("stpg")


###################################################
### code chunk number 62: jss816.Rnw:1351-1353 (eval = FALSE)
###################################################
## library(spacetime)
## vignette("sto")


