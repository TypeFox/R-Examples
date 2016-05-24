###################################################
### chunk number 1: 
###################################################
owidth <- getOption("width")
options("width"=80)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0
argX11<-formals(x11)


###################################################
### chunk number 2: afig eval=FALSE
###################################################
## if (!is.null(argX11$colortype)) {
##   graphics.off()
##   x11(colortype="gray")
## }
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
## png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
## opar <- par(no.readonly = TRUE)


###################################################
### chunk number 3: zfig eval=FALSE
###################################################
## par(opar)
## dev.null <- dev.off()
## cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 4: zfigkasc eval=FALSE
###################################################
## par(opar)
## dev.null <- dev.off()
## cat("\\includegraphics[height=12cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 5: 
###################################################
library(adehabitat)


###################################################
### chunk number 6: rand eval=FALSE
###################################################
## mat <- matrix(rnorm(10000), 100, 100)
## asc <- as.asc(mat)
## image(asc)
## box()


###################################################
### chunk number 7:  eval=FALSE
###################################################
## mat <- matrix(rnorm(10000), 100, 100)
## asc <- as.asc(mat)
## image(asc)
## box()


###################################################
### chunk number 8: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
mat <- matrix(rnorm(10000), 100, 100)
asc <- as.asc(mat)
image(asc)
box()
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 9: 
###################################################
(path.to.file <- paste(system.file(package = "adehabitat"), 
                       "ascfiles/elevation.asc", sep = "/"))


###################################################
### chunk number 10: elev eval=FALSE
###################################################
## el <- import.asc(path.to.file)
## image(el, main = "Elevation")


###################################################
### chunk number 11:  eval=FALSE
###################################################
## el <- import.asc(path.to.file)
## image(el, main = "Elevation")


###################################################
### chunk number 12: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
el <- import.asc(path.to.file)
image(el, main = "Elevation")
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 13: 
###################################################
(path.to.map <- paste( system.file(package = "adehabitat"), 
                      "ascfiles/aspect.asc", sep = "/"))
asp <- import.asc(path.to.map, type = "factor")
levels(asp)


###################################################
### chunk number 14: 
###################################################
asp <- import.asc(path.to.map, type = "factor",
                  lev = c("North", "East", "West", "South"))
levels(asp)


###################################################
### chunk number 15: 
###################################################
(path.to.table <- paste(system.file(package = "adehabitat"), 
                        "ascfiles/aspect.txt", sep = "/"))


###################################################
### chunk number 16: 
###################################################
asp <- import.asc(path.to.map, type = "factor", lev = 
                  path.to.table, levnb = 1, labnb = 3)
levels(asp)
co <- colasc(asp, North = "blue", East = "yellow", 
             West = "orange", South = "red")


###################################################
### chunk number 17: asp eval=FALSE
###################################################
## image(asp, clfac = co)
## legend(696662, 3166028, legend = levels(asp), fill = co)


###################################################
### chunk number 18:  eval=FALSE
###################################################
## image(asp, clfac = co)
## legend(696662, 3166028, legend = levels(asp), fill = co)


###################################################
### chunk number 19: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(asp, clfac = co)
legend(696662, 3166028, legend = levels(asp), fill = co)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 20: kasc eval=FALSE
###################################################
## data(puechabon)
## kasc <- puechabon$kasc
## image(kasc)


###################################################
### chunk number 21:  eval=FALSE
###################################################
## data(puechabon)
## kasc <- puechabon$kasc
## image(kasc)


###################################################
### chunk number 22: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
data(puechabon)
kasc <- puechabon$kasc
image(kasc)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=12cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 23: 
###################################################
(obj <- as.kasc(list(Elevation = el, Aspect = asp)))


###################################################
### chunk number 24: 
###################################################
data(puechabon)
puechabon$kasc
puechabon$locs[1:4,]


###################################################
### chunk number 25: prespuech eval=FALSE
###################################################
## el <- getkasc(puechabon$kasc, "Elevation")
## opar <- par(mfrow = c(2,2), mar=c(0,0,4,0))
## for (i in levels(puechabon$locs$Name)) {
##   image(el, 
##         main = paste("Wild boar named", i),
##         axes=FALSE)
##   points(puechabon$locs[puechabon$locs$Name==i,c("X","Y")], pch=16)
## }
## par(opar)


###################################################
### chunk number 26:  eval=FALSE
###################################################
## el <- getkasc(puechabon$kasc, "Elevation")
## opar <- par(mfrow = c(2,2), mar=c(0,0,4,0))
## for (i in levels(puechabon$locs$Name)) {
##   image(el, 
##         main = paste("Wild boar named", i),
##         axes=FALSE)
##   points(puechabon$locs[puechabon$locs$Name==i,c("X","Y")], pch=16)
## }
## par(opar)


###################################################
### chunk number 27: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
el <- getkasc(puechabon$kasc, "Elevation")
opar <- par(mfrow = c(2,2), mar=c(0,0,4,0))
for (i in levels(puechabon$locs$Name)) {
  image(el, 
        main = paste("Wild boar named", i),
        axes=FALSE)
  points(puechabon$locs[puechabon$locs$Name==i,c("X","Y")], pch=16)
}
par(opar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=12cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 28: 
###################################################
data(chamois)
chamois$map
chamois$locs[1:4,]


###################################################
### chunk number 29: preschart eval=FALSE
###################################################
## sl <- getkasc(chamois$map, "Slope")
## image(sl, main = "Distribution of chamois occurrences in the Chartreuse mountain")
## points(chamois$locs, pch=16)


###################################################
### chunk number 30:  eval=FALSE
###################################################
## sl <- getkasc(chamois$map, "Slope")
## image(sl, main = "Distribution of chamois occurrences in the Chartreuse mountain")
## points(chamois$locs, pch=16)


###################################################
### chunk number 31: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
sl <- getkasc(chamois$map, "Slope")
image(sl, main = "Distribution of chamois occurrences in the Chartreuse mountain")
points(chamois$locs, pch=16)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 32: 
###################################################
kasc <- puechabon$kasc
(el <- getkasc(kasc, "Elevation"))


###################################################
### chunk number 33: disfm eval=FALSE
###################################################
## image(distfacmap(getkasc(puechabon$kasc, "Aspect")))


###################################################
### chunk number 34:  eval=FALSE
###################################################
## image(distfacmap(getkasc(puechabon$kasc, "Aspect")))


###################################################
### chunk number 35: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(distfacmap(getkasc(puechabon$kasc, "Aspect")))
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 36: 
###################################################
er8 <- morphology(el, operation="erode", nt=8)
di8 <- morphology(el, operation="dilate", nt=8)


###################################################
### chunk number 37: morpho eval=FALSE
###################################################
## image(di8, col="black")
## image(el, col="gray", add=TRUE)
## image(er8, col="white", add=TRUE)
## 
## arrows(704295, 3159355, 706588, 3157294, col="red", lwd = 2, code = 1)
## text(706240, 3156738, "Boundary of the study area")


###################################################
### chunk number 38:  eval=FALSE
###################################################
## image(di8, col="black")
## image(el, col="gray", add=TRUE)
## image(er8, col="white", add=TRUE)
## 
## arrows(704295, 3159355, 706588, 3157294, col="red", lwd = 2, code = 1)
## text(706240, 3156738, "Boundary of the study area")


###################################################
### chunk number 39: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(di8, col="black")
image(el, col="gray", add=TRUE)
image(er8, col="white", add=TRUE)

arrows(704295, 3159355, 706588, 3157294, col="red", lwd = 2, code = 1)
text(706240, 3156738, "Boundary of the study area")
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 40: 
###################################################
data(puechabon)
puechabon$locs[1:4,]


###################################################
### chunk number 41: ptsel eval=FALSE
###################################################
## image(el)
## points(puechabon$locs[,c("X","Y")], pch = 16)


###################################################
### chunk number 42:  eval=FALSE
###################################################
## image(el)
## points(puechabon$locs[,c("X","Y")], pch = 16)


###################################################
### chunk number 43: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(el)
points(puechabon$locs[,c("X","Y")], pch = 16)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 44: buffel eval=FALSE
###################################################
## bu <- buffer(puechabon$locs[,c("X","Y")], el, 500)
## image(bu)
## points(puechabon$locs[,c("X","Y")], pch = 16)


###################################################
### chunk number 45:  eval=FALSE
###################################################
## bu <- buffer(puechabon$locs[,c("X","Y")], el, 500)
## image(bu)
## points(puechabon$locs[,c("X","Y")], pch = 16)


###################################################
### chunk number 46: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
bu <- buffer(puechabon$locs[,c("X","Y")], el, 500)
image(bu)
points(puechabon$locs[,c("X","Y")], pch = 16)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 47: 
###################################################
bubis <- bu * el
mean(as.vector(bubis), na.rm = TRUE)
sd(as.vector(bubis), na.rm = TRUE)


###################################################
### chunk number 48: bufbis eval=FALSE
###################################################
## image(bubis)


###################################################
### chunk number 49:  eval=FALSE
###################################################
## image(bubis)


###################################################
### chunk number 50: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(bubis)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 51: 
###################################################
vec <- join.asc(puechabon$locs[,c("X", "Y")], el)
length(vec)
nrow(puechabon$locs)
vec[1:10]


###################################################
### chunk number 52: 
###################################################
df <- join.kasc(puechabon$locs[,c("X", "Y")], puechabon$kasc)
nrow(df)
nrow(puechabon$locs)
df[1:10,]


###################################################
### chunk number 53: 
###################################################
(cp <- count.points(puechabon$locs[,c("X","Y")],  el))


###################################################
### chunk number 54: countpoints eval=FALSE
###################################################
## image(cp)
## box()


###################################################
### chunk number 55:  eval=FALSE
###################################################
## image(cp)
## box()


###################################################
### chunk number 56: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(cp)
box()
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 57: 
###################################################
(cp <- count.points.id(puechabon$locs[,c("X","Y")],  puechabon$locs$Name, el))


###################################################
### chunk number 58: cpid eval=FALSE
###################################################
## image(cp)


###################################################
### chunk number 59:  eval=FALSE
###################################################
## image(cp)


###################################################
### chunk number 60: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(cp)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 61: ascgen eval=FALSE
###################################################
## hihi <- ascgen(xy = puechabon$locs[,c("X","Y")], cellsize = 500)
## image(hihi)
## box()


###################################################
### chunk number 62:  eval=FALSE
###################################################
## hihi <- ascgen(xy = puechabon$locs[,c("X","Y")], cellsize = 500)
## image(hihi)
## box()


###################################################
### chunk number 63: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
hihi <- ascgen(xy = puechabon$locs[,c("X","Y")], cellsize = 500)
image(hihi)
box()
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 64: cpidvg eval=FALSE
###################################################
## tmpbis <- count.points.id(xy = puechabon$locs[,c("X","Y")], id = puechabon$locs$Name, hihi)
## image(tmpbis)


###################################################
### chunk number 65:  eval=FALSE
###################################################
## tmpbis <- count.points.id(xy = puechabon$locs[,c("X","Y")], id = puechabon$locs$Name, hihi)
## image(tmpbis)


###################################################
### chunk number 66: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
tmpbis <- count.points.id(xy = puechabon$locs[,c("X","Y")], id = puechabon$locs$Name, hihi)
image(tmpbis)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 67: 
###################################################
el <- getkasc(puechabon$kasc, "Elevation")
elcat <- el < 200
class(elcat)
names(attributes(elcat))


###################################################
### chunk number 68: 
###################################################
(elcat<-getascattr(el, elcat, type = "factor", lev = c("> 200 m", "< 200 m")))


###################################################
### chunk number 69: elcat eval=FALSE
###################################################
## image(elcat)
## legend(698000, 3165000, levels(elcat), fill=rainbow(2))


###################################################
### chunk number 70:  eval=FALSE
###################################################
## image(elcat)
## legend(698000, 3165000, levels(elcat), fill=rainbow(2))


###################################################
### chunk number 71: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(elcat)
legend(698000, 3165000, levels(elcat), fill=rainbow(2))
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 72: managna eval=FALSE
###################################################
## kasc <- puechabon$kasc
## el <- getkasc(kasc, "Elevation")
## sl <- getkasc(kasc, "Slope")
## el[el < 200] <- NA
## tmp <- as.kasc(list(Elevation = el, Slope = sl))
## image(tmp)


###################################################
### chunk number 73:  eval=FALSE
###################################################
## kasc <- puechabon$kasc
## el <- getkasc(kasc, "Elevation")
## sl <- getkasc(kasc, "Slope")
## el[el < 200] <- NA
## tmp <- as.kasc(list(Elevation = el, Slope = sl))
## image(tmp)


###################################################
### chunk number 74: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
kasc <- puechabon$kasc
el <- getkasc(kasc, "Elevation")
sl <- getkasc(kasc, "Slope")
el[el < 200] <- NA
tmp <- as.kasc(list(Elevation = el, Slope = sl))
image(tmp)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 75: manna eval=FALSE
###################################################
## tmp <- managNAkasc(tmp)
## image(tmp)


###################################################
### chunk number 76:  eval=FALSE
###################################################
## tmp <- managNAkasc(tmp)
## image(tmp)


###################################################
### chunk number 77: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
tmp <- managNAkasc(tmp)
image(tmp)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 78: 
###################################################
data(puechabon)
kasc <- puechabon$kasc
toto <- kasc[1:10,]
class(toto) <- "data.frame"
toto


###################################################
### chunk number 79: 
###################################################
huhu <- kasc2df(kasc)
names(huhu)
huhu$index[1:4]
huhu$tab[1:4,]


###################################################
### chunk number 80: 
###################################################
huhu$tab$Aspect <- NULL
(pc <- dudi.pca(huhu$tab, scannf =FALSE, nf=2))


###################################################
### chunk number 81: df2kasc eval=FALSE
###################################################
## map <- df2kasc(pc$li, huhu$index, kasc)
## image(map)


###################################################
### chunk number 82:  eval=FALSE
###################################################
## map <- df2kasc(pc$li, huhu$index, kasc)
## image(map)


###################################################
### chunk number 83: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
map <- df2kasc(pc$li, huhu$index, kasc)
image(map)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=12cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 84: 
###################################################
(kasc <- chamois$map)
(si1 <- object.size(kasc))


###################################################
### chunk number 85: donncham eval=FALSE
###################################################
## image(kasc)


###################################################
### chunk number 86:  eval=FALSE
###################################################
## image(kasc)


###################################################
### chunk number 87: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(kasc)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=12cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 88: 
###################################################
(m <- lowres(kasc, np = 4))


###################################################
### chunk number 89: lowres eval=FALSE
###################################################
## image(m)


###################################################
### chunk number 90:  eval=FALSE
###################################################
## image(m)


###################################################
### chunk number 91: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(m)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=12cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 92: 
###################################################
(si2 <- object.size(m))
(si1 - si2)/si1


###################################################
### chunk number 93: subset eval=FALSE
###################################################
## data(chamois)
## slope <- getkasc(chamois$map, "Slope")
## def.par <- par(no.readonly = TRUE)
## layout(matrix(c(1,1,1,1,1,1,1,1,2), ncol = 3, byrow = TRUE))
## par(mar = c(0,0,0,0))
## image(slope, axes=FALSE)
## box()
## 
## x <- c(863603.8, 867286.5)
## y <- c(2042689, 2045797)
## polygon(x = c(x[1], x[2], x[2], x[1]),
##         y = c(y[1], y[1], y[2], y[2]), lwd=2)
## 
## sl2 <- subsetmap(slope, xlim = x, ylim = y)
## par(mar = c(0,0,2,0))
## image(sl2, axes = FALSE, main = "Reduced map")
## box()
## par(def.par)


###################################################
### chunk number 94:  eval=FALSE
###################################################
## data(chamois)
## slope <- getkasc(chamois$map, "Slope")
## def.par <- par(no.readonly = TRUE)
## layout(matrix(c(1,1,1,1,1,1,1,1,2), ncol = 3, byrow = TRUE))
## par(mar = c(0,0,0,0))
## image(slope, axes=FALSE)
## box()
## 
## x <- c(863603.8, 867286.5)
## y <- c(2042689, 2045797)
## polygon(x = c(x[1], x[2], x[2], x[1]),
##         y = c(y[1], y[1], y[2], y[2]), lwd=2)
## 
## sl2 <- subsetmap(slope, xlim = x, ylim = y)
## par(mar = c(0,0,2,0))
## image(sl2, axes = FALSE, main = "Reduced map")
## box()
## par(def.par)


###################################################
### chunk number 95: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
data(chamois)
slope <- getkasc(chamois$map, "Slope")
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,1,1,1,1,1,1,1,2), ncol = 3, byrow = TRUE))
par(mar = c(0,0,0,0))
image(slope, axes=FALSE)
box()

x <- c(863603.8, 867286.5)
y <- c(2042689, 2045797)
polygon(x = c(x[1], x[2], x[2], x[1]),
        y = c(y[1], y[1], y[2], y[2]), lwd=2)

sl2 <- subsetmap(slope, xlim = x, ylim = y)
par(mar = c(0,0,2,0))
image(sl2, axes = FALSE, main = "Reduced map")
box()
par(def.par)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 96: 
###################################################
data(elec88)
ar <- elec88$area
ar[1:5,]


###################################################
### chunk number 97: area eval=FALSE
###################################################
## ar <- as.area(ar)
## plot(ar)


###################################################
### chunk number 98:  eval=FALSE
###################################################
## ar <- as.area(ar)
## plot(ar)


###################################################
### chunk number 99: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
ar <- as.area(ar)
plot(ar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 100: 
###################################################
data(puechabon)
lo <- puechabon$locs
cp <- mcp(lo[,c("X", "Y")], lo[,"Name"])
class(cp)


###################################################
### chunk number 101: convpol eval=FALSE
###################################################
## opar <- par(mar=c(0,0,0,0))
## plot(cp, colp=NULL)
## points(puechabon$locs[,c("X", "Y")], pch=16, col = as.numeric(puechabon$locs$Name))
## box()
## par(opar)


###################################################
### chunk number 102:  eval=FALSE
###################################################
## opar <- par(mar=c(0,0,0,0))
## plot(cp, colp=NULL)
## points(puechabon$locs[,c("X", "Y")], pch=16, col = as.numeric(puechabon$locs$Name))
## box()
## par(opar)


###################################################
### chunk number 103: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
opar <- par(mar=c(0,0,0,0))
plot(cp, colp=NULL)
points(puechabon$locs[,c("X", "Y")], pch=16, col = as.numeric(puechabon$locs$Name))
box()
par(opar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 104: 
###################################################
el <- getkasc(puechabon$kasc, "Elevation")
cont.el <- getcontour(el)
class(cont.el)
nlevels(cont.el[,1])


###################################################
### chunk number 105: getcontour eval=FALSE
###################################################
## image(el)
## polygon(cont.el[,2:3], lwd = 3)


###################################################
### chunk number 106:  eval=FALSE
###################################################
## image(el)
## polygon(cont.el[,2:3], lwd = 3)


###################################################
### chunk number 107: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
image(el)
polygon(cont.el[,2:3], lwd = 3)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 108: 
###################################################
lo <- puechabon$locs
kasc <- puechabon$kasc
cp <- mcp(lo[,c("X", "Y")], lo[,"Name"])


###################################################
### chunk number 109: 
###################################################
(rast <- hr.rast(cp, kasc))


###################################################
### chunk number 110: mcprast eval=FALSE
###################################################
## def.par <- par(no.readonly = TRUE)
## layout(matrix(c(1,1,2,4,3,5),2,3))
## par(mar=c(0,0,4,0))
## plot(cp, colp=NULL)
## points(puechabon$locs[,c("X", "Y")], pch=16, col = as.numeric(puechabon$locs$Name))
## box()
## for (i in names(rast)) {
## image(getkasc(rast,i), main = paste("Wild boar named", i), axes=FALSE)
## polygon(cont.el[,2:3])
## box()
## }
## par(def.par)


###################################################
### chunk number 111:  eval=FALSE
###################################################
## def.par <- par(no.readonly = TRUE)
## layout(matrix(c(1,1,2,4,3,5),2,3))
## par(mar=c(0,0,4,0))
## plot(cp, colp=NULL)
## points(puechabon$locs[,c("X", "Y")], pch=16, col = as.numeric(puechabon$locs$Name))
## box()
## for (i in names(rast)) {
## image(getkasc(rast,i), main = paste("Wild boar named", i), axes=FALSE)
## polygon(cont.el[,2:3])
## box()
## }
## par(def.par)


###################################################
### chunk number 112: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,1,2,4,3,5),2,3))
par(mar=c(0,0,4,0))
plot(cp, colp=NULL)
points(puechabon$locs[,c("X", "Y")], pch=16, col = as.numeric(puechabon$locs$Name))
box()
for (i in names(rast)) {
image(getkasc(rast,i), main = paste("Wild boar named", i), axes=FALSE)
polygon(cont.el[,2:3])
box()
}
par(def.par)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=12cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 113: polmask eval=FALSE
###################################################
## el <- getkasc(puechabon$kasc, "Elevation")
## pol <- data.frame(x = c(700658, 699222, 698342, 698643, 700427, 701029),
##                   y = c(3160768, 3160676, 3159402, 3158336, 3158869, 3159657))
## image(el)
## polygon(pol, lwd=2)


###################################################
### chunk number 114:  eval=FALSE
###################################################
## el <- getkasc(puechabon$kasc, "Elevation")
## pol <- data.frame(x = c(700658, 699222, 698342, 698643, 700427, 701029),
##                   y = c(3160768, 3160676, 3159402, 3158336, 3158869, 3159657))
## image(el)
## polygon(pol, lwd=2)


###################################################
### chunk number 115: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
el <- getkasc(puechabon$kasc, "Elevation")
pol <- data.frame(x = c(700658, 699222, 698342, 698643, 700427, 701029),
                  y = c(3160768, 3160676, 3159402, 3158336, 3158869, 3159657))
image(el)
polygon(pol, lwd=2)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 116: mask eval=FALSE
###################################################
## pr <- mcp.rast(pol, el)
## masked.kasc <- setmask(puechabon$kasc, pr)
## image(masked.kasc, xlim = c(696999, 702373), 
##       ylim = c(3156784, 3162297))


###################################################
### chunk number 117:  eval=FALSE
###################################################
## pr <- mcp.rast(pol, el)
## masked.kasc <- setmask(puechabon$kasc, pr)
## image(masked.kasc, xlim = c(696999, 702373), 
##       ylim = c(3156784, 3162297))


###################################################
### chunk number 118: 
###################################################
if (!is.null(argX11$colortype)) {
  graphics.off()
  x11(colortype="gray")
}
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".png", sep="")
png(file=file, width = 700, height = 700, pointsize = 12, bg = "white")
opar <- par(no.readonly = TRUE)
pr <- mcp.rast(pol, el)
masked.kasc <- setmask(puechabon$kasc, pr)
image(masked.kasc, xlim = c(696999, 702373), 
      ylim = c(3156784, 3162297))
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,keepaspectratio]{", file, "}\n\n", sep="")


###################################################
### chunk number 119: 
###################################################
def.pol <- function(x) {
  toto<-locator(1)
  for (i in 2:x) {
    tutu<-locator(1)
    toto$x<-c(toto$x, tutu$x)
    toto$y<-c(toto$y, tutu$y)
    lines(toto$x, toto$y)
  }
  polygon(toto)
  return(toto)
}


