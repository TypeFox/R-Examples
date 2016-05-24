### R code from vignette source 'adehabitatMA.Rnw'

###################################################
### code chunk number 1: adehabitatMA.Rnw:28-32
###################################################
oldopt <- options(width=80, warn=-1)
.PngNo <- 0
wi <- 480
pt <- 20


###################################################
### code chunk number 2: afig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
##           .PngNo, ".png", sep="")
## png(file=file, width = wi, height = wi, pointsize = pt)


###################################################
### code chunk number 3: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 4: zfigg (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 5: adehabitatMA.Rnw:103-104
###################################################
library(adehabitatMA)


###################################################
### code chunk number 6: adehabitatMA.Rnw:107-109
###################################################
set.seed(13431)
adeoptions(shortprint=TRUE)


###################################################
### code chunk number 7: adehabitatMA.Rnw:129-130
###################################################
data(lynxjura)


###################################################
### code chunk number 8: adehabitatMA.Rnw:142-143
###################################################
head(lynxjura$locs)


###################################################
### code chunk number 9: adehabitatMA.Rnw:150-151
###################################################
lynxjura$map


###################################################
### code chunk number 10: figu1 (eval = FALSE)
###################################################
## mimage(lynxjura$map)


###################################################
### code chunk number 11: adehabitatMA.Rnw:161-162 (eval = FALSE)
###################################################
## mimage(lynxjura$map)


###################################################
### code chunk number 12: adehabitatMA.Rnw:166-169
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
mimage(lynxjura$map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 13: adehabitatMA.Rnw:186-187 (eval = FALSE)
###################################################
## explore(lynxjura$map)


###################################################
### code chunk number 14: adehabitatMA.Rnw:198-199 (eval = FALSE)
###################################################
## explore(lynxjura$map, panel.last=function() points(lynxjura$locs, pch=3))


###################################################
### code chunk number 15: adehabitatMA.Rnw:209-210
###################################################
map <- lynxjura$map


###################################################
### code chunk number 16: sodssss (eval = FALSE)
###################################################
## hist(map)


###################################################
### code chunk number 17: adehabitatMA.Rnw:219-220 (eval = FALSE)
###################################################
## hist(map)


###################################################
### code chunk number 18: adehabitatMA.Rnw:224-227
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
hist(map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 19: qsfddfdd (eval = FALSE)
###################################################
## forest <- map[,1]
## forest[[1]][forest[[1]]<95] <- NA
## image(forest, col="green")


###################################################
### code chunk number 20: adehabitatMA.Rnw:242-243 (eval = FALSE)
###################################################
## forest <- map[,1]
## forest[[1]][forest[[1]]<95] <- NA
## image(forest, col="green")


###################################################
### code chunk number 21: adehabitatMA.Rnw:248-251
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
forest <- map[,1]
forest[[1]][forest[[1]]<95] <- NA
image(forest, col="green")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 22: ksdkksss (eval = FALSE)
###################################################
## lab <- labcon(forest)
## image(lab)


###################################################
### code chunk number 23: adehabitatMA.Rnw:265-266 (eval = FALSE)
###################################################
## lab <- labcon(forest)
## image(lab)


###################################################
### code chunk number 24: adehabitatMA.Rnw:270-273
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
lab <- labcon(forest)
image(lab)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 25: adehabitatMA.Rnw:282-284
###################################################
lab
max(lab[[1]])


###################################################
### code chunk number 26: adehabitatMA.Rnw:290-291
###################################################
gridparameters(lab)


###################################################
### code chunk number 27: adehabitatMA.Rnw:296-297
###################################################
table(lab[[1]])*500*500


###################################################
### code chunk number 28: adehabitatMA.Rnw:308-310
###################################################
fullgrid(lab) <- TRUE
fullgrid(map) <- TRUE


###################################################
### code chunk number 29: adehabitatMA.Rnw:315-316
###################################################
mean(map[[2]][lab[[1]]==1], na.rm=TRUE)


###################################################
### code chunk number 30: sldsss (eval = FALSE)
###################################################
## comp1 <- map[2]
## comp1[[1]][map[[1]]<95] <- NA
## comp1[[1]][lab[[1]]!=1] <- NA
## image(comp1)


###################################################
### code chunk number 31: adehabitatMA.Rnw:329-330 (eval = FALSE)
###################################################
## comp1 <- map[2]
## comp1[[1]][map[[1]]<95] <- NA
## comp1[[1]][lab[[1]]!=1] <- NA
## image(comp1)


###################################################
### code chunk number 32: adehabitatMA.Rnw:334-337
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
comp1 <- map[2]
comp1[[1]][map[[1]]<95] <- NA
comp1[[1]][lab[[1]]!=1] <- NA
image(comp1)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 33: skks (eval = FALSE)
###################################################
## image(forest, col="red")


###################################################
### code chunk number 34: adehabitatMA.Rnw:351-352 (eval = FALSE)
###################################################
## image(forest, col="red")


###################################################
### code chunk number 35: adehabitatMA.Rnw:356-359
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(forest, col="red")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 36: adehabitatMA.Rnw:366-367
###################################################
con <- getcontour(forest)


###################################################
### code chunk number 37: ssskkkq (eval = FALSE)
###################################################
## plot(con, col="green")


###################################################
### code chunk number 38: adehabitatMA.Rnw:394-395 (eval = FALSE)
###################################################
## plot(con, col="green")


###################################################
### code chunk number 39: adehabitatMA.Rnw:399-402
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(con, col="green")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 40: adehabitatMA.Rnw:425-427
###################################################
for1 <- morphology(forest, "dilate", nt=1)
for1


###################################################
### code chunk number 41: sdkskss (eval = FALSE)
###################################################
## image(for1, col="blue")
## image(forest, col="yellow", add=TRUE)


###################################################
### code chunk number 42: adehabitatMA.Rnw:439-440 (eval = FALSE)
###################################################
## image(for1, col="blue")
## image(forest, col="yellow", add=TRUE)


###################################################
### code chunk number 43: adehabitatMA.Rnw:444-447
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(for1, col="blue")
image(forest, col="yellow", add=TRUE)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 44: skkssss (eval = FALSE)
###################################################
## plot(getcontour(for1), col="green")


###################################################
### code chunk number 45: adehabitatMA.Rnw:459-460 (eval = FALSE)
###################################################
## plot(getcontour(for1), col="green")


###################################################
### code chunk number 46: adehabitatMA.Rnw:464-467
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(getcontour(for1), col="green")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 47: skdskdk (eval = FALSE)
###################################################
## map <- lynxjura$map
## mimage(map)


###################################################
### code chunk number 48: adehabitatMA.Rnw:491-492 (eval = FALSE)
###################################################
## map <- lynxjura$map
## mimage(map)


###################################################
### code chunk number 49: adehabitatMA.Rnw:496-499
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
map <- lynxjura$map
mimage(map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 50: adehabitatMA.Rnw:505-506
###################################################
gridparameters(map)


###################################################
### code chunk number 51: sssshhh (eval = FALSE)
###################################################
## mimage(lowres(map, 10))


###################################################
### code chunk number 52: adehabitatMA.Rnw:519-520 (eval = FALSE)
###################################################
## mimage(lowres(map, 10))


###################################################
### code chunk number 53: adehabitatMA.Rnw:524-527
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
mimage(lowres(map, 10))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 54: sqkdksssq (eval = FALSE)
###################################################
## map[[1]] <- as.numeric(cut(map[[1]],3))
## image(map, 1)


###################################################
### code chunk number 55: adehabitatMA.Rnw:542-543 (eval = FALSE)
###################################################
## map[[1]] <- as.numeric(cut(map[[1]],3))
## image(map, 1)


###################################################
### code chunk number 56: adehabitatMA.Rnw:547-550
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
map[[1]] <- as.numeric(cut(map[[1]],3))
image(map, 1)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 57: sdskkkk (eval = FALSE)
###################################################
## image(lowres(map, 10, which.fac=1))


###################################################
### code chunk number 58: adehabitatMA.Rnw:566-567 (eval = FALSE)
###################################################
## image(lowres(map, 10, which.fac=1))


###################################################
### code chunk number 59: adehabitatMA.Rnw:571-574
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(lowres(map, 10, which.fac=1))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 60: skssjdx (eval = FALSE)
###################################################
## image(forest, col="green")
## box()


###################################################
### code chunk number 61: adehabitatMA.Rnw:589-590 (eval = FALSE)
###################################################
## image(forest, col="green")
## box()


###################################################
### code chunk number 62: adehabitatMA.Rnw:594-597
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(forest, col="green")
box()
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 63: adehabitatMA.Rnw:604-605 (eval = FALSE)
###################################################
## for2 <- subsetmap(forest)


###################################################
### code chunk number 64: adehabitatMA.Rnw:611-613
###################################################
for2 <- subsetmap(forest, xlim=c(850254.2, 878990.2),
                  ylim=c(2128744, 2172175))


###################################################
### code chunk number 65: skksks (eval = FALSE)
###################################################
## image(for2, col="green")
## box()


###################################################
### code chunk number 66: adehabitatMA.Rnw:621-622 (eval = FALSE)
###################################################
## image(for2, col="green")
## box()


###################################################
### code chunk number 67: adehabitatMA.Rnw:626-629
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(for2, col="green")
box()
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 68: adehabitatMA.Rnw:665-670
###################################################
data(lynxjura)
map <- lynxjura$map
class(map)
locs <- lynxjura$loc
class(locs)


###################################################
### code chunk number 69: flkqskqfkjc (eval = FALSE)
###################################################
## image(map, 1)
## points(locs, pch=3)


###################################################
### code chunk number 70: adehabitatMA.Rnw:680-681 (eval = FALSE)
###################################################
## image(map, 1)
## points(locs, pch=3)


###################################################
### code chunk number 71: adehabitatMA.Rnw:685-688
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(map, 1)
points(locs, pch=3)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 72: adehabitatMA.Rnw:694-695
###################################################
cp <- count.points(locs, map)


###################################################
### code chunk number 73: ssckcc (eval = FALSE)
###################################################
## image(cp)


###################################################
### code chunk number 74: adehabitatMA.Rnw:716-717 (eval = FALSE)
###################################################
## image(cp)


###################################################
### code chunk number 75: adehabitatMA.Rnw:721-724
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(cp)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 76: adehabitatMA.Rnw:737-738
###################################################
head(locs)


###################################################
### code chunk number 77: adehabitatMA.Rnw:748-750
###################################################
cpr <- count.points(locs[,"Type"], map)
cpr


###################################################
### code chunk number 78: sdkkckck (eval = FALSE)
###################################################
## mimage(cpr)


###################################################
### code chunk number 79: adehabitatMA.Rnw:759-760 (eval = FALSE)
###################################################
## mimage(cpr)


###################################################
### code chunk number 80: adehabitatMA.Rnw:764-767
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
mimage(cpr)
dev.null <- dev.off()
cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 81: adehabitatMA.Rnw:781-782
###################################################
df <- join(locs, map)


###################################################
### code chunk number 82: adehabitatMA.Rnw:788-789
###################################################
head(df)


###################################################
### code chunk number 83: adehabitatMA.Rnw:808-810
###################################################
asc <- ascgen(locs, cellsize=5000)
asc


###################################################
### code chunk number 84: llldcvdv (eval = FALSE)
###################################################
## image(asc)


###################################################
### code chunk number 85: adehabitatMA.Rnw:821-822 (eval = FALSE)
###################################################
## image(asc)


###################################################
### code chunk number 86: adehabitatMA.Rnw:826-829
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(asc)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 87: adehabitatMA.Rnw:842-843
###################################################
po <- locs[locs[["Type"]]=="O",]


###################################################
### code chunk number 88: sfjkfc (eval = FALSE)
###################################################
## image(buffer(po, map, 3000))


###################################################
### code chunk number 89: adehabitatMA.Rnw:853-854 (eval = FALSE)
###################################################
## image(buffer(po, map, 3000))


###################################################
### code chunk number 90: adehabitatMA.Rnw:858-861
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(buffer(po, map, 3000))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 91: sdskwkl (eval = FALSE)
###################################################
## plot(con)


###################################################
### code chunk number 92: adehabitatMA.Rnw:872-873 (eval = FALSE)
###################################################
## plot(con)


###################################################
### code chunk number 93: adehabitatMA.Rnw:877-880
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(con)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 94: dfkwkfckj (eval = FALSE)
###################################################
## image(buffer(con, map, 3000))


###################################################
### code chunk number 95: adehabitatMA.Rnw:891-892 (eval = FALSE)
###################################################
## image(buffer(con, map, 3000))


###################################################
### code chunk number 96: adehabitatMA.Rnw:896-899
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(buffer(con, map, 3000))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 97: qksdk (eval = FALSE)
###################################################
## sl <- as(con, "SpatialLines")
## image(buffer(sl, map, 500))


###################################################
### code chunk number 98: adehabitatMA.Rnw:913-914 (eval = FALSE)
###################################################
## sl <- as(con, "SpatialLines")
## image(buffer(sl, map, 500))


###################################################
### code chunk number 99: adehabitatMA.Rnw:918-921
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
sl <- as(con, "SpatialLines")
image(buffer(sl, map, 500))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 100: adehabitatMA.Rnw:958-959
###################################################
options(oldopt)


