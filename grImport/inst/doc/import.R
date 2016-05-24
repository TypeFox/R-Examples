### R code from vignette source 'import.Rnw'

###################################################
### code chunk number 1: import.Rnw:56-64
###################################################
options(prompt="R> ")
options(continue = "+  ")
options(width = 60)
options(useFancyQuotes = FALSE)
strOptions(strict.width = TRUE)
library(grid)
library(lattice)



###################################################
### code chunk number 2: import.Rnw:83-111
###################################################
chess <- read.table("chessmod.txt", sep = ":", quote = "",
                    col.names = c("player1", "player2",
                      "result", "moves", "year", "place",
                      "openingDetailed"))
chess$result <- factor(ifelse(chess$result == "1-1",
                              "draw",
                              ifelse((chess$result == "1-0" &
                                      chess$player1 == "La Bourdonnais") |
                                     (chess$result == "0-1" &
                                      chess$player2 == "La Bourdonnais"),
                                     "win",
                                     "loss")),
                       levels = c("win", "draw", "loss"))
chess$opening <- reorder(factor(gsub("^... |, .+$", "", 
                                     chess$openingDetailed)),
                         chess$moves, FUN = median)
chess$draw <- ifelse(chess$result == "draw", "draw", "result")

chess.tab <- xtabs( ~ moves + result, chess)
chess.tab.df <- as.data.frame(chess.tab)
chess.tab.df$nmoves <- as.numeric(as.character(chess.tab.df$moves))

chess.df <- subset(chess.tab.df, Freq > 0)
# Fiddle to force y-scale to include 0
chess.df <- rbind(chess.df, 
                  data.frame(moves = NA, result = "win",
                             Freq = 0, nmoves = 1000))



###################################################
### code chunk number 3: chessplot (eval = FALSE)
###################################################
## xyplot(Freq ~ nmoves | result, data = chess.df, type = "h", 
##        layout = c(1, 3), xlim = c(0, 100))   
## 


###################################################
### code chunk number 4: chess
###################################################
print(
# One tiny hidden fiddle to control y-axis labels
xyplot(Freq ~ nmoves | result, data = chess.df, type = "h", 
       layout = c(1, 3), xlim = c(0, 100),
       scales = list(y = list(at = seq(0, 6, 2))))   
)



###################################################
### code chunk number 5: chesspiece
###################################################
library("grImport")
PostScriptTrace("chess_game_01.fromInkscape.eps")
chessPicture <- readPicture("chess_game_01.fromInkscape.eps.xml")
pawn <- chessPicture[205:206]
# grid.newpage()
grid.picture(pawn)



###################################################
### code chunk number 6: chesspluspieceplot (eval = FALSE)
###################################################
## xyplot(Freq ~ nmoves | result, data = chess.df, type = "h", 
##        layout = c(1, 3), xlim = c(0, 100),
##        panel = function(...) { 
##            panel.xyplot(...)
##            grid.symbols(pawn, .05, .5, use.gc = FALSE,
##                         size = unit(.5, "npc"),
##                         gp = gpar(fill = switch(which.packet(),
##                                   "white", "grey", "black")))
##        })   
## 


###################################################
### code chunk number 7: chesspluspiece
###################################################
print(
# Fiddle to force y-scale to include 0
xyplot(Freq ~ nmoves | result, data = chess.df, type = "h", 
       layout = c(1, 3), xlim = c(0, 100),
       scales = list(y = list(at = seq(0, 6, 2))),
       panel = function(...) { 
           panel.xyplot(...)
           grid.symbols(pawn, .05, .5, use.gc = FALSE,
                        size = unit(.5, "npc"),
                        gp = gpar(fill = switch(which.packet(),
                                  "white", "grey", "black")))
       })   
)


###################################################
### code chunk number 8: petaltrace
###################################################
PostScriptTrace("petal.ps")



###################################################
### code chunk number 9: petal
###################################################
petal <- readPicture("petal.ps.xml")
grid.newpage()
grid.picture(petal)



###################################################
### code chunk number 10: import.Rnw:354-357
###################################################
petalps <- readLines("petal.ps")
cat(petalps,
    sep = "\n")


###################################################
### code chunk number 11: import.Rnw:387-391
###################################################
cat(gsub("\t", "    ", 
         gsub("(source)=", "\n         \\1=",
              readLines("petal.ps.xml"))),
    sep = "\n")


###################################################
### code chunk number 12: import.Rnw:409-410 (eval = FALSE)
###################################################
## PostScriptTrace("petal.ps")
## 


###################################################
### code chunk number 13: petaloutline
###################################################
pointify <- function(object, ...) {
    # Thin out the dots for a better diagram
    n <- length(object@x)
    subset <- c(1, seq(2, n, 3), n)
    gTree(children = gList(linesGrob(object@x[subset], object@y[subset], 
                                   default = "native",
                                   gp = gpar(col = "grey"), ...),
                         pointsGrob(object@x[subset], object@y[subset], 
                                    size = unit(2, "mm"), pch = 16, ...)))
}
grid.picture(petal, FUN = pointify)



###################################################
### code chunk number 14: flowertrace
###################################################
PostScriptTrace("flower.ps")



###################################################
### code chunk number 15: flower
###################################################
PSflower <- readPicture("flower.ps.xml")
grid.newpage()
grid.picture(PSflower)



###################################################
### code chunk number 16: import.Rnw:492-494
###################################################
cat(readLines("flower.ps"),
    sep = "\n")


###################################################
### code chunk number 17: import.Rnw:560-562
###################################################
cat(readLines("convert.ps"),
    sep = "\n")


###################################################
### code chunk number 18: import.Rnw:617-618 (eval = FALSE)
###################################################
## PostScriptTrace("flower.ps")
## 


###################################################
### code chunk number 19: import.Rnw:623-638
###################################################
flowerLines <- gsub("\t", "    ", 
                    gsub("(source)=", "\n         \\1=",
                         readLines("flower.ps.xml")))
flowerLines <- flowerLines[nchar(flowerLines) > 0]
moves <- grep("<move", flowerLines)
endpaths <- grep("</path", flowerLines)
startLine <- 1
for (i in 1:5) {
    cat(flowerLines[startLine:(moves[i])],
        "    ...",
        sep = "\n")
    startLine <- endpaths[i] - 1
}
cat(flowerLines[startLine:length(flowerLines)],
    sep = "\n")


###################################################
### code chunk number 20: import.Rnw:702-706
###################################################
flowerRGML <- xmlParse("flower.ps.xml") 
xpathApply(flowerRGML, "//path//rgb", 'xmlAttrs<-', 
           value = c(r = .3, g = .6, b = .8))
saveXML(flowerRGML, "blueflower.ps.xml")


###################################################
### code chunk number 21: blueflower
###################################################
blueflower <- readPicture("blueflower.ps.xml")
grid.picture(blueflower)



###################################################
### code chunk number 22: import.Rnw:748-749
###################################################
petal <- readPicture("petal.ps.xml")


###################################################
### code chunk number 23: import.Rnw:759-760
###################################################
str(petal)


###################################################
### code chunk number 24: import.Rnw:782-784
###################################################
PSflower <- readPicture("flower.ps.xml")



###################################################
### code chunk number 25: import.Rnw:785-786
###################################################
str(PSflower@summary)


###################################################
### code chunk number 26: import.Rnw:794-795
###################################################
petals <- PSflower[2:3]


###################################################
### code chunk number 27: import.Rnw:801-802
###################################################
str(petals@summary)


###################################################
### code chunk number 28: petals
###################################################
grid.picture(petals)



###################################################
### code chunk number 29: import.Rnw:854-856
###################################################
petal@summary@xscale
petal@summary@yscale


###################################################
### code chunk number 30: import.Rnw:877-879
###################################################
library("cluster")



###################################################
### code chunk number 31: flowerplotcode (eval = FALSE)
###################################################
## xyplot(V8 ~ V7, data = flower, 
##        xlab = "Height", ylab = "Distance Apart",
##        panel = function(x, y, ...) { 
##            grid.symbols(PSflower, x, y, units = "native", 
##                         size = unit(5, "mm")) 
##        })


###################################################
### code chunk number 32: flowerplot
###################################################
trellis.device("pdf", file = "import-flowerplot.pdf",
               width = 6, height = 4, color = FALSE)
print({ .Last.value <- 
xyplot(V8 ~ V7, data = flower, 
       xlab = "Height", ylab = "Distance Apart",
       panel = function(x, y, ...) { 
           grid.symbols(PSflower, x, y, units = "native", 
                        size = unit(5, "mm")) 
       })
}); rm(.Last.value)
dev.off()



###################################################
### code chunk number 33: tiger (eval = FALSE)
###################################################
## PostScriptTrace("tiger.ps")
## tiger <- readPicture("tiger.ps.xml")
## grid.picture(tiger[-1])
## 


###################################################
### code chunk number 34: import.Rnw:947-951
###################################################
png("import-tiger.png", width=900, height=900)
PostScriptTrace("tiger.ps")
tiger <- readPicture("tiger.ps.xml")
grid.picture(tiger[-1])

dev.off()



###################################################
### code chunk number 35: import.Rnw:1006-1008
###################################################
cat(gsub("%.+$", "", petalps[grep("curveto", petalps)]))



###################################################
### code chunk number 36: bezier
###################################################
# grid.newpage()
# Main diagram in 2x1 inch viewport (extra space is for labels)
pushViewport(viewport(width = unit(2, "inches")))
pushViewport(viewport(width = .9, height = .9,
                      xscale = c(-10, 10), yscale = c(10, 20)))
x <- c(-5, -10, 10, 5)
y <- c(10, 20, 20, 10)
grid.circle(x, y, default = "native", r = unit(1, "mm"), 
            gp = gpar(col = NA, fill = "grey"))
grid.segments(x[1], y[1], x[2], y[2],
              default = "native",
              gp = gpar(col = "grey"))
grid.segments(x[3], y[3], x[4], y[4],
              default = "native",
              gp = gpar(col = "grey"))
grid.text("(-5, 10)    ", x[1], y[1], default = "native", 
          just = c("right", "bottom"), gp = gpar(cex = .5))
grid.text("(-10, 20)   ", x[2], y[2], default = "native", 
          just = c("right", "top"), gp = gpar(cex = .5))
grid.text("   (10, 20)", x[3], y[3], default = "native", 
          just = c("left", "top"), gp = gpar(cex = .5))
grid.text("    (5, 10)", x[4], y[4], default = "native", 
          just = c("left", "bottom"), gp = gpar(cex = .5))
Ms <- 1/6*rbind(c(1, 4, 1, 0),
                c(-3, 0, 3, 0),
                c(3, -6, 3, 0),
                c(-1, 3, -3, 1))
Msinv <- solve(Ms)
# Bezier control matrix
Mb <- rbind(c(1, 0, 0, 0),
            c(-3, 3, 0, 0),
            c(3, -6, 3, 0),
            c(-1, 3, -3, 1))
# Get B-spline control points from Bezier control points by
# Msinv %*% Mb %*% bezier control points
xs <- Msinv %*% Mb %*% x
ys <- Msinv %*% Mb %*% y
grid.xspline(xs, ys,
             default = "native",
             shape = 1, repEnds = FALSE,
             gp = gpar(col = "black", lwd = 2))
popViewport(2)



###################################################
### code chunk number 37: flatbezier
###################################################
PostScriptTrace("petal.ps", "petalrough.xml", setflat = 3)
petalrough <- readPicture("petalrough.xml")
# grid.newpage()
xf <- petalrough@paths[[1]]@x
xf <- xf[-1]
xf <- rev(xf)[-1]
yf <- petalrough@paths[[1]]@y
yf <- yf[-1]
yf <- rev(yf)[-1]
pushViewport(viewport(width = unit(2, "inches")))
pushViewport(dataViewport(x, y))
x <- c(-50, -100, 100, 50)
y <- c(100, 200, 200, 100)
grid.circle(x, y, default = "native", r = unit(1, "mm"), 
            gp = gpar(col = NA, fill = "grey"))
grid.segments(x[1], y[1], x[2], y[2],
              default = "native",
              gp = gpar(col = "grey"))
grid.segments(x[3], y[3], x[4], y[4],
              default = "native",
              gp = gpar(col = "grey"))
grid.text("(-5, 10)    ", x[1], y[1], default = "native", 
          just = c("right", "bottom"), gp = gpar(cex = .5))
grid.text("(-10, 20)   ", x[2], y[2], default = "native", 
          just = c("right", "top"), gp = gpar(cex = .5))
grid.text("   (10, 20)", x[3], y[3], default = "native", 
          just = c("left", "top"), gp = gpar(cex = .5))
grid.text("    (5, 10)", x[4], y[4], default = "native", 
          just = c("left", "bottom"), gp = gpar(cex = .5))
grid.lines(xf, yf, 
           default = "native")
grid.circle(xf, yf, r = unit(.5, "mm"), 
            default = "native",
            gp = gpar(fill = "black"))
popViewport(2)



###################################################
### code chunk number 38: import.Rnw:1178-1181
###################################################
hellops <- readLines("hello.ps")
cat(hellops,
    sep = "\n")


###################################################
### code chunk number 39: hellotrace
###################################################
PostScriptTrace("hello.ps", "hello.xml")


###################################################
### code chunk number 40: hello
###################################################
hello <- readPicture("hello.xml")
grid.picture(hello)


###################################################
### code chunk number 41: hellotexttrace
###################################################
PostScriptTrace("hello.ps", 
                "hellotext.xml",
                charpath = FALSE)


###################################################
### code chunk number 42: hellotext
###################################################
hellotext <- readPicture("hellotext.xml")
grid.picture(hellotext)


###################################################
### code chunk number 43: floweroutline
###################################################
grid.picture(PSflower, 
             use.gc = FALSE, 
             gp = gpar(fill = NA, col = "black"))


###################################################
### code chunk number 44: flowerplotcode2 (eval = FALSE)
###################################################
## xyplot(V8 ~ V7, data = flower, 
##        xlab = "Height", ylab = "Distance Apart",
##        panel=function(x, y, ...) { 
##            grid.symbols(PSflower, x, y, units = "native", 
##                         size = unit(5, "mm"),
##                         use.gc = FALSE, 
##                         gp = gpar(col = "white", fill = "black", lwd = .5)) 
##        })


###################################################
### code chunk number 45: flowerplot2
###################################################
trellis.device("pdf", file = "import-flowerplot2.pdf",
               width = 6, height = 4, color = FALSE)
print({ .Last.value <- 
xyplot(V8 ~ V7, data = flower, 
       xlab = "Height", ylab = "Distance Apart",
       panel=function(x, y, ...) { 
           grid.symbols(PSflower, x, y, units = "native", 
                        size = unit(5, "mm"),
                        use.gc = FALSE, 
                        gp = gpar(col = "white", fill = "black", lwd = .5)) 
       })
}); rm(.Last.value)
dev.off()



###################################################
### code chunk number 46: import.Rnw:1397-1407
###################################################
petalrgml <- readLines("petal.ps.xml")
petalrgml <- petalrgml[nchar(petalrgml) > 0]
petalrgml <- gsub("\t", "    ", 
                  gsub("(source)=", "\n         \\1=",
                       petalrgml))
cat(petalrgml[1:(grep("<move", petalrgml) + 1)], 
    "    ...",
    petalrgml[(grep("</path", petalrgml) - 1):length(petalrgml)],
    sep = "\n")



###################################################
### code chunk number 47: import.Rnw:1432-1439
###################################################
hellotextrgml <- readLines("hellotext.xml")
hellotextrgml <- hellotextrgml[nchar(hellotextrgml) > 0]
hellotextrgml <- gsub("\t", "    ", 
                      gsub("source=", "\n         source=",
                           hellotextrgml))
cat(hellotextrgml, sep = "\n")



###################################################
### code chunk number 48: import.Rnw:1467-1468
###################################################
slotNames(petal)


###################################################
### code chunk number 49: import.Rnw:1495-1496
###################################################
str(petal@paths[[1]])


###################################################
### code chunk number 50: import.Rnw:1505-1506
###################################################
str(hellotext@paths[[1]])


###################################################
### code chunk number 51: import.Rnw:1514-1515
###################################################
str(petal@summary)


###################################################
### code chunk number 52: import.Rnw:1560-1562
###################################################
cat(readLines("blueshade.R"), sep = "\n")



###################################################
### code chunk number 53: import.Rnw:1566-1567
###################################################
source("blueshade.R")


###################################################
### code chunk number 54: import.Rnw:1576-1578
###################################################
cat(readLines("blueify.R"), sep = "\n")



###################################################
### code chunk number 55: import.Rnw:1582-1583
###################################################
source("blueify.R")


###################################################
### code chunk number 56: bluetiger (eval = FALSE)
###################################################
## grid.picture(tiger[-1], 
##              FUN = blueify)
## 


###################################################
### code chunk number 57: import.Rnw:1606-1610
###################################################
png("import-bluetiger.png", width=900, height=900)
grid.picture(tiger[-1], 
             FUN = blueify)

dev.off()



###################################################
### code chunk number 58: gnulogo
###################################################
PostScriptTrace("GNU.ps", "GNU.xml")
GNU <- readPicture("GNU.xml")
grid.picture(GNU)


###################################################
### code chunk number 59: gnulogopaths
###################################################
picturePaths(GNU, nr = 1, nc = 2, label = FALSE)


###################################################
### code chunk number 60: brokengnupaths
###################################################
brokenGNU <- explodePaths(GNU)
picturePaths(brokenGNU, nr = 3, nc = 5, 
             label = FALSE, freeScales = TRUE)


###################################################
### code chunk number 61: logobody (eval = FALSE)
###################################################
## barchart(~ cit, main = "Number of Citations per Year", xlab = "",
##          panel = function(...) {
##              grid.picture(GNU)
##              grid.rect(gp = gpar(fill = rgb(1, 1, 1, .9)))
##              panel.barchart(...)
##          })


###################################################
### code chunk number 62: logo
###################################################
cit <- c("1998"=4, "1999"=15, "2000"=17, "2001"=39, 
         "2002"=119, "2003"=276, "2004"=523, 
         "2005"=945, "2006"=1475, "2007"=2015) 
trellis.device("pdf", file = "import-logo.pdf", height = 4, color = TRUE)
print({ .Last.value <- 
barchart(~ cit, main = "Number of Citations per Year", xlab = "",
         panel = function(...) {
             grid.picture(GNU)
             grid.rect(gp = gpar(fill = rgb(1, 1, 1, .9)))
             panel.barchart(...)
         })
}); rm(.Last.value)
dev.off()



###################################################
### code chunk number 63: import.Rnw:1807-1808
###################################################
PostScriptTrace("page27.ps")


###################################################
### code chunk number 64: import.Rnw:1818-1820
###################################################
page27 <- readPicture("page27.ps.xml")
survivalPlot <- page27[c(3:16, 18, 27)]


###################################################
### code chunk number 65: survivalplot
###################################################
pushViewport(viewport(gp = gpar(lex = .2)))
grid.picture(survivalPlot)
popViewport()


###################################################
### code chunk number 66: import.Rnw:1851-1853
###################################################
zeroY <- survivalPlot@paths[[9]]@y[1]
zeroY


###################################################
### code chunk number 67: import.Rnw:1859-1861
###################################################
unitY <- (survivalPlot@paths[[14]]@y[1] - zeroY)/100
unitY


###################################################
### code chunk number 68: import.Rnw:1869-1871
###################################################
greenY <- (survivalPlot@paths[[15]]@y - zeroY)/unitY
head(round(unname(greenY), 1), n = 20)


###################################################
### code chunk number 69: import.Rnw:1876-1880
###################################################
library("survival")
sfit <- survfit(Surv(time, status) ~ trt, data = veteran)
originalGreenY <- sfit$surv[1:sfit$strata[1]]
head(round(originalGreenY*100, 1), n = 9)


