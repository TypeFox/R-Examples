## ----label=chunk1--------------------------------------------------------
library(ade4)
library(adegraphics)
args(s.label)

## ----label=fig-gargsVSclass, include=TRUE, echo=FALSE, fig.width=7, fig.height=8----
source("gargsVSclass.R")

## ----label=chunk2--------------------------------------------------------
data(olympic)
pca1 <- dudi.pca(olympic$tab, scannf = FALSE)

## ----label=plot1, fig.height=4, fig.width=4------------------------------
g1 <- s1d.barchart(pca1$eig, p1d.horizontal = F, ppolygons.col = "white")

## ----label=chunk3--------------------------------------------------------
class(g1)
showClass("C1.barchart")

## ----label=chunk4--------------------------------------------------------
slotNames(g1)

## ----label=chunk5--------------------------------------------------------
g1@data

## ----label=plot2, fig.width=4, fig.height=4------------------------------
g2 <- s.corcircle(pca1$co)

## ----label=chunk6--------------------------------------------------------
class(g2)
g2@g.args

## ----label=plot3, fig.width=4, fig.height=4------------------------------
update(g2, fullcircle = FALSE)
g2@g.args

## ----label=chunk7--------------------------------------------------------
getcall(g1) ## equivalent to g1@Call

## ----label=plot4, fig.width=4, fig.height=4------------------------------
g3 <- s.label(pca1$li)
g4 <- s.arrow(5 * pca1$c1, add = TRUE)
class(g4)

## ----label=plot5, fig.width=4, fig.height=4------------------------------
zoom(g3, zoom = 2, center = c(2, -2))

## ----label=chunk8--------------------------------------------------------
fac.score <- factor(olympic$score < 8000, labels = c("MT8000", "LT8000"))

## ----label=plot6, fig.width=4, fig.height=4------------------------------
g5 <- s.class(pca1$li, fac.score, col = c("red", "blue"), chullSize = 1, ellipseSize = 0, 
              plabels.cex = 2, pbackground.col = "grey85", paxes.draw = TRUE)

## ----label=plot7, fig.width=4, fig.height=4------------------------------
g6 <- superpose(g5, g3, plot = TRUE) ## equivalent to g5 + g3
class(g6)

## ----label=chunk9, eval=FALSE--------------------------------------------
#  g5
#  s.label(pca1$li, add = TRUE)

## ----label=plot8, fig.width=6, fig.height=6------------------------------
rbindADEg(cbindADEg(g2, g3), cbindADEg(g5, g6), plot = TRUE)

## ----label=plot9, fig.width=4, fig.height=4------------------------------
g7 <- insert(g2, g6, posi = c(0.65, 0.65, 0.95, 0.95))
class(g7)

## ----label=chunk10-------------------------------------------------------
length(g7)
names(g7)
names(g7) <- c("chulls", "labels", "cor")
class(g7[1])
class(g7[[1]])
class(g7$chulls)

## ----label=plot10, fig.width=4, fig.height=4-----------------------------
pos.mat <- getpositions(g7)
pos.mat
pos.mat[3,] <- c(0.1, 0.7, 0.3, 0.9)
update(g7, positions = pos.mat)

## ----label=plot11, fig.width=4, fig.height=4-----------------------------
g7[[3]] <- g1
g7

## ----label=plot12, fig.width=4, fig.height=4-----------------------------
addhist(g3)

## ----label=plot13, fig.width=6, fig.height=3-----------------------------
ADEgS(adeglist = list(g2, g3), layout = c(1, 2))

## ----label=plot14, fig.width=5, fig.height=5-----------------------------
mlay <- matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 2), byrow = T, nrow = 3)
mlay
ADEgS(adeglist = list(g6, g2), layout = mlay)

## ----label=plot15, fig.width=5, fig.height=5-----------------------------
mpos <- rbind(c(0, 0.3, 0.7, 1), c(0.5, 0, 1, 0.5))
ADEgS(adeglist = list(g3, g5), positions = mpos)

## ----label=plot16, fig.width=4, fig.height=4-----------------------------
ADEgS(list(g5, g3), add = matrix(c(0, 1, 0, 0), byrow = TRUE, ncol = 2))

## ----label=plot17, fig.width=4, fig.height=4-----------------------------
data(jv73)
pca2 <- dudi.pca(jv73$morpho, scannf = FALSE)
s.label(pca2$li)

## ----label=plot18, fig.width=7, fig.height=5.2---------------------------
g8 <- s.label(pca2$li, facets = jv73$fac.riv)
length(g8)
names(g8)

## ----label=plot19, fig.width=7, fig.height=5.2---------------------------
s.class(pca2$li, fac = jv73$fac.riv, col = rainbow(12), facets = jv73$fac.riv)

## ----label=plot20, fig.width=6, fig.height=6-----------------------------
pca1 <- dudi.pca(olympic$tab, scannf = FALSE, nf = 4)
g9 <- s.corcircle(pca1$co, xax = 1:2, yax = 3:4)
length(g9)
names(g9)
g9@positions

## ----label=plot21, fig.width=6, fig.height=6-----------------------------
dim(pca1$li)
s1d.boxplot(pca1$li, fac.score, col = c("red", "blue"), 
  psub = list(position = "topleft", cex = 2))

## ----label=plot22, fig.width=6, fig.height=4-----------------------------
s.value(pca2$li, pca2$tab, symbol = "circle")

## ----label=plot23, fig.width=6, fig.height=6-----------------------------
data(meaudret)
pca3 <- dudi.pca(meaudret$env, scannf = FALSE)
pca4 <- dudi.pca(meaudret$spe, scale = FALSE, scannf = FALSE)
coi1 <- coinertia(pca3, pca4, scannf = FALSE, nf = 3)
g10 <- plot(coi1)
class(g10)
names(g10)
g10@Call

## ----label=chunk11-------------------------------------------------------
library(lattice)
sort(names(trellis.par.get()))

## ----label=plot24, fig.width=7, fig.height=3.5---------------------------
d <- scale(olympic$tab)
g11 <- table.image(d, plot = FALSE)
g12 <- table.image(d, axis.line = list(col = "blue"), axis.text = list(col = "red"), 
  plot = FALSE)
ADEgS(c(g11, g12), layout = c(1, 2))

## ----label=chunk12-------------------------------------------------------
names(adegpar())

## ----label=chunk13-------------------------------------------------------
adegpar("ppoints")
adegpar()$ppoints

## ----label=fig-paramVSparam, echo=FALSE, fig.width=7, fig.height=7-------
source("paramVSparam.R")

## ----label=fig-paramVSfunction, echo=FALSE, fig.width=7, fig.height=10----
source("paramVSfunction.R")

## ----label=plot25, fig.width=6, fig.height=3-----------------------------
oldadegpar <- adegpar()
adegpar("plabels")
g13 <- s.label(dfxy = pca1$li, plot = FALSE)
adegpar(plabels = list(col = "blue", cex = 1.5), paxes.draw = TRUE)
adegpar("plabels")
g14 <- s.label(dfxy = pca1$li, plot = FALSE)
ADEgS(c(g13, g14), layout = c(1, 2))

## ----label=chunk14-------------------------------------------------------
adegpar(oldadegpar)

## ----label=plot26, fig.width=4, fig.height=4-----------------------------
adegpar("ppoints")
s.label(dfxy = pca1$li, plabels.cex = 0, ppoints = list(col = c(2, 4, 5), cex = 1.5, pch = 15))
adegpar("ppoints")

## ----label=plot27, fig.width=7, fig.height=5.2---------------------------
s.label(pca2$li, facets = jv73$fac.riv, samelimits = FALSE)

## ----label=plot28, fig.width=6, fig.height=6-----------------------------
g15 <- plot(coi1, pbackground.col = "steelblue")

## ----label=plot29, fig.width=6, fig.height=6-----------------------------
names(g15)
plot(coi1, XYmatch.pbackground.col = "steelblue",  XYmatch.pgrid.col = "red", 
     eig.ppolygons.col = "orange")

## ----label=chunk15-------------------------------------------------------
tra1 <- list()
tra1$time <- runif(300)
tra1$distance <- tra1$time * 5 + rnorm(300)
class(tra1) <- "track"

## ----label=plot30, fig.width=7, fig.height=2.3---------------------------
g1 <- s1d.hist(tra1$distance, psub.text = "distance", ppolygons.col = "blue", 
               pgrid.draw = FALSE, plot = FALSE)
g2 <- s1d.hist(tra1$distance / tra1$time, psub.text = "speed", ppolygons.col = "red", 
               plot = FALSE)
g31 <- s.label(cbind(tra1$time, tra1$distance), paxes = list(aspectratio = "fill", 
               draw = TRUE), plot = FALSE)
g32 <- xyplot(tra1$distance ~ tra1$time, aspect = g31@adeg.par$paxes$aspectratio, 
              panel = function(x, y) {panel.lmline(x, y)})
g3 <- superpose(g31, g32)
G <- ADEgS(list(g1, g2, g3))

## ----label=chunk16-------------------------------------------------------
plot.track <- function(x, pos = -1, storeData = TRUE, plot = TRUE, ...) {
 
 ## step 1 : sort parameters for each graph
 graphsnames <- c("histDist", "histSpeed", "regression")
 sortparameters <- sortparamADEgS(..., graphsnames = graphsnames, 
                                  nbsubgraphs = c(1, 1, 2))
 
 ## step 2 : define default values for graphical parameters
 params <- list()
 params[[1]] <- list(psub = list(text = "distance"), ppolygons = list(col = "blue"), 
                     pgrid = list(draw = FALSE))
 params[[2]] <- list(psub = list(text = "speed"), ppolygons = list(col = "red"), 
                     pgrid = list(draw = FALSE))
 params[[3]] <- list()
 params[[3]]$l1 <- list(paxes = list(aspectratio = "fill", draw = TRUE))
 params[[3]]$l2 <- list()
 names(params) <- graphsnames
 sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
 
 ## step 3 : create each individual plot (ADEg)
 g1 <- do.call("s1d.hist", c(list(score = substitute(x$distance), plot = FALSE, 
               storeData = storeData, pos = pos - 2), sortparameters[[1]]))
 g2 <- do.call("s1d.hist", c(list(score = substitute(x$distance / x$time), 
               plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters[[2]]))
 g31 <- do.call("s.label", c(list(dfxy = substitute(cbind(x$time, x$distance)), plot = 
               FALSE, storeData = storeData, pos = pos - 2), sortparameters[[3]][[1]]))
 g32 <- xyplot(x$distance ~ x$time, aspect = g31@adeg.par$paxes$aspectratio,
               panel = function(x, y) {panel.lmline(x, y)})
 g3 <- do.call("superpose", list(g31, g32))
 g3@Call <- call("superpose", g31@Call, g32$call)
 
 
 ## step 4 : create the multiple plot (ADEgS)
 lay <- matrix(1:3, 1, 3)
 object <- new(Class = "ADEgS", ADEglist = list(g1, g2, g3), positions = 
               layout2position(lay), add = matrix(0, ncol = 3, nrow = 3), 
               Call = match.call())
 names(object) <- graphsnames
 if(plot)
   print(object)
 invisible(object)
}

## ----label=plot31, fig.width=7, fig.height=2.3---------------------------
plot(tra1)

## ----label=plot32, fig.width=7, fig.height=2.3---------------------------
plot(tra1, histDist.ppoly.col = "green", pbackground.col = "grey")

## ----label=plot33, fig.width=6, fig.height=3-----------------------------
data(meaudret)
g16 <- s.label(pca3$li, plot = FALSE)
g17 <- s.label(pca3$li, ppoints.col= "red", plabels = list(box = list(draw = FALSE), 
  optim = TRUE), plot = FALSE)
ADEgS(c(g16, g17), layout = c(1, 2))

## ----label=plot34, fig.width=6, fig.height=6-----------------------------
g18 <- s.class(pca3$li, fac = meaudret$design$season, plot = FALSE)
g19 <- s.class(pca3$li, fac = meaudret$design$season, ellipseSize = 0, 
  chullSize = 1, starSize = 0.5, col = TRUE, plot = FALSE)
g20 <- s.class(pca3$li, fac = meaudret$design$season, pellipses.lwd = 2, 
  pellipses.border = 2:5, pellipses.col = 2:5, plot = FALSE)
g21 <- s.class(pca3$li, fac = meaudret$design$season, ellipseSize = 0, 
  chullSize = 0, ppolygons.lwd = 2, plines.col = 2:5, starSize = 1.2, plot = FALSE)
ADEgS(c(g18, g19, g20, g21), layout = c(2, 2))

## ----label=plot35, fig.width=6, fig.height=6-----------------------------
data(rpjdl)
coa2 <- dudi.coa(rpjdl$fau, scannf = FALSE, nf = 3)
g22 <- s.value(coa2$li, coa2$li[,3], plot = FALSE)
g23 <- s.value(coa2$li, coa2$li[,3], method = "color", ppoints.cex = 0.8, 
  plegend.size= 0.8, plot = FALSE)
g24 <- s.value(coa2$li, coa2$li[,3], plegend.size = 0.8, ppoints.cex = 0.8, 
  symbol = "square", method = "color", key = list(columns = 1), 
  col = colorRampPalette(c("yellow", "blue"))(6), plot = FALSE)
g25 <- s.value(coa2$li, coa2$li[, 3], center = 0, method = "size", ppoints.cex = 0.6, 
  symbol = "circle", col = c("yellow", "red"), plot = FALSE)
ADEgS(c(g22, g23, g24, g25), layout = c(2, 2))

## ----label=plot36, fig.width=6, fig.height=6-----------------------------
score1 <- c(rnorm(1000, mean = -0.5, sd = 0.5), rnorm(1000, mean = 1))
fac1 <- rep(c("A", "B"), each = 1000)
g26 <- s1d.density(score1, fac1, pback.col = "grey75", plot = FALSE)
g27 <- s1d.density(score1, fac1, col = c(2, 4), plot = FALSE)
g28 <- s1d.density(score1, fac1, col = c(2, 4), p1d.reverse = TRUE, p1d.horizontal = FALSE, 
  p1d.rug.draw = FALSE, plot = FALSE)
g29 <- s1d.density(score1, fac1, col = c(2, 4), ppolygons.alpha = 0.2, 
  p1d = list(rug = list(tck = 1, line = FALSE)), plot = FALSE)
ADEgS(c(g26, g27, g28, g29), layout = c(2, 2))

## ----label=plot37, fig.width=6, fig.height=3-----------------------------
library(Guerry)
library(sp)
data(gfrance85)
region.names <- data.frame(gfrance85)[, 5]
col.region <- colors()[c(149, 254, 468, 552, 26)]
g30 <- s.class(coordinates(gfrance85), region.names, porigin.include = FALSE, plot = FALSE)
g31 <- s.class(coordinates(gfrance85), region.names, ellipseSize = 0, starSize = 0, 
  Sp = gfrance85, pgrid.draw = FALSE, pSp.col = col.region[region.names], pSp.alpha = 0.4, 
  plot = FALSE)
ADEgS(c(g30, g31), layout = c(1, 2))

## ----label=plot38, fig.width=6, fig.height=4-----------------------------
s.Spatial(gfrance85[,7:12])

## ----label=plot39, fig.width=6, fig.height=3-----------------------------
data(mafragh, package = "ade4")
g32 <- s.label(mafragh$xy, nb = mafragh$nb, plot = FALSE)
g33 <- s.label(mafragh$xy, nb = mafragh$nb, pnb.ed.col = "red", plab.cex = 0, 
  pnb.node = list(cex = 3, col = "blue"), ppoints.col = "green", plot = FALSE)
ADEgS(c(g32, g33), layout = c(1, 2))

## ----label=plot40, fig.width=6, fig.height=3-----------------------------
data(euro123, package = "ade4")
df <- rbind.data.frame(euro123$in78, euro123$in86, euro123$in97)
row.names(df) <- paste(row.names(euro123$in78), rep(c(1, 2, 3), rep(12, 3)), sep = "")
g34 <- triangle.label(df, label = row.names(df), showposition = TRUE, plot = FALSE)
g35 <- triangle.label(euro123$in78, plabels.cex = 0, ppoints.cex = 2, addmean = TRUE, 
  show = FALSE, plot = FALSE)
ADEgS(c(g34, g35), layout = c(1, 2))

