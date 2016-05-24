logi.hist.plot <- function (independ, depend, logi.mod = 1,
type = "dit", boxp = TRUE, rug = FALSE,
  ylabel= "Probability", ylabel2="Frequency", xlabel="", mainlabel="",
las.h = 1, counts=FALSE, ...){

# get the label for the x-axis
#xlabel <- paste(deparse(substitute(independ)))

# define functions:

# set the draw area if no box plots are to be drawn
logi.scater <- function (independ, depend, scater = "n",
x.lab = xlabel, las =las.h){
plot(independ, depend, cex = 1, type = scater,
ylab = ylabel, xlab = x.lab,main= mainlabel,
cex.lab = 1.2, las = las)
}

# add rug plot if desired; you could change pch.rug
# (symbol type) or cex.rug (symbol size)
logi.rug <- function (independ, depend, pch.rug = 16,
cex.rug = 1){
points(independ, depend, pch = pch.rug ,cex = cex.rug)
}

# set the draw area and add box plots; you could change
# cold.box (color of the boxes)
logi.box <- function(independ, depend, col.box = "gray",
x.lab = xlabel, las = las.h){
plot(independ, depend, cex = 1, type = "n",
ylim = c(-0.1,1.1), ylab = ylabel,    # change ylabel2 to ylabel Nov 2, 2015
xlab = x.lab, cex.lab = 1.2, las = las)
indep.1 <- independ[depend == 1]
indep.0 <- independ[depend == 0]
boxplot(indep.1, horizontal = TRUE, add = TRUE,
at = 1.05, boxwex = 0.1, col = col.box, notch = TRUE)
boxplot(indep.0, horizontal = TRUE, add = TRUE,
at = -0.05, boxwex = 0.1, col = col.box, notch = TRUE)
}

# fit binomial glm and add predicted curve; you could
# change col.cur (color of the curve) or lwd.cur(width
# of the curve)
logi.curve <- function(independ, depend, mod = logi.mod,
col.cur = "red", lwd.cur = 4){
if (mod == 1) mod3 <- glm(depend ~ independ,
family = binomial)
if (mod == 2) mod3 <- glm(depend ~ independ +
I(independ^2), family = binomial)
x.new <- seq(min(independ), max(independ), len = 100)
y.new <- predict(mod3, data.frame(independ = x.new),
type = "response")
lines(x.new, y.new, lwd = lwd.cur, col = col.cur)
}

# add dit plot; you may want to change pch.dit (type of
# points), cex.p (size of points), and incre (space
# between points)
logi.dit <- function (independ, depend, cex.p = 1,
pch.dit = 1, incre = 0.02){

indep.0 <- independ[depend == 0]
indep.1 <- independ[depend == 1]
uni.plot.0 <- function(x) length(which(indep.0 == x))
uni.plot.1 <- function(x) length(which(indep.1 == x))

# get the number of repeated values of "independ":

cosa.0 <- apply(as.matrix(unique(indep.0)), 1, uni.plot.0)
cosa.1 <- apply(as.matrix(unique(indep.1)), 1, uni.plot.1)

# start ploting:
points(independ, depend, pch = pch.dit, cex = cex.p)

for (i in 1:max(cosa.0)){
for (j in 1:i){
points(unique(indep.0)[which(cosa.0 == i+1)],
rep(0 + incre*j, length(which(cosa.0 == i+1))),
pch = pch.dit, cex = cex.p)
}
}

for (i in 1:max(cosa.1)){
for (j in 1:i){
points(unique(indep.1)[which(cosa.1 == i+1)],
rep(1 - incre*j, length(which(cosa.1 == i+1))),
pch = pch.dit, cex = cex.p)
}
}
}

# add histograms and frequency axes; you may want to change
# scale.hist (factor to scale histogram height to 0-1
# interval) or col.hist (color of histogram)
logi.hist <- function(independ, depend, scale.hist = 5,
col.hist = "blue", count.hist=TRUE,
intervalo = 0, las.h1 = las.h){

# get the position of bins
h.br <- hist(independ, plot = FALSE)$br
if (intervalo > 0) h.br <- seq(from = range(h.br)[1],
to = range(h.br)[2], by = intervalo)
h.x <- hist(independ[depend == 0], breaks = h.br,
plot = FALSE)$mid

# get counts in each bin
h.y0 <- hist(independ[depend == 0], breaks = h.br,
plot = FALSE)$counts
h.y1 <- hist(independ[depend == 1], breaks = h.br,
plot = FALSE)$counts

# scale the histogram bars to max desired length:
h.y0n <- h.y0/(max(c(h.y0,h.y1))* scale.hist)
h.y1n <- 1 - h.y1/(max(c(h.y0,h.y1))* scale.hist)

# draw bottom histogram:
for (i in 1:length(h.y0n)){
if (h.y0n[i] > 0)
polygon(c(rep(h.br[i], 2), rep(h.br[i+1], 2)),
c(0, rep(h.y0n[i], 2), 0), col = col.hist)
}

# draw top histogram:
for (i in 1:length(h.y1n)){
if (h.y1n[i] < 1)
polygon(c(rep(h.br[i], 2), rep(h.br[i+1], 2)),
c(h.y1n[i], 1, 1, h.y1n[i]), col = col.hist)
}

# add counts to bins if required:
if (counts == TRUE)
for (i in 1 : length(h.x)){
text(h.x[i], h.y1n[i], h.y1[i], cex = 1, pos = 1)
text(h.x[i], h.y0n[i], h.y0[i], cex = 1, pos = 3)
}

# plot the axes of histograms:
axis.hist <- function (h.y0, h.y1, scale.hist,
las = las.h1){
tope <- max(c(h.y0, h.y1))
label.down <- c(0, (ceiling(tope/10))*5,
(ceiling(tope/10))*10)
label.up <- c((ceiling(tope/10))*10,
(ceiling(tope/10))*5, 0)
at.down <- label.down/(tope * scale.hist)
at.up <- 1 - (label.up/(tope * scale.hist))
at.hist <- c(at.down, at.up)
label.hist <- c(label.down, label.up)
axis(side = 4, at = at.hist, labels = label.hist,
las = las)
mtext(ylabel2, side = 4, line = 2, cex = 1.2)
}
axis.hist(h.y0, h.y1, scale.hist)
axis (side = 2, las = las.h1)
}


# set the margins of plot area
old.mar <- par()$mar
par(mar = c(5.1,4.1,4.1,4.1))

# plot the combined graph
if (boxp == TRUE) logi.box(independ, depend)
if (boxp == FALSE) logi.scater(independ, depend)
if (type != "dit") logi.hist(independ, depend,...)
if (rug == TRUE) logi.rug (independ, depend)
logi.curve(independ, depend)
if (type == "dit") logi.dit(independ, depend)

# reset the margins to old margins
par(mar = old.mar)
}

