######################################################
#### ternary plot demo
#### Task: plotting data point hulls in a ternary plot
#### data provided by Manuel Dominguez-Rodrigo
######################################################

library(vcd)

## data

humans=matrix(c(18,19,17,21,7,9,8,62,70,53,69,81,73,71,20,10,30,10,12,18,19),
ncol=3)
colnames(humans)=c("young", "adult", "old")
lions=matrix(c(41,59,62,49,45,21,12,5,11,13,38,29,33,40,42), ncol=3)
colnames(lions)=c("young", "adult", "old")
site=matrix(c(9,12,15,11,70,62,69,68,21,26,16,21), ncol=3)
colnames(site)=c("young", "adult", "old")
humans=matrix(c(18,19,17,21,7,9,8,62,70,53,69,81,73,71,20,10,30,10,12,18,19),
ncol=3)

## regular ternary plot

data = rbind(humans, lions, site)
count = c(nrow(humans), nrow(lions), nrow(site))
rownames(data) = rep(c("humans", "lions", "site"), count)
cols = rep(c("red", "green", "blue"), count)

ternaryplot(data, col = cols)

## now try to draw hull

prop2xy <- function(x) {
  x <- as.matrix(x)
  x <- x / rowSums(x)
  xp <- x[,2] + x[,3] / 2
  yp <- x[,3] * sqrt(3) / 2
  cbind(x = xp, y = yp)
}

hullpoints <- function(x) {
    ind <- chull(x)
    ind <- c(ind, ind[1])
    x[ind,]
}

drawhull <- function(data, color) {
    hp <- hullpoints(prop2xy(data))
    grid.lines(hp[,"x"], hp[,"y"], gp = gpar(col = color))
}

## setup plot region without data points
ternaryplot(data, col = NA, pop = FALSE)

## grab plot viewport
downViewport("plot")

## now plot hulls
drawhull(humans, "blue")
drawhull(site, "red")
drawhull(lions, "green")

