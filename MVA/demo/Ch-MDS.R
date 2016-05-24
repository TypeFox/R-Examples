### R code from vignette source 'Ch-MDS.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
library("MVA")
set.seed(280875)
library("lattice")
lattice.options(default.theme =
    function()
        standard.theme("pdf", color = FALSE))

if (file.exists("deparse.R")) {
    if (!file.exists("figs")) dir.create("figs")
    source("deparse.R")
    options(prompt = "R> ", continue = "+  ", width = 64,
        digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

    options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           figtwo =   function() {par(mfrow = c(2,1))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           figthree = function() {par(mfrow = c(3,1))},
                           fourfig =  function() {par(mfrow = c(2,2))},
                           sixfig =   function() {par(mfrow = c(3,2))},
                           nomar = function() par("mai" = c(0, 0, 0, 0))))
}


###################################################
### code chunk number 2: ch:MDS:data
###################################################
teensex <- matrix(c(21, 8, 2, 21, 9, 3, 14, 6, 4, 13, 8, 10, 8, 2, 10), nrow = 3)
rownames(teensex) <- c("No boyfriend", "Boyfriend no sex", "Boyfriend sex")
colnames(teensex) <- c("<16", "16-17", "17-18", "18-19", "19-20")
teensex <- as.table(teensex, dnn = c("Boyfriend", "Age group"))
teensex <- as.data.frame(teensex)
names(teensex) <- c("Boyfriend", "Age", "Freq")


###################################################
### code chunk number 3: ch:MDS:X
###################################################
X <- matrix(c(
    3, 5, 6, 1, 4, 2, 0, 0, 7, 2,
    4, 1, 2, 1, 7, 2, 4, 6, 6, 1,
    4, 1, 0, 1, 3, 5, 1, 4, 5, 4,
    6, 7, 2, 0, 6, 1, 1, 3, 1, 3,
    1, 3, 6, 3, 2, 0, 1, 5, 4, 1), nrow = 10)
X


###################################################
### code chunk number 4: ch:MDS:dist
###################################################
(D <- dist(X))


###################################################
### code chunk number 5: ch:MDA:cmdscale
###################################################
cmdscale(D, k = 9, eig = TRUE)


###################################################
### code chunk number 6: ch:MDA:distcmd
###################################################
max(abs(dist(X) - dist(cmdscale(D, k = 5))))


###################################################
### code chunk number 7: ch:MDA:cmdscalePCA
###################################################
max(abs(prcomp(X)$x) - abs(cmdscale(D, k = 5)))


###################################################
### code chunk number 8: ch:MDA:man
###################################################
X_m <- cmdscale(dist(X, method = "manhattan"), 
                k = nrow(X) - 1, eig = TRUE)


###################################################
### code chunk number 9: ch:MDA:man:eigen
###################################################
(X_eigen <- X_m$eig)


###################################################
### code chunk number 10: ch:MDA:man:eigen2
###################################################
cumsum(abs(X_eigen)) / sum(abs(X_eigen))
cumsum(X_eigen^2) / sum(X_eigen^2)


###################################################
### code chunk number 11: ch:MDS:airdist:tab
###################################################
"airline.dist" <-
structure(.Data = list(c(0, 587, 1212, 701, 1936, 604, 748, 2139, 2181, 543)
, c(587, 0, 920, 940, 1745, 1188, 713, 1858, 1737, 597)
, c(1212, 920, 0, 879, 831, 1726, 1631, 949, 1021, 1494)
, c(701, 940, 879, 0, 1374, 968, 1420, 1645, 1891, 1220)
, c(1936, 1745, 831, 1374, 0, 2339, 2451, 347, 959, 2300)
, c(604, 1188, 1726, 968, 2339, 0, 1092, 2594, 2734, 923)
, c(748, 713, 1631, 1420, 2451, 1092, 0, 2571, 2408, 205)
, c(2139, 1858, 949, 1645, 347, 2594, 2571, 0, 678, 2442)
, c(218, 1737, 1021, 1891, 959, 2734, 2408, 678, 0, 2329)
, c(543, 597, 1494, 1220, 2300, 923, 205, 2442, 2329, 0)
)
, names = c("ATL", "ORD", "DEN", "HOU", "LAX", "MIA",    
            "JFK", "SFO", "SEA", "IAD")
, row.names = c("ATL", "ORD", "DEN", "HOU", "LAX", "MIA",
            "JFK", "SFO", "SEA", "IAD")     
, class = "data.frame"
)

airdist <- as.dist(as.matrix(airline.dist))

tmp <- airdist
airdist <- as.data.frame(as.matrix(airdist))
htab <- HSAURtable(airdist)
htab$data[upper.tri(htab$data)] <- " "
toLatex(htab, pcol = 1, xname = "airdist",
    caption = "Airline distances between ten US cities.",
    label = "ch:MDS:airdist:tab", rownames = TRUE)
airdist <- tmp   


###################################################
### code chunk number 12: ch:MDS:airdist
###################################################
airline_mds <- cmdscale(airdist, k = 9, eig = TRUE)
airline_mds$points


###################################################
### code chunk number 13: ch:MDS:airdistprint
###################################################
airline_mds$points[, -9]


###################################################
### code chunk number 14: ch:MDS:airdist:eigen
###################################################
(lam <- airline_mds$eig)


###################################################
### code chunk number 15: ch:MDS:airdist:crit
###################################################
cumsum(abs(lam)) / sum(abs(lam))
cumsum(lam^2) / sum(lam^2)


###################################################
### code chunk number 16: ch:MDS:airdist:plot
###################################################
lim <- range(airline_mds$points[,1] * (-1)) * 1.2
plot(airline_mds$points[,1] * (-1), airline_mds$points[,2] * (-1),
     type = "n", xlab = "Coordinate 1", ylab = "Coordinate 2",
     xlim = lim, ylim = lim)
text(airline_mds$points[,1] *(-1), airline_mds$points[,2] * (-1), 
     labels(airdist), cex = 0.7)


###################################################
### code chunk number 17: ch:MDA:skulls:tab
###################################################
data("skulls", package = "HSAUR2")
toLatex(HSAURtable(skulls), pcol = 3,
    caption = "Measurements of four variables taken from Egyptian skulls of five periods.",
    label = "ch:MDA:skulls:tab")


###################################################
### code chunk number 18: ch:MDS:skulls:maha
###################################################
skulls_var <- tapply(1:nrow(skulls), skulls$epoch, 
                     function(i) var(skulls[i,-1]))
S <- 0
for (v in skulls_var) S <- S + 29 * v
(S <- S / 149)
skulls_cen <- tapply(1:nrow(skulls), skulls$epoch, 
    function(i) apply(skulls[i,-1], 2, mean))
skulls_cen <- matrix(unlist(skulls_cen), 
    nrow = length(skulls_cen), byrow = TRUE)
skulls_mah <- apply(skulls_cen, 1, 
    function(cen) mahalanobis(skulls_cen, cen, S))
skulls_mah
cmdscale(skulls_mah, k = nrow(skulls_mah) - 1, 
         eig = TRUE)$eig
skulls_mds <- cmdscale(skulls_mah)


###################################################
### code chunk number 19: ch:MDS:skulls:plot
###################################################
lim <- range(skulls_mds) * 1.2
plot(skulls_mds, xlab = "Coordinate 1", ylab = "Coordinate 2",
     xlim = lim, ylim = lim, type = "n")
text(skulls_mds, labels = levels(skulls$epoch), cex = 0.7)


###################################################
### code chunk number 20: MDS-watervoles-tab
###################################################
data("watervoles", package = "HSAUR2")
tmp <- watervoles
colnames(watervoles) <- abbreviate(colnames(watervoles))
watervoles <- as.data.frame(watervoles)
htab <- HSAURtable(watervoles)
htab$data[upper.tri(htab$data)] <- " "
toLatex(htab, pcol = 1, xname = "watervoles",
    caption = "Water voles data-dissimilarity matrix.",
    label = "MDS-watervoles-tab",
    rownames = TRUE)
watervoles <- tmp


###################################################
### code chunk number 21: MDS-voles-cmdscale
###################################################
data("watervoles", package = "HSAUR2")
voles_mds <- cmdscale(watervoles, k = 13, eig = TRUE)
voles_mds$eig


###################################################
### code chunk number 22: MDS-voles-criterion1
###################################################
cumsum(abs(voles_mds$eig))/sum(abs(voles_mds$eig))


###################################################
### code chunk number 23: MDS-voles-criterion2
###################################################
cumsum((voles_mds$eig)^2)/sum((voles_mds$eig)^2)


###################################################
### code chunk number 24: MDS-watervoles-plot
###################################################
x <- voles_mds$points[,1]
y <- voles_mds$points[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
     xlim = range(x)*1.2, type = "n")
text(x, y, labels = colnames(watervoles), cex = 0.7)


###################################################
### code chunk number 25: MDS-watervoles-mst
###################################################
library("ape")
st <- mst(watervoles)
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
     xlim = range(x)*1.2, type = "n")
for (i in 1:nrow(watervoles)) {
    w1 <- which(st[i, ] == 1)
    segments(x[i], y[i], x[w1], y[w1])
}
text(x, y, labels = colnames(watervoles), cex = 0.7)


###################################################
### code chunk number 26: MDS-voting
###################################################
library("MASS")
data("voting", package = "HSAUR2")
voting_mds <- isoMDS(voting)


###################################################
### code chunk number 27: MDS-voting-tab
###################################################
data("voting", package = "HSAUR2")
tmp <- voting
tmpname <- gsub("\\(.*", "", colnames(voting))
colnames(voting) <- abbreviate(tmpname, 3)
voting <- as.data.frame(voting)
htab <- HSAURtable(voting)
htab$data[upper.tri(htab$data)] <- " "
toLatex(htab, pcol = 1, xname = "voting",
    caption = paste("House of Representatives voting data;",
                    "(R) is short for Republican, (D) for Democrat."), 
    label = "MDS-voting-tab",
    rownames = TRUE)
voting <- tmp


###################################################
### code chunk number 28: MDS-voting-plot
###################################################
x <- voting_mds$points[,1]
y <- voting_mds$points[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
     xlim = range(voting_mds$points[,1])*1.2, type = "n")
text(x, y, labels = colnames(voting), cex = 0.6)
voting_sh <- Shepard(voting[lower.tri(voting)],
                     voting_mds$points)


###################################################
### code chunk number 29: MDS-voting-Shepard
###################################################
plot(voting_sh, pch = ".", xlab = "Dissimilarity",
     ylab = "Distance", xlim = range(voting_sh$x), 
     ylim = range(voting_sh$x))
lines(voting_sh$x, voting_sh$yf, type = "S")


###################################################
### code chunk number 30: ch:MDS:WWIIleaders:tab
###################################################
WWIIleaders <- c(
3,
4, 6, 
7, 8, 4, 
3, 5, 6, 8,
8, 9, 3, 9, 8, 
3, 2, 5, 7, 6, 7,
4, 4, 3, 5, 6, 5, 4,  
8, 9, 8, 9, 6, 9, 8, 7,       
9, 9, 5, 4, 7, 8, 8, 4, 4,
4, 5, 5, 4, 7, 2, 2, 5, 9, 5,
7, 8, 2, 4, 7, 8, 3, 2, 4, 5, 7)
tmp <- matrix(0, ncol = 12, nrow = 12)
tmp[upper.tri(tmp)] <- WWIIleaders
tmp <- tmp + t(tmp)
rownames(tmp) <- colnames(tmp) <- c("Hitler", "Mussolini", "Churchill",
    "Eisenhower", "Stalin", "Attlee", "Franco", "De Gaulle", "Mao Tse-Tung",
    "Truman", "Chamberlin", "Tito")
WWIIleaders <- as.dist(tmp)

tmp <- WWIIleaders
WWIIleaders <- as.data.frame(as.matrix(WWIIleaders))
colnames(WWIIleaders) <- abbreviate(colnames(WWIIleaders), 3)
htab <- HSAURtable(WWIIleaders)
htab$data[upper.tri(htab$data)] <- " "
toLatex(htab, pcol = 1, xname = "WWIIleaders",
    caption = "Subjective distances between WWII leaders.",
    label = "ch:MDS:WWIIleaders:tab", rownames = TRUE)
WWIIleaders <- tmp   


###################################################
### code chunk number 31: ch:MDA:WWIIleadersMDS
###################################################
(WWII_mds <- isoMDS(WWIIleaders))


###################################################
### code chunk number 32: ch:MDS:WWIIleaders:plot
###################################################
x <- WWII_mds$points[,1]
y <- WWII_mds$points[,2]
lim <- range(c(x, y)) * 1.2
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
     xlim = lim, ylim = lim, type = "n")
text(x, y, labels = labels(WWIIleaders), cex = 0.7)


###################################################
### code chunk number 33: ch:MDS:teensex:tab
###################################################
teensex <- xtabs(Freq ~ Boyfriend + Age, data = teensex)
toLatex(HSAURtable(teensex), pcol = 1,
    caption = "The influence of age on relationships with boyfriends.",
    label = "ch:MDS:teensex:tab", rownames = FALSE)


###################################################
### code chunk number 34: ch:MDS:teensex:CA
###################################################
D <- function(x) {
    a <- t(t(x) / colSums(x))
    ret <- sqrt(colSums((a[,rep(1:ncol(x), ncol(x))] - 
        a[, rep(1:ncol(x), rep(ncol(x), ncol(x)))])^2 * 
            sum(x) / rowSums(x)))
    matrix(ret, ncol = ncol(x))
}
(dcols <- D(teensex))
(drows <- D(t(teensex)))


###################################################
### code chunk number 35: ch:MDS:teensex:plot
###################################################
r1 <- cmdscale(dcols, eig = TRUE)
c1 <- cmdscale(drows, eig = TRUE)
plot(r1$points, xlim = range(r1$points[,1], c1$points[,1]) * 1.5,
     ylim = range(r1$points[,1], c1$points[,1]) * 1.5, type = "n",  
     xlab = "Coordinate 1", ylab = "Coordinate 2", lwd = 2)
text(r1$points, labels = colnames(teensex), cex = 0.7)
text(c1$points, labels = rownames(teensex), cex = 0.7)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)


###################################################
### code chunk number 36: MDS-gardenflowers-tab
###################################################
data("gardenflowers", package = "HSAUR2")
gfnames <- attr(gardenflowers, "Labels")
attr(gardenflowers, "Labels") <- gsub(" \\(.*", "", gfnames)
tmp <- as.matrix(gardenflowers)
colnames(tmp) <- abbreviate(colnames(tmp), 3)
tmp <- as.data.frame(tmp)
tmp <- HSAURtable(tmp, xname = "gardenflowers")
tmp$data[upper.tri(tmp$data)] <- " "
toLatex(tmp, pcol = 1,
    caption = "Dissimilarity matrix of $18$ species of gardenflowers.",
    label = "MDS-gardenflowers-tab",
    rownames = TRUE)


