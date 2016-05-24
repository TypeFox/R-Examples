### R code from vignette source 'Ch-CA.Rnw'

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
### code chunk number 2: ch:CA:data
###################################################
measure <-
structure(list(V1 = 1:20, V2 = c(34L, 37L, 38L, 36L, 38L, 43L, 
40L, 38L, 40L, 41L, 36L, 36L, 34L, 33L, 36L, 37L, 34L, 36L, 38L, 
35L), V3 = c(30L, 32L, 30L, 33L, 29L, 32L, 33L, 30L, 30L, 32L, 
24L, 25L, 24L, 22L, 26L, 26L, 25L, 26L, 28L, 23L), V4 = c(32L, 
37L, 36L, 39L, 33L, 38L, 42L, 40L, 37L, 39L, 35L, 37L, 37L, 34L, 
38L, 37L, 38L, 37L, 40L, 35L)), .Names = c("V1", "V2", "V3", 
"V4"), class = "data.frame", row.names = c(NA, -20L))
measure <- measure[,-1]
names(measure) <- c("chest", "waist", "hips")
measure$gender <- gl(2, 10)
levels(measure$gender) <- c("male", "female")


###################################################
### code chunk number 3: ch:CA-scplot
###################################################
library("mvtnorm")
dat <- rbind(rmvnorm(25, mean = c(3,3)),
             rmvnorm(20, mean = c(10, 8)),
             rmvnorm(10, mean = c(20, 1)))
plot(abs(dat), xlab = expression(x[1]), ylab = expression(x[2]))


###################################################
### code chunk number 4: ch:CA-dist
###################################################
set.seed(29)
x1 <- c(0.7, 0.8, 0.85, 0.9, 1.1, 1, 0.95)
x <- c(x1, x1 + 1.5)
y1 <- sample(x1)
y <- c(y1, y1 + 1)
plot(x, y, main = "single")
lines(c(1, 0.7 + 1.5), c(1.1, 0.7 + 1), col = "grey")
set.seed(29)
x1 <- c(0.7, 0.8, 0.85, 0.9, 1.1, 1, 0.95)
x <- c(x1, x1 + 1.5)
y1 <- sample(x1)    
y <- c(y1, y1 + 1)
plot(x, y, main = "complete")
lines(c(0.7, 2.5), c(0.7, 1.1 + 1), col = "grey")
set.seed(29)
x1 <- c(0.7, 0.8, 0.85, 0.9, 1.1, 1, 0.95)
x <- c(x1, x1 + 1.5)
y1 <- sample(x1)    
y <- c(y1, y1 + 1)
plot(x, y, main = "average")
for (i in 1:7) {
    for (j in 8:14) lines(x[c(i, j)], y[c(i, j)], col = rgb(0.1, 0.1, 0.1, 0.1))
}


###################################################
### code chunk number 5: ch:CA:measure (eval = FALSE)
###################################################
## (dm <- dist(measure[, c("chest", "waist", "hips")]))


###################################################
### code chunk number 6: ch:CA:measure
###################################################
dm <- dist(measure[, c("chest", "waist", "hips")])
round(dm, 2)


###################################################
### code chunk number 7: ch:CA:measure:plot
###################################################
plot(cs <- hclust(dm, method = "single"))
plot(cc <- hclust(dm, method = "complete"))
plot(ca <- hclust(dm, method = "average"))


###################################################
### code chunk number 8: ch:CA:measure:plotplot
###################################################
body_pc <- princomp(dm, cor = TRUE)
layout(matrix(1:6, nr = 2), height = c(2, 1))
plot(cs <- hclust(dm, method = "single"), main = "Single")
abline(h = 3.8, col = "lightgrey")
xlim <- range(body_pc$scores[,1])
plot(body_pc$scores[,1:2], type = "n", xlim = xlim, ylim = xlim,
     xlab = "PC1", ylab = "PC2")
lab <- cutree(cs, h = 3.8)
text(body_pc$scores[,1:2], labels = lab, cex=0.6)
plot(cc <- hclust(dm, method = "complete"), main = "Complete")
abline(h = 10, col = "lightgrey")
plot(body_pc$scores[,1:2], type = "n", xlim = xlim, ylim = xlim,
     xlab = "PC1", ylab = "PC2")
lab <- cutree(cc, h = 10)  
text(body_pc$scores[,1:2], labels = lab, cex=0.6)     
plot(ca <- hclust(dm, method = "average"), main = "Average")
abline(h = 7.8, col = "lightgrey")
plot(body_pc$scores[,1:2], type = "n", xlim = xlim, ylim = xlim,
     xlab = "PC1", ylab = "PC2")
lab <- cutree(ca, h = 7.8)                             
text(body_pc$scores[,1:2], labels = lab, cex=0.6)     


###################################################
### code chunk number 9: ch:CA:measure:pca (eval = FALSE)
###################################################
## body_pc <- princomp(dm, cor = TRUE)
## xlim <- range(body_pc$scores[,1])
## plot(body_pc$scores[,1:2], type = "n", 
##      xlim = xlim, ylim = xlim)
## lab <- cutree(cs, h = 3.8)
## text(body_pc$scores[,1:2], labels = lab, cex = 0.6)


###################################################
### code chunk number 10: ch:CA:jet:tab
###################################################
jet <-
structure(list(V1 = c(82L, 89L, 101L, 107L, 115L, 122L, 127L,
137L, 147L, 166L, 174L, 175L, 177L, 184L, 187L, 189L, 194L, 197L,
201L, 204L, 255L, 328L), V2 = c(1.468, 1.605, 2.168, 2.054, 2.467,
1.294, 2.183, 2.426, 2.607, 4.567, 4.588, 3.618, 5.855, 2.898,
3.88, 0.455, 8.088, 6.502, 6.081, 7.105, 8.548, 6.321), V3 = c(3.3,
3.64, 4.87, 4.72, 4.11, 3.75, 3.97, 4.65, 3.84, 4.92, 3.82, 4.32,
4.53, 4.48, 5.39, 4.99, 4.5, 5.2, 5.65, 5.4, 4.2, 6.45), V4 = c(0.166,
0.154, 0.177, 0.275, 0.298, 0.15, 0, 0.117, 0.155, 0.138, 0.249,
0.143, 0.172, 0.178, 0.101, 0.008, 0.251, 0.366, 0.106, 0.089,
0.222, 0.187), V5 = c(0.1, 0.1, 2.9, 1.1, 1, 0.9, 2.4, 1.8, 2.3,
3.2, 3.5, 2.8, 2.5, 3, 3, 2.64, 2.7, 2.9, 2.9, 3.2, 2.9, 2),
    V6 = c(0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L,
    0L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L)), .Names = c("V1", "V2",
"V3", "V4", "V5", "V6"), class = "data.frame", row.names = c(NA,
-22L))                                                   
colnames(jet) <- c("FFD", "SPR", "RGF", "PLF", "SLF", "CAR")
rownames(jet) <- c("FH-1", "FJ-1", "F-86A", "F9F-2", "F-94A", "F3D-1", "F-89A",
                   "XF10F-1", "F9F-6", "F-100A", "F4D-1", "F11F-1",
                   "F-101A", "F3H-2", "F-102A", "F-8A", "F-104B",
                   "F-105B", "YF-107A", "F-106A", "F-4B", "F-111A")
jet$CAR <- factor(jet$CAR, labels = c("no", "yes"))

toLatex(HSAURtable(jet), pcol = 1,
    caption = "Jet fighters data.",
    label = "ch:CA:jet:tab")


###################################################
### code chunk number 11: ch:CA:jet:hclust
###################################################
X <- scale(jet[, c("SPR", "RGF", "PLF", "SLF")],
           center = FALSE, scale = TRUE) 
dj <- dist(X)
plot(cc <- hclust(dj), main = "Jets clustering")
cc


###################################################
### code chunk number 12: ch:CA:jet:hclust:plot
###################################################
X <- scale(jet[, c("SPR", "RGF", "PLF", "SLF")],
           center = FALSE, scale = TRUE) 
dj <- dist(X)
plot(cc <- hclust(dj), main = "Jets clustering")
cc


###################################################
### code chunk number 13: ch:CA:jet:PCA
###################################################
pr <- prcomp(dj)$x[, 1:2]
plot(pr, pch = (1:2)[cutree(cc, k = 2)],
     col = c("black", "darkgrey")[jet$CAR], 
     xlim = range(pr) * c(1, 1.5))
legend("topright", col = c("black", "black", 
                           "darkgrey", "darkgrey"), 
       legend = c("1 / no", "2 / no", "1 / yes", "2 / yes"), 
       pch = c(1:2, 1:2), title = "Cluster / CAR", bty = "n")


###################################################
### code chunk number 14: ch:CA:crime:tab
###################################################
`crime` <-
structure(c(2, 2.2, 2, 3.6, 3.5, 4.6, 10.7, 5.2, 5.5, 5.5, 6,
8.9, 11.3, 3.1, 2.5, 1.8, 9.2, 1, 4, 3.1, 4.4, 4.9, 9, 31, 7.1,
5.9, 8.1, 8.6, 11.2, 11.7, 6.7, 10.4, 10.1, 11.2, 8.1, 12.8,
8.1, 13.5, 2.9, 3.2, 5.3, 7, 11.5, 9.3, 3.2, 12.6, 5, 6.6, 11.3,
8.6, 4.8, 14.8, 21.5, 21.8, 29.7, 21.4, 23.8, 30.5, 33.2, 25.1,
38.6, 25.9, 32.4, 67.4, 20.1, 31.8, 12.5, 29.2, 11.6, 17.7, 24.6,
32.9, 56.9, 43.6, 52.4, 26.5, 18.9, 26.4, 41.3, 43.9, 52.7, 23.1,
47, 28.4, 25.8, 28.9, 40.1, 36.4, 51.6, 17.3, 20, 21.9, 42.3,
46.9, 43, 25.3, 64.9, 53.4, 51.1, 44.9, 72.7, 31, 28, 24, 22,
193, 119, 192, 514, 269, 152, 142, 90, 325, 301, 73, 102, 42,
170, 7, 16, 51, 80, 124, 304, 754, 106, 41, 88, 99, 214, 367,
83, 208, 112, 65, 80, 224, 107, 240, 20, 21, 22, 145, 130, 169,
59, 287, 135, 206, 343, 88, 106, 102, 92, 103, 331, 192, 205,
431, 265, 176, 235, 186, 434, 424, 162, 148, 179, 370, 32, 87,  
184, 252, 241, 476, 668, 167, 99, 354, 525, 319, 605, 222, 274, 
408, 172, 278, 482, 285, 354, 118, 178, 243, 329, 538, 437, 180,
354, 244, 286, 521, 401, 103, 803, 755, 949, 1071, 1294, 1198,
1221, 1071, 735, 988, 887, 1180, 1509, 783, 1004, 956, 1136,
385, 554, 748, 1188, 1042, 1296, 1728, 813, 625, 1225, 1340,
1453, 2221, 824, 1325, 1159, 1076, 1030, 1461, 1787, 2049, 783,
1003, 817, 1792, 1845, 1908, 915, 1604, 1861, 1967, 1696, 1162,
1339, 2347, 2208, 2697, 2189, 2568, 2758, 2924, 2822, 1654, 2574,
2333, 2938, 3378, 2802, 2785, 2801, 2500, 2049, 1939, 2677, 3008,
3090, 2978, 4131, 2522, 1358, 2423, 2846, 2984, 4373, 1740, 2126,
2304, 1845, 2305, 3417, 3142, 3987, 3314, 2800, 3078, 4231, 3712,
4337, 4074, 3489, 4267, 4163, 3384, 3910, 3759, 164, 228, 181,
906, 705, 447, 637, 776, 354, 376, 328, 628, 800, 254, 288, 158,
439, 120, 99, 168, 258, 272, 545, 975, 219, 169, 208, 277, 430,
598, 193, 544, 267, 150, 195, 442, 649, 714, 215, 181, 169, 486,
343, 419, 223, 478, 315, 402, 762, 604, 328), .Dim = c(51L, 7L
), .Dimnames = list(c("ME", "NH", "VT", "MA", "RI", "CT", "NY",
"NJ", "PA", "OH", "IN", "IL", "MI", "WI", "MN", "IA", "MO", "ND",
"SD", "NE", "KS", "DE", "MD", "DC", "VA", "WV", "NC", "SC", "GA",
"FL", "KY", "TN", "AL", "MS", "AR", "LA", "OK", "TX", "MT", "ID",
"WY", "CO", "NM", "AZ", "UT", "NV", "WA", "OR", "CA", "AK", "HI"
), c("Murder", "Rape", "Robbery", "Assault", "Burglary", "Theft",
"Vehicle")))
crime <- as.data.frame(crime)

toLatex(HSAURtable(crime), pcol = 1, rownames = TRUE,
    caption = "Crime data.",
    label = "ch:CA:crime:tab")


###################################################
### code chunk number 15: ch:CA:crime:plot
###################################################
plot(crime, pch = ".", cex = 1.5)


###################################################
### code chunk number 16: ch:CA:crime:outlier
###################################################
subset(crime, Murder > 15)


###################################################
### code chunk number 17: ch:CA:crime:plot2
###################################################
plot(crime, pch = c(".", "+")[(rownames(crime) == "DC") + 1], cex = 1.5)


###################################################
### code chunk number 18: ch:CA:crime:var
###################################################
sapply(crime, var)


###################################################
### code chunk number 19: ch:CA:crime:stand
###################################################
rge <- sapply(crime, function(x) diff(range(x)))
crime_s <- sweep(crime, 2, rge, FUN = "/")
sapply(crime_s, var)


###################################################
### code chunk number 20: ch:CA:crime:wss
###################################################
n <- nrow(crime_s)
wss <- rep(0, 6)
wss[1] <- (n - 1) * sum(sapply(crime_s, var))
for (i in 2:6)
    wss[i] <- sum(kmeans(crime_s,
                         centers = i)$withinss)
plot(1:6, wss, type = "b", xlab = "Number of groups",
     ylab = "Within groups sum of squares")


###################################################
### code chunk number 21: ch:CA:crimes:k2
###################################################
kmeans(crime_s, centers = 2)$centers * rge


###################################################
### code chunk number 22: ch:CA:crime:PCA
###################################################
crime_pca <- prcomp(crime_s)
plot(crime_pca$x[, 1:2], 
     pch = kmeans(crime_s, centers = 2)$cluster)


###################################################
### code chunk number 23: ca:CA:pottery:dist
###################################################
pottery_dist <- dist(pots <- scale(pottery[, colnames(pottery) != "kiln"], 
                                   center = FALSE))
library("lattice")
levelplot(as.matrix(pottery_dist), xlab = "Pot Number",
          ylab = "Pot Number")


###################################################
### code chunk number 24: ch:CA:pottery:distplot
###################################################
trellis.par.set(standard.theme(color = FALSE))
plot(levelplot(as.matrix(pottery_dist), xlab = "Pot Number", ylab = "Pot Number",
     scales = list(x = list(draw = FALSE), y = list(draw = FALSE))))


###################################################
### code chunk number 25: ch:CA:pottery:wss
###################################################
n <- nrow(pots)
wss <- rep(0, 6)
wss[1] <- (n - 1) * sum(sapply(pots, var))
for (i in 2:6)
    wss[i] <- sum(kmeans(pots,
                         centers = i)$withinss)
plot(1:6, wss, type = "b", xlab = "Number of groups",   
     ylab = "Within groups sum of squares")


###################################################
### code chunk number 26: ch:CA:pots:PCA
###################################################
pots_pca <- prcomp(pots)
plot(pots_pca$x[, 1:2], 
     pch = kmeans(pots, centers = 3)$cluster)


###################################################
### code chunk number 27: ch:CA:pottery:cluster
###################################################
set.seed(29)
pottery_cluster <- kmeans(pots, centers = 3)$cluster
xtabs(~ pottery_cluster + kiln, data = pottery)


###################################################
### code chunk number 28: ch:CA:thompson
###################################################
cnt <- c("Iceland", "Norway", "Sweden", "Finland", "Denmark", "UK", "Eire",
         "Germany", "Netherlands", "Belgium", "Switzerland", "France", "Spain",
         "Portugal", "Italy", "Greece", "Yugoslavia", "Albania", "Bulgaria", "Romania",
         "Hungary", "Czechia", "Slovakia", "Poland", "CIS", "Lithuania", "Latvia", "Estonia")

thomson <- expand.grid(answer = factor(c("no", "yes")),
                       question = factor(paste("Q", 1:6, sep = "")),
                       country = factor(cnt, levels = cnt))  
thomson$Freq <- c(
0, 5, 0, 5, 0, 4, 0, 5, 0, 5, 0, 5,
1, 6, 1, 5, 0, 6, 0, 5, 0, 4, 1, 4,
0, 11, 4, 7, 0, 7, 0, 11, 5, 5, 3, 6,
0, 6, 2, 4, 0, 6, 0, 6, 1, 5, 2, 4,
1, 12, 4, 9, 0, 12, 3, 9, 7, 4, 6, 7,
7, 12, 2, 16, 0, 20, 1, 19, 9, 10, 0, 17,
0, 1, 1, 2, 0, 3, 2, 0, 2, 0, 0, 3,
0, 14, 0, 13, 0, 13, 2, 12, 11, 2, 1, 13,
0, 8, 0, 8, 0, 8, 1, 7, 2, 5, 1, 7,
2, 0, 0, 2, 0, 2, 1, 1, 2, 0, 0, 2,
0, 5, 0, 5, 0, 4, 2, 2, 5, 0, 0, 4,
7, 3, 1, 7, 3, 5, 8, 2, 10, 0, 1, 7,
11, 1, 0, 12, 2, 8, 5, 6, 11, 0, 0, 11,
5, 1, 0, 6, 2, 4, 3, 3, 6, 0, 0, 6,
8, 7, 0, 15, 1, 13, 9, 6, 13, 2, 0, 15,
7, 1, 0, 8, 3, 5, 7, 1, 8, 0, 0, 7,
11, 4, 0, 15, 7, 8, 11, 4, 15, 0, 0, 14,
3, 2, 2, 3, 3, 2, 3, 2, 3, 3, 2, 3,
3, 0, 0, 3, 2, 1, 3, 0, 3, 0, 0, 3,
7, 0, 0, 6, 6, 1, 6, 1, 6, 1, 0, 7,
4, 1, 0, 5, 1, 4, 5, 0, 5, 0, 0, 5,
18, 2, 0, 20, 17, 3, 20, 0, 20, 0, 0, 20,
13, 0, 1, 14, 14, 0, 16, 0, 13, 0, 15, 0,
18, 0, 0, 19, 13, 5, 17, 2, 17, 0, 0, 19,
7, 0, 1, 6, 5, 2, 7, 0, 7, 0, 1, 6,
8, 0, 0, 8, 8, 0, 8, 0, 8, 0, 0, 8,
5, 0, 0, 5, 5, 0, 5, 0, 5, 0, 0, 5,
2, 2, 0, 3, 0, 3, 3, 0, 3, 0, 0, 3)
ttab <- xtabs(Freq ~ country + answer + question, data = thomson)

thomsonprop <- prop.table(ttab, c(1,3))[,"yes",]

plot(1:(22 * 6), rep(-1, 22 * 6), 
     ylim = c(-nlevels(thomson$country), -1), type = "n",
     axes = FALSE, xlab = "", ylab = "")
for (q in 1:6) {   
   tmp <- ttab[,,q]
   xstart <- (q - 1) * 22 + 1
   y <- -rep(1:nrow(tmp), rowSums(tmp))
   x <- xstart + unlist(sapply(rowSums(tmp), function(i) 1:i))
   pch <- unlist(apply(tmp, 1, function(x) c(rep(19, x[2]), rep(1, x[1]))))
   points(x, y, pch = pch)
}
axis(2, at = -(1:nlevels(thomson$country)), labels = levels(thomson$country),
     las = 2, tick = FALSE, line = 0)
mtext(text = paste("Question", 1:6), 3, at = 22 * (0:5), adj = 0)


###################################################
### code chunk number 29: ch:CA:thompsonMC
###################################################
library("mclust")
(mc <- Mclust(thomsonprop))


###################################################
### code chunk number 30: ch:CA:thompsonMC:plot
###################################################
plot(mc, thomsonprop, what = "BIC", col = "black")


###################################################
### code chunk number 31: ch:CA:thompsonMC
###################################################
cl <- mc$classification
nm <- unlist(sapply(1:3, function(i) names(cl[cl == i])))
ttab <- ttab[nm,,]
plot(1:(22 * 6), rep(-1, 22 * 6), 
     ylim = c(-nlevels(thomson$country), -1), type = "n",
     axes = FALSE, xlab = "", ylab = "")
for (q in 1:6) {   
   tmp <- ttab[,,q]
   xstart <- (q - 1) * 22 + 1
   y <- -rep(1:nrow(tmp), rowSums(tmp))
   x <- xstart + unlist(sapply(rowSums(tmp), function(i) 1:i))
   pch <- unlist(apply(tmp, 1, function(x) c(rep(19, x[2]), rep(1, x[1]))))
   points(x, y, pch = pch)
}
axis(2, at = -(1:nlevels(thomson$country)), labels = dimnames(ttab)[[1]],
     las = 2, tick = FALSE, line = 0)
mtext(text = paste("Question", 1:6), 3, at = 22 * (0:5), adj = 0)
abline(h = -cumsum(table(cl))[-3] - 0.5, col = "grey")
text(-c(0.75, 0.75, 0.75), -cumsum(table(cl)) + table(cl)/2,
     label = paste("Cluster", 1:3), srt = 90, pos = 1)


###################################################
### code chunk number 32: ch:CA:neighbor
###################################################
library("flexclust")
library("mvtnorm")
set.seed(290875)
x <- rbind(rmvnorm(n = 20, mean = c(0, 0), 
                   sigma = diag(2)),
           rmvnorm(n = 20, mean = c(3, 3), 
                   sigma = 0.5 * diag(2)),
           rmvnorm(n = 20, mean = c(7, 6), 
                   sigma = 0.5 * (diag(2) + 0.25)))
k <- cclust(x, k = 5, save.data = TRUE)
plot(k, hull = FALSE, col = rep("black", 5), xlab = "x", ylab = "y")


###################################################
### code chunk number 33: ch:CA:pottery:neighbor
###################################################
k <- cclust(pots, k = 3, save.data = TRUE)
plot(k, project = prcomp(pots), hull = FALSE, col = rep("black", 3),
     xlab = "PC1", ylab = "PC2")     


###################################################
### code chunk number 34: ch:CA:art:stripes1
###################################################
set.seed(912345654)
x <- rbind(matrix(rnorm(100, sd = 0.5), ncol= 2 ),
           matrix(rnorm(100, mean =4, sd = 0.5), ncol = 2),
           matrix(rnorm(100, mean =7, sd = 0.5), ncol = 2),
           matrix(rnorm(100, mean =-1.0, sd = 0.7), ncol = 2),
           matrix(rnorm(100, mean =-4.0, sd = 1.0), ncol = 2))
c5 <- cclust(x, 5, save.data = TRUE)
stripes(c5, type = "second", col = 1)


###################################################
### code chunk number 35: ch:CA:art:stripes2
###################################################
set.seed(912345654)
x <- rbind(matrix(rnorm(100, sd = 2.5), ncol = 2),
           matrix(rnorm(100, mean = 3, sd = 0.5), ncol = 2),
           matrix(rnorm(100, mean = 5, sd = 0.5), ncol = 2),
           matrix(rnorm(100, mean = -1.0, sd = 1.5), ncol = 2),
           matrix(rnorm(100, mean = -4.0, sd = 2.0), ncol = 2))
c5 <- cclust(x, 5, save.data = TRUE)
stripes(c5, type = "second", col = 1)


###################################################
### code chunk number 36: ch:CA:pottery:stripes
###################################################
set.seed(15)
c5 <- cclust(pots, k = 3, save.data = TRUE)
stripes(c5, type = "second", col = "black")


