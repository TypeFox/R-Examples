### R code from vignette source 'Ch-PCA.Rnw'

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
### code chunk number 2: ch:PCA:data
###################################################
bc <- c(
 0.290,           
 0.202,  0.415,       
-0.055,  0.285,  0.419,       
-0.105, -0.376, -0.521, -0.877,      
-0.252, -0.349, -0.441, -0.076,  0.206,
-0.229, -0.164, -0.145,  0.023,  0.034,  0.192,
 0.058, -0.129, -0.076, -0.131,  0.151,  0.077,  0.423)
blood_sd <- c(rblood = 0.371, plate = 41.253,  wblood = 1.935,
              neut = 0.077, lymph = 0.071, bilir = 4.037,
              sodium = 2.732, potass = 0.297)
blood_corr <- diag(length(blood_sd)) / 2
blood_corr[upper.tri(blood_corr)] <- bc     
blood_corr <- blood_corr + t(blood_corr)
blood_cov <- blood_corr * outer(blood_sd, blood_sd, "*")


###################################################
### code chunk number 3: ch:PCA:blood_corr
###################################################
blood_corr


###################################################
### code chunk number 4: ch:PCA:blood_sd
###################################################
blood_sd


###################################################
### code chunk number 5: ch:PCA:blood:PCA
###################################################
blood_pcacov <- princomp(covmat = blood_cov)
summary(blood_pcacov, loadings = TRUE)
blood_pcacor <- princomp(covmat = blood_corr)
summary(blood_pcacor, loadings = TRUE)


###################################################
### code chunk number 6: ch:PCA:blood:plot1
###################################################
plot(blood_pcacor$sdev^2, xlab = "Component number",
     ylab = "Component variance", type = "l", main = "Scree diagram")
plot(log(blood_pcacor$sdev^2), xlab = "Component number",
     ylab = "log(Component variance)", type="l",
     main = "Log(eigenvalue) diagram")


###################################################
### code chunk number 7: ch:PCA:headsize:tab
###################################################
"headsize" <-
matrix(c(191, 195, 181, 183, 176, 208, 189, 197, 188, 192, 179, 183, 174, 190, 188, 163, 195, 186, 181, 175, 192, 174,
        176, 197, 190, 155, 149, 148, 153, 144, 157, 150, 159, 152, 150, 158, 147, 150, 159, 151, 137, 155, 153,
        145, 140, 154, 143, 139, 167, 163, 179, 201, 185, 188, 171, 192, 190, 189, 197, 187, 186, 174, 185, 195,
        187, 161, 183, 173, 182, 165, 185, 178, 176, 200, 187, 145, 152, 149, 149, 142, 152, 149, 152, 159, 151,
        148, 147, 152, 157, 158, 130, 158, 148, 146, 137, 152, 147, 143, 158, 150)
, nrow = 25, ncol = 4
,  dimnames = list(character(0)
, c("head1", "breadth1", "head2", "breadth2")))
x <- headsize
headsize <- as.data.frame(headsize)
toLatex(HSAURtable(headsize), pcol = 2,
    caption = "Head Size Data.",
    label = "ch:PCA:headsize:tab", rownames = FALSE)
headsize <- x


###################################################
### code chunk number 8: ch:PCA:head
###################################################
head_dat <- headsize[, c("head1", "head2")]
colMeans(head_dat)
cov(head_dat)


###################################################
### code chunk number 9: ch:PCA:head:PCA
###################################################
head_pca <- princomp(x = head_dat)
head_pca
print(summary(head_pca), loadings = TRUE)


###################################################
### code chunk number 10: ch:PCA:head:PCAint
###################################################
s1 <- round(diag(cov(head_pca$scores))[1], 3)
s2 <- round(diag(cov(head_pca$scores))[2], 3)
s <- summary(head_pca)
l1 <- round(s$loadings[,1], 2)
l2 <- round(s$loadings[,2], 2)


###################################################
### code chunk number 11: ch:PCA:head:PCAvar
###################################################
diag(cov(head_pca$scores))


###################################################
### code chunk number 12: ch:PCA:head:plot1
###################################################
a1<-183.84-0.721*185.72/0.693
b1<-0.721/0.693
a2<-183.84-(-0.693*185.72/0.721)
b2<--0.693/0.721
plot(head_dat, xlab = "First son's head length (mm)",
     ylab = "Second son's head length")
abline(a1, b1)
abline(a2, b2, lty = 2)


###################################################
### code chunk number 13: ch:PCA:head:plot2
###################################################
xlim <- range(head_pca$scores[,1])
plot(head_pca$scores, xlim = xlim, ylim = xlim)


###################################################
### code chunk number 14: ch:PCA:heptathlon:tab
###################################################
data("heptathlon",package="HSAUR2")
toLatex(HSAURtable(heptathlon), pcol = 1,
    caption = "Results of Olympic heptathlon, Seoul, 1988.",
    label = "ch:PCA:heptathlon:tab",
    rownames = TRUE)


###################################################
### code chunk number 15: ch:PCA:heptathlon:recode
###################################################
heptathlon$hurdles <- with(heptathlon, max(hurdles)-hurdles)
heptathlon$run200m <- with(heptathlon, max(run200m)-run200m)
heptathlon$run800m <- with(heptathlon, max(run800m)-run800m)
score <- which(colnames(heptathlon) == "score")
round(cor(heptathlon[,-score]), 2)
plot(heptathlon[,-score])


###################################################
### code chunk number 16: ch:PCA:heptathlon:scatter
###################################################
plot(heptathlon[,-score], pch = ".", cex = 1.5)


###################################################
### code chunk number 17: ch:PCA:heptathlon:PNG
###################################################
heptathlon <- heptathlon[-grep("PNG", rownames(heptathlon)),]
score <- which(colnames(heptathlon) == "score")
round(cor(heptathlon[,-score]), 2)


###################################################
### code chunk number 18: ch:PCA:heptathlon:scatter2
###################################################
plot(heptathlon[,-score], pch = ".", cex = 1.5)


###################################################
### code chunk number 19: PCA-opt
###################################################
op <- options(digits = 2)


###################################################
### code chunk number 20: PCA-heptathlon-pca
###################################################
heptathlon_pca <- prcomp(heptathlon[, -score], scale = TRUE)
print(heptathlon_pca)


###################################################
### code chunk number 21: PCA-heptathlon-summary
###################################################
summary(heptathlon_pca)


###################################################
### code chunk number 22: PCA-heptathlon-a1
###################################################
a1 <- heptathlon_pca$rotation[,1]
a1


###################################################
### code chunk number 23: PCA-heptathlon-scaling
###################################################
center <- heptathlon_pca$center
scale <- heptathlon_pca$scale


###################################################
### code chunk number 24: PCA-heptathlon-s1
###################################################
hm <- as.matrix(heptathlon[,-score])
drop(scale(hm, center = center, scale = scale) %*% 
     heptathlon_pca$rotation[,1])


###################################################
### code chunk number 25: PCA-heptathlon-s1
###################################################
predict(heptathlon_pca)[,1]


###################################################
### code chunk number 26: PCA-heptathlon-sdev
###################################################
sdev <- heptathlon_pca$sdev
prop12 <- round(sum(sdev[1:2]^2)/sum(sdev^2)*100, 0)


###################################################
### code chunk number 27: PCA-heptathlon-pca-plot (eval = FALSE)
###################################################
## plot(heptathlon_pca)


###################################################
### code chunk number 28: PCA-heptathlon-pca-plot
###################################################
plot(heptathlon_pca, main = "")


###################################################
### code chunk number 29: PCA-scorecor
###################################################
cor(heptathlon$score, heptathlon_pca$x[,1])


###################################################
### code chunk number 30: PCA-heptathlonscore
###################################################
plot(heptathlon$score, heptathlon_pca$x[,1])


###################################################
### code chunk number 31: ch:PCA:USairpollution:scatter
###################################################
data("USairpollution", package = "HSAUR2")
panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
}
USairpollution$negtemp <- USairpollution$temp * (-1)
USairpollution$temp <- NULL
pairs(USairpollution[,-1], diag.panel = panel.hist, 
      pch = ".", cex = 1.5)


###################################################
### code chunk number 32: ch:PCA:USairpollution:pca
###################################################
cor(USairpollution[,-1])
usair_pca <- princomp(USairpollution[,-1], cor = TRUE)


###################################################
### code chunk number 33: ch:PCA:USairpollution:pcasummary
###################################################
summary(usair_pca, loadings = TRUE)


###################################################
### code chunk number 34: ch:PCA:USairpollution:pcaplot
###################################################
pairs(usair_pca$scores[,1:3], ylim = c(-6, 4), xlim = c(-6, 4),
      panel = function(x,y, ...) {
          text(x, y, abbreviate(row.names(USairpollution)), 
               cex = 0.6)
          bvbox(cbind(x,y), add = TRUE)
      })


###################################################
### code chunk number 35: ch:PCA:USairpollution:lmplot
###################################################
out <- sapply(1:6, function(i) {
    plot(USairpollution$SO2,usair_pca$scores[,i],
         xlab = paste("PC", i, sep = ""), 
         ylab = "Sulphur dioxide concentration")
    })


###################################################
### code chunk number 36: ch:PCA:USairpollution:lm
###################################################
usair_reg <- lm(SO2 ~ usair_pca$scores, 
                data = USairpollution)
summary(usair_reg)


###################################################
### code chunk number 37: PCA-heptathlon-biplot (eval = FALSE)
###################################################
## biplot(heptathlon_pca, col = c("gray", "black"))


###################################################
### code chunk number 38: PCA-heptathlon-biplot
###################################################
tmp <- heptathlon[, -score]
rownames(tmp) <- abbreviate(gsub(" \\(.*", "", rownames(tmp)))
biplot(prcomp(tmp, scale = TRUE), col = c("black", "darkgray"), xlim =
c(-0.5, 0.7), cex = 0.7)


###################################################
### code chunk number 39: ch:PCA:headsize
###################################################
headsize.std <- sweep(headsize, 2, 
                      apply(headsize, 2, sd), FUN = "/")
R <- cor(headsize.std)
r11 <- R[1:2, 1:2]
r22 <- R[-(1:2), -(1:2)]
r12 <- R[1:2, -(1:2)]
r21 <- R[-(1:2), 1:2]
(E1 <- solve(r11) %*% r12 %*% solve(r22) %*%r21)
(E2 <- solve(r22) %*% r21 %*% solve(r11) %*%r12)
(e1 <- eigen(E1))
(e2 <- eigen(E2))


###################################################
### code chunk number 40: ch:PCA:dummy
###################################################
p <- function(x) formatC(x, format = "f", digits = 2)
f <- function(x, add = 0) paste(ifelse(x < 0, "-", "+"), p(abs(x)), "x_", 1:length(x) + add, 
                       collapse = "")
ff <- function(x, xname) paste(ifelse(x < 0, "-", "+"), p(abs(x)), "\\\\text{", xname, "}", 
                       collapse = "")


###################################################
### code chunk number 41: ch:PCA:headsize-cor
###################################################
girth1 <- headsize.std[,1:2] %*% e1$vectors[,1]
girth2 <- headsize.std[,3:4] %*% e2$vectors[,1]
shape1 <- headsize.std[,1:2] %*% e1$vectors[,2]
shape2 <- headsize.std[,3:4] %*% e2$vectors[,2]
(g <- cor(girth1, girth2))
(s <- cor(shape1, shape2))


###################################################
### code chunk number 42: ch:PCA:headsize:plot
###################################################
plot(girth1, girth2)
plot(shape1, shape2)


###################################################
### code chunk number 43: ch:PCA:LAdepr:tab
###################################################
depr <- c(
 0.212,
 0.124,  0.098,
-0.164,  0.308,  0.044,
-0.101, -0.207, -0.106, -0.208,
-0.158, -0.183, -0.180, -0.192, 0.492)
LAdepr <- diag(6) / 2
LAdepr[upper.tri(LAdepr)] <- depr
LAdepr <- LAdepr + t(LAdepr)
rownames(LAdepr) <- colnames(LAdepr) <- c("CESD", "Health", "Gender", "Age", "Edu", "Income")
x <- LAdepr
LAdepr <- as.data.frame(LAdepr)
toLatex(HSAURtable(LAdepr), 
    caption = "Los Angeles Depression Data.",
    label = "ch:PCA:LAdepr:tab", rownames = FALSE)
LAdepr <- x


###################################################
### code chunk number 44: ch:PCA:LAdepr:CCA
###################################################
r11 <- LAdepr[1:2, 1:2]
r22 <- LAdepr[-(1:2), -(1:2)]
r12 <- LAdepr[1:2, -(1:2)]
r21 <- LAdepr[-(1:2), 1:2]
(E1 <- solve(r11) %*% r12 %*% solve(r22) %*%r21)
(E2 <- solve(r22) %*% r21 %*% solve(r11) %*%r12)
(e1 <- eigen(E1))
(e2 <- eigen(E2))


