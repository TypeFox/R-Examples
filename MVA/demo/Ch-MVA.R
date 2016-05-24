### R code from vignette source 'Ch-MVA.Rnw'

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
### code chunk number 2: ch:MVA:tab:hypo
###################################################
hypo <-
  structure(list(individual = 1:10, sex = structure(c(2L, 2L, 2L,
    2L, 2L, 1L, 1L, 1L, 1L, 1L), .Label = c("Female", "Male"), class = "factor"),
    age = c(21L, 43L, 22L, 86L, 60L, 16L, NA, 43L, 22L, 80L),
    IQ = c(120L, NA, 135L, 150L, 92L, 130L, 150L, NA, 84L, 70L
    ), depression = structure(c(2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L,
    1L, 1L), .Label = c("No", "Yes"), class = "factor"), health = structure(c(3L,
    3L, 1L, 4L, 2L, 2L, 3L, 1L, 1L, 2L), .Label = c("Average",
    "Good", "Very good", "Very poor"), class = "factor"), weight = c(150L,
    160L, 135L, 140L, 110L, 110L, 120L, 120L, 105L, 100L)), .Names = c("individual",
    "sex", "age", "IQ", "depression", "health", "weight"), class = "data.frame", row.names = c(NA, -10L))
toLatex(HSAURtable(hypo), pcol = 1,
    caption = "Hypothetical Set of Multivariate Data.",
    label = "ch:MVA:tab:hypo", rownames = FALSE)


###################################################
### code chunk number 3: ch:MVA:hypo:subset
###################################################
hypo[1:2, c("health", "weight")]


###################################################
### code chunk number 4: ch:MVA:tab:measure
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

toLatex(HSAURtable(measure), pcol = 2,
    caption = "Chest, waist, and hip measurements on 20 individuals (in inches).",
    label = "ch:MVA:tab:measure", rownames = FALSE)


###################################################
### code chunk number 5: ch:MVA:tab:pottery
###################################################
data("pottery", package = "HSAUR2")
toLatex(HSAURtable(pottery), pcol = 1,
    caption = "Romano-British pottery data.",
    label = "ch:MVA:tab:pottery", rownames = FALSE)


###################################################
### code chunk number 6: ch:MVA:tab:exam
###################################################
exam <-
  structure(list(subject = 1:5, math = c(60L, 80L, 53L, 85L, 45L
  ), english = c(70L, 65L, 60L, 79L, 80L), history = c(75L, 66L,
  50L, 71L, 80L), geography = c(58L, 75L, 48L, 77L, 84L), chemistry = c(53L,
  70L, 45L, 68L, 44L), physics = c(42L, 76L, 43L, 79L, 46L)), .Names = c("subject",
  "maths", "english", "history", "geography", "chemistry", "physics"
  ), class = "data.frame", row.names = c(NA, -5L))
toLatex(HSAURtable(exam), pcol = 1,
    caption = "Exam scores for five psychology students.",
    label = "ch:MVA:tab:exam", rownames = FALSE)


###################################################
### code chunk number 7: ch:MVA:USairpollution:tab
###################################################
data("USairpollution", package = "HSAUR2")
toLatex(HSAURtable(USairpollution), pcol = 1,
    caption = paste("Air pollution in $41$ US cities."),
    label = "ch:MVA:USairpollution:tab", rownames = TRUE)


###################################################
### code chunk number 8: ch:MVA:measure:cov
###################################################
cov(measure[, c("chest", "waist", "hips")])


###################################################
### code chunk number 9: ch:MVA:measure:cov
###################################################
cov(subset(measure, gender == "female")[, 
           c("chest", "waist", "hips")])
cov(subset(measure, gender == "male")[, 
           c("chest", "waist", "hips")])


###################################################
### code chunk number 10: ch:MVA:hypo:cor
###################################################
cor(measure[, c("chest", "waist", "hips")])


###################################################
### code chunk number 11: ch:MVA:measure:dist (eval = FALSE)
###################################################
## dist(scale(measure[, c("chest", "waist", "hips")], 
##      center = FALSE))


###################################################
### code chunk number 12: ch:MVA:measure:dist
###################################################
x <- dist(scale(measure[, c("chest", "waist", "hips")], 
     center = FALSE))
as.dist(round(as.matrix(x), 2)[1:12, 1:12])
cat("...")


###################################################
### code chunk number 13: ch:MVA:fig:dmvnorm
###################################################
library("mvtnorm")
x <- y <- seq(from = -3, to = 3, length = 50)
dat <- as.matrix(expand.grid(x, y))
d <- dmvnorm(dat, mean = c(0, 0), 
             sigma = matrix(c(1, 0.5, 0.5, 1), ncol = 2))
d <- matrix(d, ncol = length(x))
persp(x = x, y = y, z = d, xlab = "x1", ylab = "x2",
      zlab = "f(x)")


###################################################
### code chunk number 14: ch:MVA:fig:cdf
###################################################
x <- seq(from = -3, to = 3, length = 1000)
Fx <- pnorm(x)
Fy <- pnorm(x, mean = -1)
plot(x, Fx, type = "l", axes = FALSE, xlab = "",
     ylab = "Cumulative distribution function") 
lines(x, Fy, type = "l")
x0 <- which.min((x - 1.2)^2)
x05 <- which.min((x + 0.5)^2)
x08 <- which.min((x + 0.9)^2)
xx <- which.min(abs(Fy - Fx[x0]))
arrows(x0 = c(min(x), x[x0], x[xx], x[x08], x[x08], x[x08]),
       y0 = c(Fx[x0], Fx[x0], Fy[xx], 0, Fx[x08], Fy[x08]), 
       x1 = c(x[x0], x[x0], x[xx], x[x08], -3, -3), 
       y1 = c(Fx[x0], 0, 0, Fy[x08], Fx[x08], Fy[x08]))
mtext(at = c(x[x08], x[xx], x[x0]), side = 1, line = 1, text =
      c(expression(q), expression(q[2](p)), expression(q[1](p))))
mtext(at = c(0, Fx[x08], Fy[x08], Fx[x0], 1), line = 1, side = 2, text =
      c(0, expression(p[1](q)), expression(p[2](q)), expression(p), 1)) 
box()


###################################################
### code chunk number 15: ch:MVA:fig:measure:chisq:setup1
###################################################
x <- measure[, c("chest", "waist", "hips")]


###################################################
### code chunk number 16: ch:MVA:fig:measure:chisq:setup2
###################################################
cm <- colMeans(x)
S <- cov(x)


###################################################
### code chunk number 17: ch:MVA:fig:measure:chisq:setup3
###################################################
d <- apply(x, MARGIN = 1, function(x) 
           t(x - cm) %*% solve(S) %*% (x - cm))


###################################################
### code chunk number 18: ch:MVA:fig:measure:qq
###################################################
qqnorm(measure[,"chest"], main = "chest"); qqline(measure[,"chest"])
qqnorm(measure[,"waist"], main = "waist"); qqline(measure[,"waist"])
qqnorm(measure[,"hips"], main = "hips"); qqline(measure[,"hips"])


###################################################
### code chunk number 19: ch:MVA:fig:measure:chisq
###################################################
plot(qchisq((1:nrow(x) - 1/2) / nrow(x), df = 3), sort(d),
     xlab = expression(paste(chi[3]^2, " Quantile")), 
     ylab = "Ordered distances")
abline(a = 0, b = 1)


###################################################
### code chunk number 20: ch:MVA:fig:USairpollution:qq:setup (eval = FALSE)
###################################################
## layout(matrix(1:8, nc = 2))
## sapply(colnames(USairpollution), function(x) {
##     qqnorm(USairpollution[[x]], main = x)
##     qqline(USairpollution[[x]])
## })


###################################################
### code chunk number 21: ch:MVA:fig:USairpollution:qq
###################################################
layout(matrix(1:8, nc = 2))
sapply(colnames(USairpollution), function(x) {
    qqnorm(USairpollution[[x]], main = x)
    qqline(USairpollution[[x]])
})


###################################################
### code chunk number 22: ch:MVA:fig:USairpollution:chisq
###################################################
x <- USairpollution
cm <- colMeans(x)
S <- cov(x)
d <- apply(x, 1, function(x) t(x - cm) %*% solve(S) %*% (x - cm))
plot(qc <- qchisq((1:nrow(x) - 1/2) / nrow(x), df = 7), 
     sd <- sort(d),
     xlab = expression(paste(chi[7]^2, " Quantile")), 
     ylab = "Ordered distances", xlim = range(qc) * c(1, 1.1))
oups <- which(rank(abs(qc - sd), ties = "random") > nrow(x) - 3)
text(qc[oups], sd[oups] - 1.5, names(oups))
abline(a = 0, b = 1)


###################################################
### code chunk number 23: ex
###################################################
s <- c(3.8778,
       2.8110,  2.1210,
       3.1480,  2.2669,  2.6550,
       3.5062,  2.5690,  2.8341,   3.2352)
S <- diag(4)
S[!lower.tri(S)] <- s
S <- S + t(S)
diag(S) <- diag(S) / 2
writeLines(apply(S, 1, function(x) paste(paste(formatC(x, format = "f"), collapse = "&"), "\\\\")))


###################################################
### code chunk number 24: ex
###################################################
X <- matrix(
 c(3, 4, 4, 6, 1,
   5, 1, 1, 7, 3,
   6, 2, 0, 2, 6,
   1, 1, 1, 0, 3,
   4, 7, 3, 6, 2,
   2, 2, 5, 1, 0,
   0, 4, 1, 1, 1,
   0, 6, 4, 3, 5,
   7, 6, 5, 1, 4,
   2, 1, 4, 3, 1), ncol = 5)
writeLines(apply(X, 1, function(x) paste(paste(x, collapse = "&"), "\\\\")))


