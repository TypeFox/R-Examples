### R code from vignette source 'betareg-ext.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width = 70, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)
library("betareg")
combine <- function(x, sep, width) {
  cs <- cumsum(nchar(x))
  remaining <- if (any(cs[-1] > width)) combine(x[c(FALSE, cs[-1] > width)], sep, width)
  c(paste(x[c(TRUE, cs[-1] <= width)], collapse= sep), remaining)
}
prettyPrint <- function(x, sep = " ", linebreak = "\n\t", width = getOption("width")) {
  x <- strsplit(x, sep)[[1]]
  paste(combine(x, sep, width), collapse = paste(sep, linebreak, collapse = ""))
}
cache <- FALSE
enumerate <- function(x) paste(paste(x[-length(x)], collapse = ", "), x[length(x)], sep = " and ")
betamix_methods <-
  enumerate(paste("\\\\fct{", gsub("\\.betamix", "", as.character(methods(class = "betamix"))), "}", sep = ""))


###################################################
### code chunk number 2: betareg-ext.Rnw:794-797
###################################################
cat(prettyPrint(prompt(extraComponent, filename = NA)$usage[[2]], sep = ", ",
  linebreak = paste("\n", paste(rep(" ", 2), collapse = ""), sep= ""),
  width = 60))


###################################################
### code chunk number 3: betareg-ext.Rnw:804-811
###################################################
data("ReadingSkills", package = "betareg")
mean_accuracy <-
  format(round(with(ReadingSkills, tapply(accuracy, dyslexia, mean)), digits = 3),
         nsmall = 3)
mean_iq <-
  format(round(with(ReadingSkills, tapply(iq, dyslexia, mean)), digits = 3),
         nsmall = 3)


###################################################
### code chunk number 4: ReadingSkills
###################################################
data("ReadingSkills", package = "betareg")
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills, hessian = TRUE)
cl1 <- hcl(c(260, 0), 90, 40)
cl2 <- hcl(c(260, 0), 10, 95)
plot(accuracy ~ iq, data = ReadingSkills, col = cl2[as.numeric(dyslexia)],
  main = "Reading skills data", xlab = "IQ score", ylab = "Reading accuracy",
  pch = c(19, 17)[as.numeric(dyslexia)], cex = 1.5)
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5,
  pch = (1:2)[as.numeric(dyslexia)], col = cl1[as.numeric(dyslexia)])
nd <- data.frame(dyslexia = "no", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[1], lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = cl1[1], lty = 2, lwd = 2)
nd <- data.frame(dyslexia = "yes", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[2], lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = cl1[2], lty = 2, lwd = 2)
legend("topleft", c("control", "dyslexic", "betareg", "lm"),
  lty = c(NA, NA, 1:2), pch = c(19, 17, NA, NA), lwd = 2,
  col = c(cl2, 1, 1), bty = "n")
legend("topleft", c("control", "dyslexic", "betareg", "lm"),
  lty = c(NA, NA, 1:2), pch = c(1, 2, NA, NA),
  col = c(cl1, NA, NA), bty = "n")


###################################################
### code chunk number 5: ReadingSkills-bias
###################################################
data("ReadingSkills", package = "betareg")
rs_f <- accuracy ~ dyslexia * iq | dyslexia * iq
rs_ml <- betareg(rs_f, data = ReadingSkills, type = "ML")
rs_bc <- betareg(rs_f, data = ReadingSkills, type = "BC")
rs_br <- betareg(rs_f, data = ReadingSkills, type = "BR")


###################################################
### code chunk number 6: ReadingSkills-bias-table
###################################################
rs_list <- list(rs_ml, rs_bc, rs_br)
cf <- paste("$", format(round(sapply(rs_list, coef), digits = 3), nsmall = 3), "$\\phantom{)}", sep = "")
se <- paste("(", format(round(sapply(rs_list, function(x) sqrt(diag(vcov(x)))), digits = 3), nsmall = 3), ")", sep = "")
ll <- paste("$", format(round(sapply(rs_list, logLik), digits = 3), nsmall = 3), "$\\phantom{)}", sep = "")
cfse <- matrix(as.vector(rbind(cf, se)), ncol = 3)
cfse <- cbind(
  c("Mean", rep("", 7), "Precision", rep("", 7)),
  rep(as.vector(rbind(c("(Intercept)", "\\code{dyslexia}", "\\code{iq}", "\\code{dyslexia:iq}"), "")), 2),
  cfse[, 1:2],
  paste(cfse[,3],
    c(rep("\\\\", 7), "\\\\ \\hline", rep("\\\\", 7), "\\\\ \\hline")))
cfse <- rbind(cfse, c("Log-likelihood", "", ll[1:2], paste(ll[3], "\\\\ \\hline")))
writeLines(apply(cfse, 1, paste, collapse = " & "))


###################################################
### code chunk number 7: ReadingSkills-phi-plot
###################################################
pr_phi <- sapply(list("Maximum likelihood" = rs_ml,
                      "Bias correction" = rs_bc,
                      "Bias reduction" = rs_br), predict, type = "precision")
pairs(log(pr_phi), panel = function(x, y, ...) {
  panel.smooth(x, y, ...)
  abline(0, 1, lty = 2)
  })


###################################################
### code chunk number 8: ReadingSkills-noise
###################################################
set.seed(1071)
n <- nrow(ReadingSkills)
ReadingSkills$x1 <- rnorm(n)
ReadingSkills$x2 <- runif(n)
ReadingSkills$x3 <- factor(sample(0:1, n, replace = TRUE))


###################################################
### code chunk number 9: ReadingSkills-tree (eval = FALSE)
###################################################
## rs_tree <- betatree(accuracy ~ iq | iq, ~ dyslexia + x1 + x2 + x3,
##   data = ReadingSkills, minsplit = 10)


###################################################
### code chunk number 10: ReadingSkills-tree0
###################################################
if(cache & file.exists("betareg-ext-betatree.rda")) {
  load("betareg-ext-betatree.rda")
} else {
rs_tree <- betatree(accuracy ~ iq | iq, ~ dyslexia + x1 + x2 + x3,
  data = ReadingSkills, minsplit = 10)
if(cache) {
  save(rs_tree, file = "betareg-ext-betatree.rda")
} else {
  if(file.exists("betareg-ext-betatree.rda")) file.remove("betareg-ext-betatree.rda")
}
}


###################################################
### code chunk number 11: ReadingSkills-tree2 (eval = FALSE)
###################################################
## rs_tree <- betatree(accuracy ~ iq | iq | dyslexia + x1 + x2 + x3,
##   data = ReadingSkills, minsplit = 10)


###################################################
### code chunk number 12: ReadingSkills-tree3
###################################################
plot(rs_tree)


###################################################
### code chunk number 13: ReadingSkills-tree-plot
###################################################
plot(rs_tree)


###################################################
### code chunk number 14: ReadingSkills-tree-coef
###################################################
coef(rs_tree)


###################################################
### code chunk number 15: ReadingSkills-tree3
###################################################
rs_tree


###################################################
### code chunk number 16: ReadingSkills-tree-sctest
###################################################
sctest(rs_tree)


###################################################
### code chunk number 17: ReadingSkills-mix (eval = FALSE)
###################################################
## rs_mix <- betamix(accuracy ~ iq, data = ReadingSkills, k = 3,
##   extra_components = extraComponent(type = "uniform",
##     coef = 0.99, delta = 0.01), nstart = 10)


###################################################
### code chunk number 18: ReadingSkills-mix2
###################################################
if(cache & file.exists("betareg-ext-betamix.rda")) {
 load("betareg-ext-betamix.rda")
} else {
rs_mix <- betamix(accuracy ~ iq, data = ReadingSkills, k = 3,
  extra_components = extraComponent(type = "uniform",
    coef = 0.99, delta = 0.01), nstart = 10)
if(cache) {
  save(rs_mix, file = "betareg-ext-betamix.rda")
} else {
  if(file.exists("betareg-ext-betamix.rda")) file.remove("betareg-ext-betamix.rda")
}
}


###################################################
### code chunk number 19: ReadingSkills-mix3
###################################################
rs_mix


###################################################
### code chunk number 20: ReadingSkills-mix4
###################################################
summary(rs_mix)


###################################################
### code chunk number 21: ReadingSkills-betamix-plot1 (eval = FALSE)
###################################################
## ix <- as.numeric(ReadingSkills$dyslexia)
## col1 <- hcl(c(260, 0), 90, 40)[ix]
## col2 <- hcl(c(260, 0), 10, 95)[ix]
## plot(accuracy ~ iq, data = ReadingSkills, col = col2, pch = 19,
##   cex = 1.5, xlim = c(-2, 2), main = "Partitioned model (dyslexia observed)")
## points(accuracy ~ iq, data = ReadingSkills, cex = 1.5, pch = 1, col = col1)


###################################################
### code chunk number 22: ReadingSkills-betamix-plot3
###################################################
par(mfrow = c(1, 2))
ix <- as.numeric(ReadingSkills$dyslexia)
prob <- 2 * (posterior(rs_mix)[cbind(seq_along(ix), clusters(rs_mix))] - 0.5)
col3 <- hcl(c(0, 260, 130), 65, 45, fixup = FALSE)
col1 <- col3[clusters(rs_mix)]
col2 <- hcl(c(0, 260, 130)[clusters(rs_mix)], 65 * abs(prob)^1.5, 95 - 50 * abs(prob)^1.5, fixup = FALSE)
plot(accuracy ~ iq, data = ReadingSkills, col = col2, pch = 19, cex = 1.5,
  xlim = c(-2, 2), main = "Mixture model (dyslexia unobserved)")
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5, pch = 1, col = col1)
iq <- -30:30/10
cf <- rbind(coef(rs_mix, model = "mean", component = 1:2), c(qlogis(0.99), 0))
for(i in 1:3) lines(iq, plogis(cf[i, 1] + cf[i, 2] * iq), lwd = 2, col = col3[i])
ix <- as.numeric(ReadingSkills$dyslexia)
col1 <- hcl(c(260, 0), 90, 40)[ix]
col2 <- hcl(c(260, 0), 10, 95)[ix]
plot(accuracy ~ iq, data = ReadingSkills, col = col2, pch = 19,
  cex = 1.5, xlim = c(-2, 2), main = "Partitioned model (dyslexia observed)")
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5, pch = 1, col = col1)
cf <- coef(rs_tree, model = "mean")
col3 <- hcl(c(260, 0), 90, 40)
for(i in 1:2) lines(iq, plogis(cf[i, 1] + cf[i, 2] * iq), lwd = 2, col = col3[i])


###################################################
### code chunk number 23: ReadingSkills-mix5
###################################################
table(clusters(rs_mix), ReadingSkills$dyslexia)


###################################################
### code chunk number 24: GasolineYield-bias
###################################################
data("GasolineYield", package = "betareg")
gy <- lapply(c("ML", "BC", "BR"), function(x)
  betareg(yield ~ batch + temp, data = GasolineYield, type = x))


###################################################
### code chunk number 25: ReadingSkills-bias-table
###################################################
cf <- matrix(paste("$", format(round(sapply(gy, coef), digits = 5), nsmall = 5), "$\\phantom{)}", sep = ""), ncol = 3)
se <- matrix(gsub(" ", "",
  paste("(", format(round(sapply(gy, function(x) sqrt(diag(vcov(x)))), digits = 5), nsmall = 5), ")", sep = ""),
  fixed = TRUE), ncol = 3)
cfse <- cbind(cf[,1], se[,1], cf[,2], se[,2], cf[,3], se[,3])
cfse <- cbind(
  c(paste("$\\beta_{", 1:11, "}$", sep = ""), "$\\phi$"),
  cfse[, 1:5],
  paste(cfse[,6],
    c(rep("\\\\", 11), "\\\\ \\hline")))
writeLines(apply(cfse, 1, paste, collapse = " & "))


###################################################
### code chunk number 26: GasolineYield-phi
###################################################
sapply(gy, coef, model = "precision")


###################################################
### code chunk number 27: GasolineYield-phi
###################################################
sapply(gy, logLik)


###################################################
### code chunk number 28: GasolineYield-bias2
###################################################
data("GasolineYield", package = "betareg")
gy2 <- lapply(c("ML", "BC", "BR"), function(x)
  betareg(yield ~ batch + temp | 1, data = GasolineYield, type = x))
sapply(gy2, logLik)


###################################################
### code chunk number 29: ReadingSkills-bias-table
###################################################
cf <- matrix(paste("$", format(round(sapply(gy2, coef), digits = 5), nsmall = 5), "$\\phantom{)}", sep = ""), ncol = 3)
se <- matrix(gsub(" ", "",
  paste("(", format(round(sapply(gy2, function(x) sqrt(diag(vcov(x)))), digits = 5), nsmall = 5), ")", sep = ""),
  fixed = TRUE), ncol = 3)
cfse <- cbind(cf[,1], se[,1], cf[,2], se[,2], cf[,3], se[,3])
cfse <- cbind(
  c(paste("$\\beta_{", 1:11, "}$", sep = ""), "$\\log\\phi$"),
  cfse[, 1:5],
  paste(cfse[,6],
    c(rep("\\\\", 11), "\\\\ \\hline")))
writeLines(apply(cfse, 1, paste, collapse = " & "))


