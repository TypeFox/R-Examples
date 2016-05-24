### R code from vignette source 'MLCM/MLCM.Stex'

###################################################
### code chunk number 1: MLCM.Stex:8-27
###################################################
options(width=65, digits = 3)
pkg <- search()[2] 
while (search()[2] !=  if(.Platform$GUI == "AQUA") "tools:RGUI" else "package:stats") {
#detach(pos = match(pkg, search()))
spkg <- strsplit( pkg, ":"  )[[1]][2] 
if (packageHasNamespace(spkg, .libPaths()[1]) )
unloadNamespace(spkg ) else detach(pos = match(pkg, search()))
 pkg <- search()[2]
 }
rm(list = ls())
gc()
library(MPDiR)
library(lattice) 
ltheme <- canonical.theme("pdf", color = FALSE) ## in-built B&W theme 
ltheme$strip.background$bg <- "grey85" ## change strip bg 
lattice.options(default.theme = ltheme) ## set as default 
formals(deparse)$width.cutoff <- 50L
assignInNamespace("deparse", deparse, "base")
rm(deparse)


###################################################
### code chunk number 2: MLCM.Stex:30-31
###################################################
library(MPDiR)


###################################################
### code chunk number 3: MLCM.Stex:437-460
###################################################
RespFun <- function(x, n) x^n
x <- seq(0, 1, len = 200)
xy <- expand.grid(x = x, y = x)
xy$Diff <- with(xy, RespFun(x, 0.25) + RespFun(y, 1.75))
xp <- c(0.4, 0.6)
yp <- c((1 - xp^0.25)^(1/1.75), (1.25 - xp^0.25)^(1/1.75))
print(
contourplot(Diff ~ x * y, xy, at = c(1, 1.25),
	aspect = "iso", subscripts = TRUE,
	scale = list(x = list(at = xp, labels = expression(varphi^1, varphi*minute^1 )), 
				y = list(at = yp, 
				labels = expression(varphi^2, varphi*minute^2, varphi*minute*minute^2,
					varphi*minute*minute*minute^2) ), cex = 1.5
				),
	labels = FALSE,
	xlab = list(label = "Scale 1", cex = 1.3),
	ylab = list(label = "Scale 2", cex = 1.3),
	panel = function(x, y, z, subscripts,...){
		panel.contourplot(x, y, z, subscripts, lwd = 2, ...)
		panel.abline(v = xp, h = yp, lty = 2, lwd = 1)
		lpoints(c(xp, xp), yp, cex = 1.25, col = "black", lwd = 1.5)
		})  
)


###################################################
### code chunk number 4: MLCM.Stex:557-559
###################################################
data(BumpyGlossy, package = "MLCM")
head(BumpyGlossy)


###################################################
### code chunk number 5: MLCM.Stex:608-625
###################################################
library(MLCM)
library(plotrix)
mr <- par("mar")
op <- par(mar = mr + c(0, 0, 0, 4.5), pty = "s")
plot(BumpyGlossy, 
  xlab = expression(paste("Surface ", S[ij], 
  " Gloss Level (i) and Roughness level (j)")),
   ylab = expression(paste("Surface ", S[kl], 
  " Gloss Level (k) and Roughness level (l)")) 
)
par(xpd = NA)
gradient.rect(5.9, 1, 6.3, 5, 
	col = grey.colors(4, 0, 0.83), gradient = "y")
text(6.375, seq(1.5, 4.5, 1), 
    formatC(round(seq(0, 1, len = 4), 2), 2, format = "f"), 
    adj = 0)
par(op)


###################################################
### code chunk number 6: MLCM.Stex:705-716
###################################################
lmlcm <- function(p, d){
    n <- max(d[, 2:3]) - 1
    pcur <- c(0, p[1:n], 0, p[-(1:n)])
    rv <- pcur[d[, 2]] -  pcur[d[, 3]] +
        pcur[d[, 4] + 5] -  pcur[d[, 5] + 5]
    -sum(ifelse(d[, 1] == "1", pnorm(rv, log.p = TRUE), 
        pnorm(rv, lower.tail = FALSE, log.p = TRUE)))
    }
p <- c(seq(0.1, 1, len = 4), seq(1, 10, len = 4))
bg.opt <- optim(p, lmlcm, d = BumpyGlossy, 
    method = "BFGS", hessian = TRUE)


###################################################
### code chunk number 7: MLCM.Stex:722-724
###################################################
bg.opt$par
bg.opt$value


###################################################
### code chunk number 8: MLCM.Stex:736-737
###################################################
sqrt(diag(solve(bg.opt$hessian)))


###################################################
### code chunk number 9: MLCM.Stex:746-758
###################################################
opt.est <- c(0, bg.opt$par[1:4], 0, bg.opt$par[-(1:4)])
ci <- 1.96 * sqrt(diag(solve(bg.opt$hess)))
ci <- c(0, ci[1:4], 0, ci[-(1:4)])
plot(opt.est[1:5], type = "l", ylim = c(0, 6), 
    lty = 2, lwd = 2, pch = 2,
    xlab = "Physical Scale Values",
    ylab = "Perceptual Scale Values")
lines(opt.est[6:10], lwd = 2, type = "l", pch = 1)
segments(1:5, opt.est[1:5] + ci[1:5], 1:5, opt.est[1:5] - ci[1:5])
segments(1:5, opt.est[-(1:5)] + ci[-(1:5)], 1:5, opt.est[-(1:5)] - ci[-(1:5)])
points(opt.est[1:5], pch = 22, bg = "white")
points(opt.est[6:10], pch = 21, bg = "white")


###################################################
### code chunk number 10: MLCM.Stex:808-825
###################################################
make.wide <- function(d){
    nr <- nrow(d)
    wts <- rep(c(1, -1), each = nr)
    ix.mat <- matrix(0, ncol = max(d), nrow = nr)
    ix.mat[matrix(c(rep(1:nr, 2), as.vector(unlist(d))), 
        ncol = 2)] <- wts
    ix.mat <- t(apply(ix.mat, 1, function(x) 
        if (sum(x) == 0) x else
            rep(0, max(d))))
    ix.mat[, -1]
}
bg.lst <- lapply(seq(2, length(BumpyGlossy) - 1, 2), 
    function(x, d){ make.wide(d[, x:(x+1)]) }, 
    d = BumpyGlossy)
X <- do.call(cbind, bg.lst)
colnames(X) <- c(paste("G", 2:5, sep = ""), 
    paste("B", 2:5, sep = ""))


###################################################
### code chunk number 11: MLCM.Stex:831-832
###################################################
head(X)


###################################################
### code chunk number 12: MLCM.Stex:857-859
###################################################
bg.df <- data.frame(resp = BumpyGlossy$Resp, X)
bg.glm <- glm(resp ~ . -1, binomial(probit), bg.df)


###################################################
### code chunk number 13: MLCM.Stex:868-877
###################################################
opt.est <- c(0, bg.opt$par[1:4], 0, bg.opt$par[-(1:4)])
glm.est <- c(0, bg.glm$coef[1:4], 0, bg.glm$coef[-(1:4)])
plot(opt.est[1:5], type = "l", ylim = c(0, 6), lty = 2, lwd = 2,
    xlab = "Physical Scale Values",
    ylab = "Perceptual Scale Values",
    cex.lab = 1.5, cex.axis = 1.3)
lines(opt.est[6:10], lwd = 2)
points(glm.est[1:5], pch = 22, bg = "white", cex = 1.5)
points(glm.est[6:10], pch = 21, bg = "white", cex = 1.5)


###################################################
### code chunk number 14: MLCM.Stex:904-906
###################################################
bg.Ind <- update(bg.glm, . ~ . - G2 - G3 - G4 - G5)
anova(bg.Ind, bg.glm, test = "Chisq")


###################################################
### code chunk number 15: MLCM.Stex:931-951
###################################################
make.wide.full <- function(d){
  cn <- function(y, mx) (y[, 2] - 1) * mxd[2] + y[, 1]
  nr <- nrow(d)
  mxd <- c(max(d[, 1:2]), max(d[, 3:4]))
  nc <- prod(mxd)	
  nms <- sapply(seq(1, ncol(d), 2), function(x) 
     substring(names(d)[x], 1, nchar((names(d)[x])) - 1))
  fnm <- mapply(paste, nms, 
      list(seq_len(mxd[1]), seq_len(mxd[2])),
      sep = "", SIMPLIFY = FALSE)
	  nms.f <- interaction(do.call(expand.grid, fnm), 
	  sep = ":")
  ix.mat <- matrix(0, ncol = nc, nrow = nr)
  ix.mat[cbind(seq_len(nr), cn(d[, c(1, 3)], mxd))] <- 1
  ix.mat[cbind(seq_len(nr), cn(d[, c(2, 4)], mxd))] <- -1
  ix.mat <- t(apply(ix.mat, 1, function(x) if (sum(x) == 0) 
       x  else rep(0, nc)))
  colnames(ix.mat) <- levels(nms.f)
  ix.mat[, -1]
}


###################################################
### code chunk number 16: MLCM.Stex:960-964
###################################################
Xf <- make.wide.full(BumpyGlossy[, -1])
bg.dff <- data.frame(resp = BumpyGlossy[, 1], Xf)
bg.dff[1, ]
 BumpyGlossy[1, ]


###################################################
### code chunk number 17: MLCM.Stex:973-975
###################################################
bg.saturated <- update(bg.glm, data = bg.dff)
anova(bg.glm, bg.saturated, test = "Chisq")


###################################################
### code chunk number 18: MLCM.Stex:1001-1002
###################################################
rm(make.wide, make.wide.full)


###################################################
### code chunk number 19: MLCM.Stex:1004-1006
###################################################
library(MLCM)
bg.add <- mlcm(BumpyGlossy)


###################################################
### code chunk number 20: MLCM.Stex:1011-1012
###################################################
methods(class = "mlcm")


###################################################
### code chunk number 21: MLCM.Stex:1015-1017
###################################################
bg.ind <- mlcm(BumpyGlossy, model = "ind", whichdim = 2)
bg.saturated <- mlcm(BumpyGlossy, model = "full")


###################################################
### code chunk number 22: MLCM.Stex:1024-1025
###################################################
bg.add$pscale


###################################################
### code chunk number 23: MLCM.Stex:1046-1058
###################################################
opar <- par(mfrow = c(1,2), cex.lab = 1.5, cex.axis = 1.3)
plot(bg.saturated, type = "b", lty = 1, lwd = 2,
	col = c("black", paste("grey", 10 + 15 * (1:4))),
	xlab = "Glossiness Level",
	ylab = "Roughness Estimates")
mtext("a", 3, adj = 0, cex = 2, line = 1)
plot(bg.saturated, transpose = TRUE, type = "b", lty = 1, lwd = 2,
	col = c("black", paste("grey", 10 + 15 * (1:4))),
	xlab = "Roughness Level",
	ylab = "Glossiness Estimates")
mtext("b", 3, adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 24: MLCM.Stex:1095-1127
###################################################
mr <- par("mar")
opar <- par(mfrow = c(1, 2), mar = mr + c(0, 0, 0, 3))
bg.add.pred <- with(bg.add, outer(pscale[, 1], pscale[, 2], "+"))
plot(bg.saturated, standard.scale = TRUE,
	type = "b", lty = 1, lwd = 2,
	col = c("black", paste("grey", 10 + 15 * (1:4))),
	ylim = c(-0.2, 1.2),
	xlab = "Gloss Level",
	ylab = "Roughness model Estimates")
cf <- coef(lm(as.vector(bg.saturated$pscale/bg.saturated$pscale[5, 5]) ~ 
	as.vector(bg.add.pred) + 0))
matplot(cf * bg.add.pred, type = "l", lty = 2, lwd = 2,
 	col = c("black", paste("grey", 10 + 15 * (1:4))),
	add = TRUE)
bg.saturated.sc <- bg.saturated$pscale/bg.saturated$pscale[5, 5]
bg.add.adj <- cf * bg.add.pred
bg.res <- (bg.add.adj - bg.saturated.sc) + 0.5
mtext("a", 3, adj = 0, cex = 2, line = 1)
par(xpd = NA)
gcl <- grey.colors(100, min(bg.res), max(bg.res))
image(1:5, 1:5, bg.res, 
        col = gcl,
        xlab = "Gloss Level", ylab = "Roughness Level" 
        )
dgcl <- grey.colors(100, 0, 1)
gradient.rect(5.8, 1, 6.2, 5, col = dgcl,
	gradient = "y")
ylab <- -1:1
ypos <- seq(1, 5, len = 3)
text(6.7, ypos, formatC(ylab, 1, format = "f"), adj = 1)
mtext("b", 3, adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 25: MLCM.Stex:1196-1199
###################################################
bg.frm <- mlcm(~ p[1] * (x - 1)^p[2] + p[3] * (y - 1)^p[4], 
	p = c(0.1, 1, 1.5, 0.8),
	data = BumpyGlossy)


###################################################
### code chunk number 26: MLCM.Stex:1208-1209
###################################################
bg.frm$par


###################################################
### code chunk number 27: MLCM.Stex:1219-1224 (eval = FALSE)
###################################################
## xx <- seq(1, 5, len = 100)
## plot(bg.add, xlab = "Physical Scale Values",
##   ylab = "Perceptual Scale Values")
## lines(xx, predict(bg.frm, newdata = xx)[seq_along(xx)])
## lines(xx, predict(bg.frm, newdata = xx)[-seq_along(xx)])


###################################################
### code chunk number 28: MLCM.Stex:1241-1246
###################################################
xx <- seq(1, 5, len = 100)
plot(bg.add, xlab = "Physical Scale Values",
  ylab = "Perceptual Scale Values")
lines(xx, predict(bg.frm, newdata = xx)[seq_along(xx)])
lines(xx, predict(bg.frm, newdata = xx)[-seq_along(xx)])


###################################################
### code chunk number 29: MLCM.Stex:1263-1264
###################################################
anova(bg.add, bg.frm)


###################################################
### code chunk number 30: MLCM.Stex:1303-1304 (eval = FALSE)
###################################################
## bg.add.boot <- boot.mlcm(bg.add, nsim = 10000)


###################################################
### code chunk number 31: MLCM.Stex:1360-1381
###################################################
dd <- structure(list(N = c(30, 150, 300, 63, 315, 630, 201.6, 1008, 
2016, 495, 2475, 4950), SD = c(0.198804903944735, 0.0494169735922849, 
0.0483475961477474, 0.11144484767624, 0.057941390390986, 0.0498386086620115, 
0.0687411704173731, 0.0489132254288009, 0.035810703580535, 0.0710439540913592, 
0.0282549870235406, 0.0152091769143304), ScaleLevels = c(5, 5, 
5, 6, 6, 6, 8, 8, 8, 10, 10, 10), PerTrial = c(10, 50, 100, 10, 
50, 100, 10, 50, 100, 10, 50, 100)), .Names = c("N", "SD", "ScaleLevels", 
"PerTrial"), row.names = c(NA, -12L), class = "data.frame") #read.table("MLCM/SD_sampling_pwr4.txt", TRUE)
dd.lm <- lm(log10(SD) ~  offset(-0.5 * log10(N)), dd)
plot(SD ~ N, dd, xlim = c(10, 5000), ylim = c(0.005, 0.3), log = "xy",
	pch = 20 + unclass(factor(dd$ScaleLevels)), cex = 1.3,
	bg = c("white", "grey", "black")[unclass(factor(dd$PerTrial))])
x <- 10^seq(1, 3.7, len = 100)
lines(x, 10^coef(dd.lm) / sqrt(x), lwd = 2)
points(SD ~ N, dd, 
	pch = 20 + unclass(factor(dd$ScaleLevels)), cex = 1.3,
	bg = c("white", "grey", "black")[unclass(factor(dd$PerTrial))])
legend(10, 0.021, unique(dd$ScaleLevels^2), pch = 21:25, bty = "n", 
	title = "Number of Stimuli" )
legend(100, 0.021, unique(dd$PerTrial), pch = 21, bty = "n", 
	title = "%subsampling", pt.bg =  c("white", "grey", "black") )


