## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("boot")) {
  boot.present <- FALSE
  cat("Package boot not available - some code may not run.\nInstall it by typing 'install.packages(\"boot\")'")
} else {
  boot.present <- TRUE
}
if (!require("leaps")) {
  leaps.present <- FALSE
  cat("Package leaps not available - some code may not run.\nInstall it by typing 'install.packages(\"leaps\")'")
} else {
  leaps.present <- TRUE
}
if (!require("subselect")) {
  subselect.present <- FALSE
  cat("Package subselect not available - some code may not run.\nInstall it by typing 'install.packages(\"subselect\")'")
} else {
  subselect.present <- TRUE
}
if (!require("lars")) {
  lars.present <- FALSE
  cat("Package lars not available - some code may not run.\nInstall it by typing 'install.packages(\"lars\")'")
} else {
  lars.present <- TRUE
}
if (!require("elasticnet")) {
  elasticnet.present <- FALSE
  cat("Package elasticnet not available - some code may not run.\nInstall it by typing 'install.packages(\"elasticnet\")'")
} else {
  elasticnet.present <- TRUE
}

data(wines, package = "ChemometricsWithRData")
odd <- seq(1, nrow(wines), 2)
X <- wines[odd,]                     
C <- classvec2classmat(vintages[odd])
wines.lm <- lm(C ~ X)                
wines.lm.summ <- summary(wines.lm)                    
wines.lm.summ[[3]]

sapply(wines.lm.summ, function(x) which(x$coefficients[,4] < .1))

if (boot.present) {
  data(gasoline, package = "pls")
  wavelengths <- seq(900, 1700, by = 2)
  
  set.seed(7)
  gas.pcr.bootCI <-
    boot(gasoline,
         function(x, ind) {
           c(coef(pcr(octane ~ ., data=x, subset = ind, ncomp = 4)))
         },
         R = 999)
  gas.BCACI <-
    sapply(1:ncol(gasoline$NIR),
           function(i, x) {
             boot.ci(x, index = i, type = "bca")$bca[,4:5]},
           gas.pcr.bootCI)
  coefs <- gas.pcr.bootCI$t0
  matplot(wavelengths, t(gas.BCACI), type = "n",
          xlab = "Wavelength (nm)", ylab = "Regression coefficient",
          main = "Gasoline data: PCR (4 PCs)")
  abline(h = 0, col = "gray")
  lines(wavelengths, coefs, col = "gray")
  matlines(wavelengths, t(gas.BCACI), col = 2, lty = 1)
  insignif <- apply(gas.BCACI, 2, prod) < 0
  coefs[insignif] <- NA
  lines(wavelengths, coefs, lwd = 2)
  
  odd <- seq(1, nrow(gasoline$NIR), by = 2)
  smallmod <- pcr(octane ~ NIR[,!insignif], data = gasoline,
                  subset = odd, ncomp = 4, validation = "LOO")
  RMSEP(smallmod, intercept = FALSE, estimate = "CV")
} 

C <- as.numeric(vintages[vintages != "Barolo"])
X <- wines[vintages != "Barolo",]
wines2.df <- data.frame(vintages = C, wines = X)
wines2.lm0 <- lm(vintages ~ 1, data = wines2.df)
add1(wines2.lm0, scope = names(wines2.df)[-1])

wines2.lmfull <- lm(vintages ~ ., data = wines2.df)
drop1(wines2.lmfull)

step(wines2.lmfull)

if (leaps.present) {
  wines2.leaps <- regsubsets(vintages ~ ., data = wines2.df)
  wines2.leaps.sum <- summary(wines2.leaps)
  names(which(wines2.leaps.sum$which[8,]))
}

data(wines, package = "ChemometricsWithRData")
X <- wines[vintages != "Barolo", ]
vint <- factor(vintages[vintages != "Barolo"])
wines.counts <- table(vint)
wines.groups <- split(as.data.frame(X), vint)
wines.covmats <- lapply(wines.groups, cov)
wines.wcovmats <- lapply(1:length(wines.groups),
                         function(i, x, y) x[[i]]*y[i],
                         wines.covmats, wines.counts)
wines.pcov12 <- Reduce("+", wines.wcovmats) / (length(vint) - 2)
WSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x) {
                       crossprod(scale(x, scale = FALSE))}))
BSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x, y) {
                       nrow(x) * tcrossprod(colMeans(x) - y)},
                     colMeans(X)))
Tii <- solve(WSS + BSS)
MLLDA <-
  solve(wines.pcov12,
        apply(sapply(wines.groups, colMeans), 1, diff))
Ddist <- mahalanobis(colMeans(wines.groups[[1]]),
                     colMeans(wines.groups[[2]]),
                     wines.pcov12)
m <- sum(sapply(wines.groups, nrow)) - 2
p <- ncol(wines)
c <- prod(sapply(wines.groups, nrow)) / 
  sum(sapply(wines.groups, nrow))
Fcal <- (MLLDA^2 / diag(Tii)) * 
  (m - p + 1) * c^2 / (m * (m + c^2 * Ddist))
Fcal ## error: output missing in book

which(Fcal > qf(.95, 1, m-p+1))

## LASSO
if (lars.present & elasticnet.present) {
  odd <- seq(1, nrow(gasoline$NIR), by = 2)
  even <- seq(2, nrow(gasoline$NIR), by = 2)
  gas.lasso <- lars(gasoline$NIR[odd,], gasoline$octane[odd])
  plot(gas.lasso, xlim = c(0, .2))
  plot(gas.lasso, breaks = FALSE)
  plot(gas.lasso, breaks = FALSE, 
       xvar = "step", xlim = c(0, 20))
  
  set.seed(7)
  gas.lasso.cv <- cv.lars(gasoline$NIR[odd,], 
                          gasoline$octane[odd])
  best <- which.min(gas.lasso.cv$cv)
  abline(v = gas.lasso.cv$fraction[best], col = "red")
  
  gas.lasso.pred <- predict(gas.lasso, gasoline$NIR[even,], 
                            s = best)
  rms(gas.lasso.pred$fit, gasoline$octane[even])
  
  gas.lasso.coefs <- predict(gas.lasso, gasoline$NIR[even,],
                             s = best, type = "coef")
  gas.lasso.coefs$coefficients[gas.lasso.coefs$coeff != 0]
  
  ## elastic nets
  gas.enet <- enet(gasoline$NIR, gasoline$octane, lambda = .5)
  plot(gas.enet, "step")
  gas.enet.cv <- cv.enet(gasoline$NIR, gasoline$octane,
                         lambda = .5, s=1:50, mode="step")
  
  gas.enet.coef <- predict(gas.enet, gasoline$NIR[even,], s = best,
                           type = "coef")
  plot(wavelengths, gas.lasso.coefs$coefficients, type = "n",
       ylab = "Coefficients", xlab = "Wavelength (nm)",
       ylim = c(-60, 35), axes = FALSE)
  box()
  axis(1)
  segments(wavelengths, gas.enet.coef$coefficients - 40,
           wavelengths, -40, lwd = 2, col = "red")
  lines(wavelengths, gas.lasso.coefs$coefficients, lwd = 2,
        type = "h")
  legend("topleft", legend = c("Lasso", "Elastic net"),
         lty = 1, col = 1:2, bty = "n") 
}

## Simulated annealing
sbst <- SAstep(NULL, 13)
sbst

(sbst <- SAstep(sbst, 13))

(sbst <- SAstep(sbst, 13))

C <- factor(vintages[vintages != "Barolo"])
X <- wines[vintages != "Barolo",]

set.seed(7)
SAobj <- SAfun(X, C, lda.loofun, Tinit = 1)                 
SAobj

set.seed(7)
SAobj2 <- SAfun2(X, C, lda.loofun, Tinit = 1)
niter <- 100
cols <- rep(1, niter)
cols[!SAobj2$accepts] <- 2
plot(0:niter, SAobj2$qualities, col = "gray",
     type = "b", pch = cols,
     ylab = "Quality", xlab = "Iteration",
     main = "SA subset selection")
lines((0:niter)[SAobj2$accepts], SAobj2$qualities[SAobj2$accepts], type = "b")

SAobj <- SAfun(gasoline$NIR, gasoline$octane,
               eval.fun = pls.cvfun, Tinit = 3,
               fraction = .02, niter = 1000, ncomp = 2)
length(SAobj$best)

sqrt(-SAobj$best.q)

if (subselect.present) {
  winesHmat <- ldaHmat(X, C)
  wines.anneal <- 
    anneal(winesHmat$mat, kmin = 3, kmax = 3,
           H = winesHmat$H, criterion = "ccr12", r = 1)
  wines.anneal$bestsets
  
  wines.anneal$bestvalues
  
  ccr12.coef((nrow(X) - 1) * var(X), winesHmat$H, 
             r = 1, c(7, 10, 11))
}

lda.loofun(X, C, c(2, 7, 10))

## Genetic algorithms
set.seed(7)
pop1 <- GA.init.pop(pops = 5, nvar = 13, kmin = 2, kmax = 4)
pop1

pop1.q <- sapply(pop1, function(subset) lda.loofun(X, C, subset))
pop1.q

set.seed(7)
GA.select(pop1, 2, pop1.q, qlt.exp = 0)

GA.select(pop1, 2, pop1.q)

GA.XO(pop1[[1]], pop1[[2]])

GA.mut(pop1[[1]], 13, 1)

C <- factor(vintages[vintages != "Barolo"])
X <- wines[vintages != "Barolo",]
set.seed(117)
GAobj <- GAfun(X, C, lda.loofun, kmin = 3, kmax = 3)
GAobj

set.seed(7)
GAobj <- GAfun(gasoline$NIR, gasoline$octane, ncomp = 2,
               eval.fun = pls.cvfun, kmin = 3, kmax = 25)
GAobj$best
sqrt(-GAobj$best.q)

set.seed(7)
GAobj2 <- GAfun2(gasoline$NIR, gasoline$octane,
                 eval.fun = pls.cvfun,
                 kmin = 3, kmax = 25, ncomp = 2)
ngen <- nrow(GAobj2$qualities)
matplot(1:ngen - 1, GAobj2$qualities[1:ngen,], pch = 1:3,
        col = 1:3,
        ylab = "Quality", xlab = "Generation", type = "b",
        ylim = c(-0.25, 0),
        main = "GA subset selection")
grid()
legend("top", legend = c("best", "median", "worst"), bty = "n",
       lty = 3:1, pch = 3:1, col = 1:3)

wines.genetic <-
  genetic(winesHmat$mat, kmin = 3, kmax = 5, nger = 20,
          popsize = 50, maxclone = 0,
          H = winesHmat$H, criterion = "ccr12", r = 1)
wines.genetic$bestvalues

wines.genetic$bestsets
