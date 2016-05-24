## 
## 9.2.1 A Single-Factor Model
##

# Note: S-Plus code in book works on R except that
# input is from .rda files in package FinTS and these 
# changes:
#
# S-Plus: solve(x, y) # non-square x
# R: solve.qr(x, y)
#
# S-Plus: var(x, SumSquares=T)
# R:  (nrow(x) - 1) * var(x)
library(FinTS)

# input returns file
data(m.fac9003)
x <- m.fac9003

# Table 9.1 showing means and standard deviations
# of stocks
cbind(mean = colMeans(x[,-14]), sd = sd(x[,-14]))

# Code at bottom of page 408 producing
# table middle of p 409 showing beta.hat, sigma
# and R squared.
xmtx <- cbind(1, as.numeric(x[, 14]))
rtn <- x[, -14]
xit.hat <- qr.solve(xmtx, rtn)
beta.hat <- t(xit.hat[2,])
E.hat <- rtn - xmtx %*% xit.hat
D.hat <- diag(crossprod(E.hat)/(168-2))
r.square <- 1 - (168 - 2) * D.hat / diag((168 -1 ) * var(rtn))
tab <- cbind(beta.hat = c(beta.hat), `sigma(i)` = sqrt(D.hat), r.square)
print(tab, digits = 3)

# Same as last snippet but making use of lm() and summary()
f <- function(y) {
    L <- lm(y ~ x[,14])
    S <- summary(L)
    c(beta.hat = coef(L)[[2]], sigma = S$sigma, r.squared = S$r.squared)
}
tab <- t(apply(x[,-14], 2, f))
print(tab, digits = 3)

# Figure 9.1 Bar chart of columns 1 and 3 of tab
opar <- par(mfrow = c(2, 1))
barplot(tab[, 1], main = quote(hat(beta)))
barplot(tab[, 3], main = quote(R^2))
par(opar)

# Code on bottom of page 409 and on page 410
cov.r <- var(x[, 14]) * (t(beta.hat) %*% beta.hat) + diag(D.hat)
# sd.r <- sqrt(diag(cov.r))
# corr.r <- cov.r / outer(sd.r, sd.r)
corr.r <- cov2cor(cov.r)
print(corr.r, digits = 1, width = 2)
print(cor(rtn), digits = 1, width = 2)

# Code on middle of page 411
# w.gmin.model <- solve(cov.r) %*% rep(1, nrow(cov.r))
# w.gmin.model <- w.gmin.model / sum(w.gmin.model)
w.gmin.model <- prop.table(solve(cov.r) %*% rep(1, nrow(cov.r)))
print(t(w.gmin.model), 4)

# w.gmin.data <- solve(var(rtn)) %*% rep(1, nrow(cov.r))
# w.gmin.data <- w.gmin.data / sum(w.gmin.data)
w.gmin.data <- prop.table(solve(var(rtn)) %*% rep(1, nrow(cov.r)))
print(t(w.gmin.data), width = 4)

# Code at bottom of page 411 and top of page 412
resi.cov <- t(E.hat) %*% E.hat/(168 - 2)
# resi.sd <- sqrt(diag(resi.cov))
# resi.cor <- resi.cov / outer(resi.sd, resi.sd)
resi.cor <- cov2cor(resi.cov)
print(resi.cor, digits = 1, width = 2)

## 
## 9.2.2 Multifactor Models
##

### ??
### pages 412-413
### note that the beta.hat and r.square values
###  do not match Figure 9.2
### According to comment near top of:
###   http://www.estima.com/textbooks/tsayp412.prg
### beta.hat also does not match when using RATS software.

data(m.cpice16.dp7503)
y1 <- m.cpice16.dp7503
library(vars)
# the VARselect() value of SC and 
# S-Plus VAR() value of BIC are 
# translated and scaled versions of each other
vs <- VARselect(y1, lag.max = 13)
print(vs)
var3.fit <- VAR(y1, 3)
res <- resid(var3.fit)[166:333,]
data(m.fac9003)
da <- m.fac9003
# xmtx <- cbind(1, res)
rtn <- da[, -14]
# xit.hat <- qr.solve(xmtx, rtn)
reg <- lm(rtn ~ res)
xit.hat <- coef(reg)
beta.hat <- t(xit.hat[2:3,])
# E.hat <- rtn - xmtx %*% xit.hat
E.hat <- resid(reg)
D.hat <- diag(crossprod(E.hat)/(168-3))
# r.square <- 1 - (168-3) * D.hat / diag((nrow(rtn) - 1) * var(rtn))
r.square <- sapply(summary(reg), "[[", "r.squared")

cov.rtn <- beta.hat %*% var(res) %*% t(beta.hat) + diag(D.hat)
# sd.rtn <- sqrt(diag(cov.rtn))
# cor.rtn <- cov.rtn / outer(sd.rtn, sd.rtn)
cor.rtn <- cov2cor(cov.rtn)
print(cor.rtn, digits = 1, width = 2)

# cov.resi <- t(E.hat) %*% E.hat / (168-3)
cov.resi <- crossprod(E.hat) / (168-3)
# sd.resi <- sqrt(diag(cov.resi))
# cor.resi <- cov.resi / outer(sd.resi, sd.resi)
cor.resi <- cov2cor(cov.resi)
print(cor.resi, digits = 1, width = 2)

### ??
### there is a comment in the book that cor.rtn and cor.resi are close
### but they don't seem to be here

##
# 9.3 Fundamental Factor Models
# code on page 417
##
data(m.barra.9003)
da <- m.barra.9003

# demean
rtn.rm <- scale(da, scale = FALSE) 
# fin <- c(rep(1, 4), rep(0, 6))
# tech <- c(rep(0, 4), rep(1, 3), rep(0, 3))
# oth <- c(rep(0, 7), rep(1, 3))
# ind.dum <- cbind(fin, tech, oth)
inds <- c("ind", "tech", "oth")
ind.fac <- factor(rep(inds, c(4, 3, 3)), levels = inds)
ind.dum <- model.matrix(~ ind.fac - 1)
# rtn <- t(rtn.rm)
# cov.rtn <- var(rtn)
# sd.rtn <- sqrt(diag(cov.rtn))
# corr.rtn <- cov.rtn / outer(sd.rtn, sd.rtn)
corr.rtn <- cov2cor(var(rtn.rm))
print(corr.rtn, digit = 1, width = 2)

# code on page 418
# F.hat.o <- solve(crossprod(ind.dum)) %*% t(ind.dum) %*% t(rtn.rm)
# E.hat.o <- t(rtn.rm) - ind.dum %*% F.hat.o
F.hat.o <- lm(t(rtn.rm) ~ ind.dum - 1)
E.hat.o <- resid(F.hat.o)

diagD.hat.o <- apply(E.hat.o, 1, var)

Dinv.hat <- diag(1/diagD.hat.o)
H1 <- t(ind.dum) %*% Dinv.hat %*% ind.dum
Hmtx <- solve(H1) %*% t(ind.dum) %*% Dinv.hat
F.hat.g <- Hmtx %*% t(rtn.rm)
F.hat.gt <- t(F.hat.g)
E.hat.g <- t(rtn.rm) - ind.dum %*% F.hat.g
diagD.hat.g <- apply(E.hat.g, 1, var)
print(t(Hmtx), digit = 4)

cov.ind <- ind.dum %*% var(F.hat.gt) %*% t(ind.dum) + diag(diagD.hat.g)
# sd.ind <- sqrt(diag(cov.ind))
# corr.ind <- cov.ind(outer(sd.ind, sd.ind))
corr.ind <- cov2cor(cov.ind)
print(corr.ind, digits = 1, width = 2)


##
# 9.4 Principal Component Analysis
##

# page 423

data(m.5cln)
rtn <- m.5cln
round(rtn[1,], 2)
round(colMeans(rtn), 2)
round(cov(rtn), 2)
round(cor(rtn), 2)
eigen(cor(rtn))

# Figure 9.4 on page 423

plot(rtn, nc = 1)

# table 9.3 on page 425

# eigenvalues and eigenvectors
# Eigenvectors are only determined up to scale so the ratio of any
# eigenvector to the one shown in the book should be constant.

e.cov <- eigen(cov(rtn))
prop.cov <- prop.table(e.cov$values)
t(cbind(Eigenvalue = e.cov$values, Proportion = prop.cov, 
    Cumulative = cumsum(prop.cov), e.cov$vectors))

e.cor <- eigen(cor(rtn))
prop.cor <- prop.table(e.cor$values)
t(cbind(Eigenvalue = e.cor$values, Proportion = prop.cor, 
    Cumulative = cumsum(prop.cor), e.cov$vectors))

pca.cor <- princomp(rtn, cor = TRUE)
pca.cor$sdev^2 # eigenvalues
pca.cor$loadings # columns are eigenvectors

# page 425
pca.cov <- princomp(rtn)
names(pca.cov)
summary(pca.cov)
pca.cov$loadings
screeplot(pca.cov, type = "l")

pca.cor <- princomp(rtn, cor = TRUE)
names(pca.cor)
summary(pca.cor)
pca.cor$loadings
screeplot(pca.cor, type = "l")

##
# 9.5 - Statistical Factor Analysis
##

# example 9.3 on pages 431-433 
# Note reference in text to example 9.2 should read 8.2

data(m.bnd)
cor(m.bnd)

# the next 4 give PC unrotated factor loadings as in table 9.5

# e$vec is normalized so that its length is 1.  Multiply that
# by sqrt of eigenvalue and multiply each column by -1 if needed
# to keep diagonal elements +ve.
e <- eigen(cor(m.bnd))
scale(e$vec, center = FALSE, scale = sign(diag(e$vec))/sqrt(e$val))[, 1:2]

# same
e <- eigen(cor(m.bnd))
e$vec %*% diag(sign(diag(e$vec)) * sqrt(e$val))

# same
p <- princomp(m.bnd, cor = TRUE)
ld <- loadings(p)  # these are eigenvectors
scale(ld, center = FALSE, scale = sign(diag(ld)) / p$sdev) 

# same
p <- princomp(m.bnd, cor = TRUE, center = FALSE) 
ld <- loadings(p)  # these are eigenvectors
ld %*% diag(sign(diag(ld)) * p$sdev)

# Variances
colSums(ld^2)

# Proportion Variances
colMeans(ld^2)

# cumumalative proportion variances
cumsum(colMeans(ld^2))

# this gives the same results as upper portion of Table 9.5
# First factor.pa gives unrotated and second factor.pa gives rotated
library(psych)
bnd.fac <- factor.pa(m.bnd, nfactors = 2, rot = "none")
bnd.fac
bnd.fac <- factor.pa(m.bnd, nfactors = 2)
bnd.fac

# commonalities -- result does not depend on rotation
rowSums(loadings(bnd.fac)[, 1:2]^2)

# NOTE !!!
# this does not give same result as lower left portion of Table 9.5
factanal(m.bnd, factors = 2, rotation = "none")
# but this does give same results as lower right portion of Table 9.5
factanal(m.bnd, factors = 2)

# commonalities -- result does not depend on rotation
rowSums(loadings(bnd.fac)[, 1:2]^2)

# example 9.4 page 433
data(m.barra.9003)
rtn <- m.barra.9003
stat.fac <- factanal(rtn, factors = 2)
stat.fac
names(stat.fac)

stat.fac <- factanal(rtn, factors = 3)
stat.fac

# Figure 9.6 on page 435
# plot top n loadings of the first k factors
plot(loadings(stat.fac))

# bottom of page 434
library(GPArotation)
stat.fac2 <- quartimax(loadings(stat.fac), normalize = TRUE)
round(loadings(stat.fac2), 3)

# page 435
factor.real <- factanal(rtn, scores = "Bartlett", factors = 3)$scores
print(stat.fac$correlation, digits = 1, width = 2)

##
# 9.6 - Asymptotic Principal Component Analysis
##

# page 438
data(m.apca0103)
M.apca0103 <- xtabs(return ~ date + CompanyID, m.apca0103)
nf1 <- apca(M.apca0103, 1)
nf6 <- apca(M.apca0103, 6)

# screeplot - Figure 9.7 on page 439
prop <- prop.table(nf6$eig)[1:7]
bp <- barplot(prop, ylim = c(0, .5))
text(bp, prop, round(prop, 3), pos = 3)


