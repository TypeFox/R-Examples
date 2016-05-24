### R code from vignette source 'wjd.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: wjd.Rnw:23-24
###################################################
  if(exists(".orig.enc")) options(encoding = .orig.enc)


###################################################
### code chunk number 2: libraryCHNOSZ
###################################################
library(CHNOSZ)
data(thermo)


###################################################
### code chunk number 3: wjd
###################################################
w <- wjd()


###################################################
### code chunk number 4: abundant
###################################################
# the order of species abundance
oX <- order(w$X, decreasing=TRUE)
# the stoichiometries of the two most abundant species
w$A[head(oX,2),]
# the number of iterations
niter <- length(w$lambda)
niter


###################################################
### code chunk number 5: abundantnames
###################################################
# the formulas of the two most abundant species
# use CHNOSZ's as.chemical.formula function, after 
# dropping elements with zero stoichiometry
f1 <- as.chemical.formula(w$A[oX[1],w$A[oX[1],]!=0])
f2 <- as.chemical.formula(w$A[oX[2],w$A[oX[2],]!=0])
f4 <- as.chemical.formula(w$A[oX[4],w$A[oX[4],]!=0])


###################################################
### code chunk number 6: iterfun
###################################################
iterfun <- function(imax, i) {
  w <- wjd(imax=imax)
  x <- w$X[i]/sum(w$X)
  return(x)
}


###################################################
### code chunk number 7: wjd.Rnw:132-132
###################################################



###################################################
### code chunk number 8: iterplot
###################################################
par(mfrow=c(1, 2))
sa <- sapply(0:niter, iterfun, i=oX[1])
plot(0:niter, sa, xlab="iteration", ylab=paste("x", f1))
sa <- sapply(0:niter, iterfun, i=oX[4])
plot(0:niter, sa, xlab="iteration", ylab=paste("x", f4))


###################################################
### code chunk number 9: Gdown
###################################################
all(diff(w$F.Y) < 0)


###################################################
### code chunk number 10: Gdownslower
###################################################
identical(diff(w$F.Y), sort(diff(w$F.Y)))


###################################################
### code chunk number 11: Gfrac
###################################################
diff(w$F.Y)/w$F.Y[1:7]


###################################################
### code chunk number 12: wjd.Rnw:207-207
###################################################



###################################################
### code chunk number 13: w3_ep_plot
###################################################
w3 <- wjd(imax=3)
ep3 <- element.potentials(w3, plot.it=TRUE)


###################################################
### code chunk number 14: wjd.Rnw:220-220
###################################################



###################################################
### code chunk number 15: w_ep_plot
###################################################
ep <- element.potentials(w, plot.it=TRUE)


###################################################
### code chunk number 16: w_ep_plot
###################################################
is.near.equil(w3) # 3 iterations
is.near.equil(w)  # 7 iterations
is.near.equil(w, tol=0.00001)


###################################################
### code chunk number 17: guessw12
###################################################
as.list(args(wjd))$Y
Y1 <- guess(method="stoich")
Y1
#Y2 <- guess(method="central")
#Y2


###################################################
### code chunk number 18: wjd.Rnw:285-285
###################################################



###################################################
### code chunk number 19: stoich_guess
###################################################
wY1 <- wjd(Y=Y1)
niterY1 <- length(wY1$lambda)
niterY1
is.near.equil(wY1, tol=0.0001)


###################################################
### code chunk number 20: wjd.Rnw:305-305
###################################################



###################################################
### code chunk number 21: compare_guesses
###################################################
plot(1:10, w$X, xlab="species", ylab="mole fraction")
points(1:10, wY1$X, pch=0)
#points(1:10, wY2$X, pch=2)


###################################################
### code chunk number 22: allguesses
###################################################
Ys <- guess(iguess=NULL, method="stoich")
# total number of species combinations
length(Ys)
# species combinations that have all-positive mole numbers
iYs <- which(!is.na(Ys))
nguess <- length(iYs)
nguess


###################################################
### code chunk number 23: allguess.equil
###################################################
sapply(iYs,function(i) is.near.equil(wjd(Y=Ys[[i]])))


###################################################
### code chunk number 24: allguess.equil.2
###################################################
sapply(iYs,function(i) is.near.equil(wjd(Y=Ys[[i]]),tol=0.0001))


###################################################
### code chunk number 25: allguess.equil.3
###################################################
sapply(iYs, function(i) {
  is.near.equil(wjd(Y=Ys[[i]], Gfrac=1e-9), tol=0.0001)
})


###################################################
### code chunk number 26: sd.species
###################################################
Xs <- sapply(iYs, function(i) wjd(Y=Ys[[i]], Gfrac=1e-9)$X)
apply(Xs, 1, sd)


###################################################
### code chunk number 27: read_dlen
###################################################
# read formulas and Gibbs energies
file <- system.file("extdata/thermo/DLEN67.csv", package="CHNOSZ")
dlen <- read.csv(file, stringsAsFactors=FALSE, row.names=1)
t(dlen[, 1, drop=FALSE])


###################################################
### code chunk number 28: setup_atmos
###################################################
# turn formulas into a stoichiometric matrix
A <- i2A(dlen$formula)
# assemble Gibbs energies/RT at 500 K
G0.RT <- 1000 * dlen$G500 / thermo$opt$R / thermo$opt$Tr
# a function to minimize Gibbs energy for system with 
# given mole fraction of carbon (xC)
min.atmos <- function(xC) {
  # the bulk composition C:H:N:O
  B <- c(xC, 100-40-xC, xC, 40)
  # guess the initial composition
  Y <- guess(A, B)
  w <- wjd(A=A, G0.RT=G0.RT, Y=Y, P=1, imax=90, Gfrac=1e-14)
  if(!is.near.equil(w)) cat(paste("not near equilibrium for xC=", xC, "\n"))
  return(w)
}


###################################################
### code chunk number 29: atmos_C15
###################################################
sort(min.atmos(15)$X, decreasing=TRUE)


###################################################
### code chunk number 30: atmos_figure
###################################################
xCs <- seq(8, 47, 3)
Xs <- sapply(xCs, function(xC) min.atmos(xC)$X)
# normalize the mole numbers to mole fractions
Xs <- t(t(Xs)/colSums(Xs))
plot(-10, 0, xlim=c(0, 55), ylim=c(-25, 1), xlab="mole percent C", ylab="log10 mole fraction")
for(i in 1:nrow(Xs)) lines(xCs, log10(Xs[i, ]))
text(48, log10(Xs[, length(xCs)]), dlen$formula, adj=0)
text(7, log10(Xs[, 1]), dlen$formula, adj=1)


###################################################
### code chunk number 31: alkanes_aromatics
###################################################
alkanes <- c("n-hexane", "n-heptane", "n-octane", "n-nonane")
ialk <- info(alkanes, "liq")
aromatics <- c("naphthalene", "anthracene", "phenanthrene", "pyrene")
iaro <- info(aromatics, "liq")


###################################################
### code chunk number 32: run_aa
###################################################
waa <- run.wjd(c(ialk, iaro), Y=rep(1, 8))
waa$elements
is.near.equil(waa)


###################################################
### code chunk number 33: run_aa_equil
###################################################
waa <- run.wjd(c(ialk, iaro), Y=rep(1, 8), imax=20, Gfrac=1e-14, nlambda=501)
is.near.equil(waa)


###################################################
### code chunk number 34: barplot_aa
###################################################
bp <- barplot(waa$X, ylab="moles")
text(bp, rep(0.2,8), c(alkanes, aromatics), srt=90, adj=0)


###################################################
### code chunk number 35: element_potentials_C84H106
###################################################
ep <- equil.potentials(waa)
print(ep)


###################################################
### code chunk number 36: basis_potentials_C84H106
###################################################
basis(c("graphite", "H2"), c("cr", "gas"))
basis.logact(ep)


###################################################
### code chunk number 37: run_aa_C100H70
###################################################
waa <- run.wjd(c(ialk, iaro), B="C100H70")
bp <- barplot(waa$X, ylab="moles")
text(bp, rep(0.2,8), c(alkanes, aromatics), srt=90, adj=0)


###################################################
### code chunk number 38: basis_potentials_C84H106
###################################################
print(ep <- equil.potentials(waa))
basis.logact(ep)


###################################################
### code chunk number 39: sessionInfo
###################################################
sessionInfo()


