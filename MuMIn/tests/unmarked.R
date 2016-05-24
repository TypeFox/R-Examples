if(length(find.package("unmarked", quiet = TRUE)) == 1L) {

library(MuMIn)
library(stats4)
library(unmarked)
options(na.action = "na.fail")

# Simulate occupancy data

set.seed(1)
umfOccu <- local({
	nSites <- 100
	nReps <- 5
	covariates <- data.frame(veght=rnorm(nSites),
		habitat = factor(c(rep('A', 50), rep('B', 50))))
	psipars <- c(-1, 1, -1)
	ppars <- c(1, -1, 0)
	X <- model.matrix(~veght+habitat, covariates) # design matrix
	psi <- plogis(X %*% psipars)
	p <- plogis(X %*% ppars)
	y <- matrix(NA, nSites, nReps)
	Z <- rbinom(nSites, 1, psi)       # true occupancy state
	for(i in 1:nSites) y[i,] <- rbinom(nReps, 1, Z[i]*p[i])
	unmarkedFrameOccu(y = y, siteCovs = covariates)
})

# Fit some models
fm1oc <- occu(~1 ~1, umfOccu)
fm2oc <- occu(~veght+habitat ~veght*habitat, umfOccu)
fm3oc <- occu(~habitat ~veght+habitat, umfOccu)
fm4oc <- occu(~veght ~veght+habitat, umfOccu)

#ms <- model.sel(fm1oc, fm2oc, fm3oc, fm4oc)
#summary(print(model.avg(ms)))

data(linetran)
ltUMF <- with(linetran, {
   unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
   siteCovs = data.frame(Length, area, habitat),
   dist.breaks = c(0, 5, 10, 15, 20),
   tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
   })

fm2 <- distsamp(~area * habitat ~ habitat, ltUMF)

dd <- dredge(fm2, fixed = ~p(sigmaarea))
mo <- get.models(dd, T)
avg <- model.avg(mo)

#formula(fm2)
#names(fm2)
#getAllTerms(fm2)
#coef(fm2)
#linearComb(fm2, matrix(c(1,1,1, 1, 2,2,2,2), ncol = 4, nrow = 2, byrow = T), "det")


#avgpred(avg, use.lincomb = TRUE, type = "link")

#MuMIn:::avgpred(avg)

#getMethod("predict", "unmarkedFit")

#yall <- lapply(mo, std_predict)
#z <- yall[[1]]
#yy1 <- avgterms(yall, Weights(dd), revised.var = T, full = F, dfs = NULL)
#yy2 <- avgsefit(lapply(yall, lapply,  "[", T , 1), Weights(dd), revised.var = T, full = F, dfs = NULL)

# Bring in Data
library(unmarked)
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)),sitevar2 = rnorm(numSites(pferUMF)))

global <- occu(~ sitevar1 + I(sitevar1^2) + sitevar2 ~ 1, pferUMF)
dd1 <- dredge(global) # <- has models containing quadratic 'I(sitecar1^2)' without corresponding 'sitevar1' (not what I want)

getAllTerms(global)

msubset <- expression(`p(sitevar1)` || `p(I(sitevar1^2))`)
dd2 <- dredge(global, subset = msubset)

# dredge(fm2oc, eval=F, fixed=~psi(habitat))

#(dd <- dredge(fm2oc, fixed = ~psi(habitat), trace = T))
(dd <- dredge(fm2oc, fixed = ~psi(habitat)))

#MuMIn:::avgpred(model.avg(dd, fit = T))

# mod <- get.models(dd, 1:5)

# ma <- model.avg(mod)
## predict(ma, umfOccu, type = "state")

# Checking if 'dc' works properly in both subset.model.selection and dredge.
dd2 <- dredge(fm2oc, subset = `psi(habitat)` & dc(`p(habitat)`, `p(veght)`))
dd1a <- subset(dd, has(psi(habitat)) & dc(p(habitat), p(veght)), recalc.delta = TRUE)
stopifnot(all.equal(dd2, dd1a, check.attributes = FALSE))
stopifnot(!any(is.na(dd2[, "p(habitat)"]) & !is.na(dd2[, "p(veght)"])))



model.sel(dd, rank = "AIC")
models <- get.models(dd, subset = 1:3)

summary(ma1 <- model.avg(models))
summary(ma2 <- model.avg(dd[1:3]))
summary(ma3 <- model.avg(model.sel(model.sel(dd, rank = "AIC"), rank = "AICc")[1:3]))

predict(ma1, type = "det")
# MuMIn:::avgpred(ma1)$fit[, "det"]

stopifnot(!any(is.na(coefTable(ma1)[, 1:2])))
stopifnot(!any(is.na(coefTable(ma2)[, 1:2])))
stopifnot(!any(is.na(coefTable(ma3)[, 1:2])))
stopifnot(isTRUE(all.equal(coefTable(ma1), coefTable(ma2))))
stopifnot(isTRUE(all.equal(coefTable(ma2), coefTable(ma3))))
stopifnot(!any(is.na(dd[, "psi(habitat)"])))

summary(model.avg(dd, delta <= 4))

#family.unmarkedFit <- function (object, ...) NA

# Model selection
print(model.sel(fm1oc, fm2oc, fm3oc))


#models <- list(fm1oc, fm2oc, fm3oc)
#traceback()
#model.avg(fm1oc, fm2oc, fm3oc)
#model.avg(fm1oc, fm2oc)
#getAllTerms(fm1oc)

print(summary(model.avg(fm1oc, fm2oc, fm3oc)))

}
