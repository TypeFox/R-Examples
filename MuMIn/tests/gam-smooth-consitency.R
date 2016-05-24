
library(MuMIn)
suppressPackageStartupMessages(library(mgcv))
RNGkind("Mersenne")


set.seed(0) ## simulate some data...
dat <- gamSim(1, n = 200, dist = "binary", scale = 2)


gamsmoothwrap <-
function(formula, k1 = NA, ...) { 
    cl <- origCall <- match.call()
    cl[[1]] <- as.name("gam")
    cl$formula <- exprApply(formula, "s", function(e, k, x) {
        i <- which(e[[2]] == x)[1]
        if(!is.na(i) && !is.na(k[i])) e[["k"]] <- k[i]
        e
    }, k = c(k1), x = c("x0"))
    cl$k1 <- NULL
	fit <- eval(cl, parent.frame())
    fit$call <- origCall # replace the stored call
    fit
}

glo <- gamsmoothwrap(y ~ s(x0, fx = TRUE), k1 = 3, data = dat, family = binomial)


print(ms <- model.sel(models <- lapply(dredge(glo, varying = list(k1 = 3:16), fixed = TRUE,
	evaluate = FALSE), eval)))


# should throw a warning
stopifnot(isTRUE(tryCatch(model.avg(models), warning = function(e) TRUE)))

#_______________________________________________________________________________

set.seed(0) ## simulate some data...
dat <- gamSim(1, n = 400, dist = "binary", scale = 2)

testSmoothKConsistency <- MuMIn:::testSmoothKConsistency


glm1 <- glm(y ~ x1 * x2, data = dat, family = binomial)
gam53 <- gam(y ~ s(x0, k = 5) + s(x1, k = 3), data = dat, family = binomial)
gam37 <- update(gam53, . ~ s(x0, k = 3) + s(x1, k = 7))
gam3e44 <- update(gam53, . ~ s(x0, k = 3) + te(x2, x1, k = 4))
gam3e54 <- update(gam53, . ~ s(x0, k = 3) + te(x2, x1, k = c(5, 4)))
gam3i54 <- update(gam53, . ~ s(x0, k = 3) + ti(x2, x1, k = c(5, 4)))

gam3e54 <- update(gam53, . ~ s(x0, k = 3) + te(x2, x1, k = c(5, 4)))
gam3i34 <- update(gam53, . ~ s(x0, k = 3) + ti(x2, x1, k = c(3, 4)))


ms <- model.sel(gam53, gam37, gam3e44, gam3e54, gam3i54, glm1)
testSmoothKConsistency(ms)
testSmoothKConsistency(list(gam3e54, gam3i34))



#_______________________________________________________________________________


