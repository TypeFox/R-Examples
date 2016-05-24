
 #### copied from help page for function library

 pkg <- "aster"
 library(pkg, character.only = TRUE)

 set.seed(42)

 #### copied from help page for function aster in package aster

 data(echinacea)
 vars <- c("ld02", "ld03", "ld04", "fl02", "fl03", "fl04",
     "hdct02", "hdct03", "hdct04")
 redata <- reshape(echinacea, varying = list(vars), direction = "long",
     timevar = "varb", times = as.factor(vars), v.names = "resp")
 redata <- data.frame(redata, root = 1)
 pred <- c(0, 1, 2, 1, 2, 3, 4, 5, 6)
 fam <- c(1, 1, 1, 1, 1, 1, 3, 3, 3)
 hdct <- grep("hdct", as.character(redata$varb))
 hdct <- is.element(seq(along = redata$varb), hdct)
 redata <- data.frame(redata, hdct = as.integer(hdct))
 aout.foo <- aster(resp ~ varb + nsloc + ewloc + pop : hdct,
     pred, fam, varb, id, root, data = redata, type = "conditional")
 summary(aout.foo)

 beta <- aout.foo$coefficients
 dbeta <- rnorm(length(beta))

 pout <- predict(aout.foo, parm.type = "canonical",
     model.type = "conditional", se.fit = TRUE)
 theta <- pout$fit
 dtheta <- as.vector(pout$gradient %*% dbeta)

 pout <- predict(aout.foo, parm.type = "canonical",
     model.type = "unconditional", se.fit = TRUE)
 phi <- pout$fit
 dphi <- as.vector(pout$gradient %*% dbeta)

 pout <- predict(aout.foo, parm.type = "mean.value",
     model.type = "unconditional", se.fit = TRUE)
 mu <- pout$fit
 dmu <- as.vector(pout$gradient %*% dbeta)

 phony <- matrix(1, nrow = nrow(aout.foo$x), ncol = ncol(aout.foo$x))
 pout <- predict.aster(aout.foo, x = phony, root = aout.foo$root,
     modmat = aout.foo$modmat, parm.type = "mean.value",
     model.type = "conditional", se.fit = TRUE)
 xi <- pout$fit
 dxi <- as.vector(pout$gradient %*% dbeta)

 offset <- as.vector(aout.foo$origin)
 modmat <- matrix(aout.foo$modmat, ncol = length(beta))

 #### copied from help page for function library

 detach(pos = match(paste("package", pkg, sep=":"), search()))

 #### end of stuff from old aster package

 rm(list = setdiff(ls(), c("beta", "theta", "phi", "xi", "mu", "tau",
     "offset", "modmat", "dbeta", "dtheta", "dphi", "dxi", "dmu", "dtau")))

 library(aster2)

 data(echinacea)

 #### saturated

 myphi <- transformSaturated(theta, echinacea, from = "theta", to = "phi")
 all.equal(phi, myphi)

 mytheta <- transformSaturated(phi, echinacea, from = "phi", to = "theta")
 all.equal(theta, mytheta)

 myxi <- transformSaturated(theta, echinacea, from = "theta", to = "xi")
 all.equal(xi, myxi)

 mymu <- transformSaturated(xi, echinacea, from = "xi", to = "mu")
 all.equal(mu, mymu)

 #### unconditional from == "beta"

 phi.foo <- transformConditional(beta, modmat, echinacea,
     from = "beta", to = "phi", offset = offset)
 all.equal(phi, phi.foo)

 theta.foo <- transformConditional(beta, modmat, echinacea,
     from = "beta", to = "theta", offset = offset)
 all.equal(theta, theta.foo)

 xi.foo <- transformConditional(beta, modmat, echinacea,
     from = "beta", to = "xi", offset = offset)
 all.equal(xi, xi.foo)

 mu.foo <- transformConditional(beta, modmat, echinacea,
     from = "beta", to = "mu", offset = offset)
 all.equal(mu, mu.foo)

 #### unconditional from == "beta" (differential)

 dphi.foo <- transformConditional(beta, modmat, echinacea,
     from = "beta", to = "phi", offset = offset, differential = dbeta)
 all.equal(dphi, dphi.foo)

 dtheta.foo <- transformConditional(beta, modmat, echinacea,
     from = "beta", to = "theta", offset = offset, differential = dbeta)
 all.equal(dtheta, dtheta.foo)

 dxi.foo <- transformConditional(beta, modmat, echinacea,
     from = "beta", to = "xi", offset = offset, differential = dbeta)
 all.equal(dxi, dxi.foo)

 dmu.foo <- transformConditional(beta, modmat, echinacea,
     from = "beta", to = "mu", offset = offset, differential = dbeta)
 all.equal(dmu, dmu.foo)

