
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
 aout4 <- aster(resp ~ varb + nsloc + ewloc + pop * hdct - pop,
     pred, fam, varb, id, root, data = redata)
 summary(aout4)

 beta <- aout4$coefficients
 dbeta <- rnorm(length(beta))

 pout <- predict(aout4, parm.type = "canonical", model.type = "conditional",
     se.fit = TRUE)
 theta <- pout$fit
 dtheta <- as.vector(pout$gradient %*% dbeta)
 jack.theta <- pout$gradient

 pout <- predict(aout4, parm.type = "canonical", model.type = "unconditional",
     se.fit = TRUE)
 phi <- pout$fit
 dphi <- as.vector(pout$gradient %*% dbeta)
 jack.phi <- pout$gradient

 pout <- predict(aout4, parm.type = "mean.value", model.type = "unconditional",
     se.fit = TRUE)
 mu <- pout$fit
 dmu <- as.vector(pout$gradient %*% dbeta)
 jack.mu <- pout$gradient

 phony <- matrix(1, nrow = nrow(aout4$x), ncol = ncol(aout4$x))
 pout <- predict.aster(aout4, x = phony, root = aout4$root,
     modmat = aout4$modmat, parm.type = "mean.value",
     model.type = "conditional", se.fit = TRUE)
 xi <- pout$fit
 dxi <- as.vector(pout$gradient %*% dbeta)
 jack.xi <- pout$gradient

 offset <- as.vector(aout4$origin)
 modmat <- matrix(aout4$modmat, ncol = length(beta))
 tau <- as.numeric(t(modmat) %*% mu)
 dtau <- as.numeric(t(modmat) %*% dmu)

 modmat.orig <- model.matrix(resp ~ varb + nsloc + ewloc + pop * hdct - pop,
     redata)
 beta.orig <- rep(0, ncol(modmat.orig))
 beta.orig[match(names(beta), colnames(modmat.orig))] <- beta
 names(beta.orig) <- colnames(modmat.orig)
 tau.orig <- as.numeric(t(modmat.orig) %*% mu)

 #### copied from help page for function library

 detach(pos = match(paste("package", pkg, sep=":"), search()))

 #### end of stuff from old aster package

 rm(list = setdiff(ls(), c("beta", "theta", "phi", "xi", "mu", "tau",
     "offset", "modmat", "dbeta", "dtheta", "dphi", "dxi", "dmu", "dtau",
     "jack.theta", "jack.phi", "jack.xi", "jack.mu",
     "modmat.orig", "beta.orig", "tau.orig")))

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

 phi.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "phi", offset = offset)
 all.equal(phi, phi.foo)

 theta.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "theta", offset = offset)
 all.equal(theta, theta.foo)

 xi.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "xi", offset = offset)
 all.equal(xi, xi.foo)

 mu.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "mu", offset = offset)
 all.equal(mu, mu.foo)

 tau.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "tau", offset = offset)
 all.equal(tau, tau.foo)

 #### unconditional from == "beta" (differential)

 dphi.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "phi", offset = offset, differential = dbeta)
 all.equal(dphi, dphi.foo)

 dtheta.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "theta", offset = offset, differential = dbeta)
 all.equal(dtheta, dtheta.foo)

 dxi.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "xi", offset = offset, differential = dbeta)
 all.equal(dxi, dxi.foo)

 dmu.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "mu", offset = offset, differential = dbeta)
 all.equal(dmu, dmu.foo)

 dtau.foo <- transformUnconditional(beta, modmat, echinacea,
     from = "beta", to = "tau", offset = offset, differential = dbeta)
 all.equal(dtau, dtau.foo)

 #### Jacobian matrices

 my.jack.theta <- jacobian(beta, echinacea, from = "beta", to = "theta",
     offset = offset, modmat = modmat, transform = "unconditional")
 all.equal(jack.theta, my.jack.theta)

 my.jack.phi <- jacobian(beta, echinacea, from = "beta", to = "phi",
     offset = offset, modmat = modmat, transform = "unconditional")
 all.equal(jack.phi, my.jack.phi)

 my.jack.xi <- jacobian(beta, echinacea, from = "beta", to = "xi",
     offset = offset, modmat = modmat, transform = "unconditional")
 all.equal(jack.xi, my.jack.xi)

 my.jack.mu <- jacobian(beta, echinacea, from = "beta", to = "mu",
     offset = offset, modmat = modmat, transform = "unconditional")
 all.equal(jack.mu, my.jack.mu)

 my.jack.tau <- jacobian(beta, echinacea, from = "beta", to = "tau",
     offset = offset, modmat = modmat, transform = "unconditional")
 all.equal(t(modmat) %*% jack.mu, my.jack.tau)

 #### unconditional from == "tau"

 beta.foo <- transformUnconditional(tau, modmat, echinacea,
     from = "tau", to = "beta", offset = offset)
 all.equal(as.vector(beta), beta.foo)

 xi.foo <- transformUnconditional(tau, modmat, echinacea,
     from = "tau", to = "xi", offset = offset)
 all.equal(as.vector(xi), xi.foo)

 beta.orig.foo <- transformUnconditional(tau.orig, modmat.orig, echinacea,
     from = "tau", to = "beta", offset = offset)
 all.equal(as.vector(beta.orig), beta.orig.foo)

 #### unconditional from == "tau" (differential)

 dbeta.foo <- transformUnconditional(tau, modmat, echinacea,
     from = "tau", to = "beta", offset = offset, differential = dtau)
 all.equal(dbeta, dbeta.foo)

 dxi.foo <- transformUnconditional(tau, modmat, echinacea,
     from = "tau", to = "xi", offset = offset, differential = dtau)
 all.equal(as.vector(dxi), dxi.foo)

 dbeta.orig <- rnorm(length(tau.orig))
 dtau.orig <- transformUnconditional(beta.orig, modmat.orig, echinacea,
     from = "beta", to = "tau", offset = offset, differential = dbeta.orig)
 dbeta.orig.foo <- transformUnconditional(tau.orig, modmat.orig, echinacea,
     from = "tau", to = "beta", offset = offset, differential = dtau.orig)
 # since betas are meaningless, cannot compare them
 dphi.orig <- as.numeric(modmat.orig %*% dbeta.orig)
 dphi.orig.foo <- as.numeric(modmat.orig %*% dbeta.orig.foo)
 is.same(dphi.orig, dphi.orig.foo, echinacea, parm.type = "phi")

