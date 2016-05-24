
 library(aster2)

 set.seed(42)

 data(test1)

 fred <- asterdata(test1,
     vars = c("m1", "m2", "m3", "n1", "n2", "b1", "p1", "z1"),
     pred = c(0, 0, 0, 1, 1, 2, 3, 6), group = c(0, 1, 2, 0, 4, 0, 0, 0),
     code = c(1, 1, 1, 2, 2, 3, 4, 5),
     families = list(fam.multinomial(3), "normal.location.scale",
     "bernoulli", "poisson", "zero.truncated.poisson"))

 ########## valid theta ##########

 theta <- rnorm(nrow(fred$redata)) * 0.1
 try(validtheta(fred, theta))

 must.be.negative <- as.character(fred$redata$varb) == "n2"
 theta[must.be.negative] <- (- abs(theta[must.be.negative]))
 is.validtheta(fred, theta)

 ########## theta to phi ##########

 phi <- transformSaturated(theta, fred, from = "theta", to = "phi")

 group <- fred$regroup
 group.idx <- seq(along = group)
 group.grp <- cbind(group.idx)
 repeat {
    group.idx <- c(0, group)[group.idx + 1]
    if (all(group.idx == 0)) break
    group.grp <- cbind(group.grp, group.idx)
 }
 outies <- sort(unique(as.vector(group.grp[ , -1])))
 outies <- outies[outies != 0]
 if (length(outies) > 0)
     group.grp <- group.grp[- outies, ]
 group.grp <- as.data.frame(t(group.grp))
 group.grp <- as.list(group.grp)
 names(group.grp) <- NULL
 group.grp <- lapply(group.grp, function(x) x[x != 0])
 group.grp <- lapply(group.grp, rev)
 identical(as.numeric(seq(along = group)), as.numeric(sort(unlist(group.grp))))

 delta <- fred$redelta
 families <- fred$families[fred$recode]
 z <- sapply(group.grp, function(i) cumulant(theta[i],
     families[[i[1]]], delta = delta[i])$zeroth)

 group.pred <- fred$repred[sapply(group.grp, function(x) x[1])]

 myphi <- theta
 for (i in seq(along = group.grp)) {
      j <- group.pred[i]
      if (j != 0) {
          myphi[j] <- myphi[j] - z[i]
      }
 }

 all.equal(phi, myphi)

 ########## theta to phi derivative ##########

 dtheta <- rnorm(length(theta))

 dphi <- transformSaturated(theta, fred, from = "theta", to = "phi",
     differential = dtheta)

 xi <- lapply(group.grp, function(i) cumulant(theta[i],
     families[[i[1]]], delta = delta[i], deriv = 1)$first)

 mydphi <- dtheta
 for (i in rev(seq(along = group.grp))) {
      j <- group.pred[i]
      if (j != 0) {
          mydphi[j] <- mydphi[j] - sum(xi[[i]] * dtheta[group.grp[[i]]])
      }
 }

 all.equal(dphi, mydphi)

 epsilon <- 1e-8
 mymydphi <- (transformSaturated(theta + epsilon * dtheta, fred,
     from = "theta", to = "phi") -
     transformSaturated(theta - epsilon * dtheta, fred,
     from = "theta", to = "phi")) / (2 * epsilon)

  all.equal(dphi, mymydphi)

 ########## phi to theta ##########

 newtheta <- transformSaturated(phi, fred, from = "phi", to = "theta")

 all.equal(theta, newtheta)

 ########## dphi to dtheta ##########

 newdtheta <- transformSaturated(phi, fred, from = "phi", to = "theta",
     differential = dphi)

 all.equal(dtheta, newdtheta)

 ########## theta to xi ##########

 xi <- transformSaturated(theta, fred, from = "theta", to = "xi")
 is.validxi(fred, xi)

 myxi <- rep(NA, length(xi))
 for (i in seq(along = group.grp)) {
     j = group.grp[[i]]
     out <- cumulant(theta[j], families[[j[1]]], delta = delta[j], deriv = 1)
     myxi[j] <- out$first
 }

 all.equal(xi, myxi)

 ########## dtheta to dxi ##########

 dxi <- transformSaturated(theta, fred, from = "theta", to = "xi",
     differential = dtheta)

 mydxi <- rep(NA, length(dxi))
 for (i in seq(along = group.grp)) {
     j = group.grp[[i]]
     out <- cumulant(theta[j], families[[j[1]]], delta = delta[j], deriv = 2)
     mydxi[j] <- as.numeric(out$second %*% dtheta[j])
 }

 all.equal(dxi, mydxi)

 ########## xi to theta ##########

 mytheta <- transformSaturated(xi, fred, from = "xi", to = "theta")
 is.same(theta, mytheta, fred, parm.type = "theta")
 ## mymyxi <- transformSaturated(mytheta, fred, from = "theta", to = "xi")
 ## all.equal(xi, mymyxi)

 ########## dxi to dtheta ##########

 mydtheta <- transformSaturated(xi, fred, from = "xi", to = "theta",
     differential = dxi)

 mymydtheta <- (
     transformSaturated(xi + epsilon * dxi, fred, from = "xi", to = "theta") -
     transformSaturated(xi - epsilon * dxi, fred, from = "xi", to = "theta")
     ) / (2 * epsilon)

 all.equal(mydtheta, mymydtheta, tol = 1e-7)
 is.same(dtheta, mydtheta, fred, parm.type = "theta")

 ########## xi to mu ##########

 mu <- transformSaturated(xi, fred, from = "xi", to = "mu")

 mymu <- rep(NA, length(mu))
 while (any(is.na(mymu))) {
     mymu.pred <- ifelse(fred$repred == 0, fred$initial,
         c(0, mymu)[fred$repred + 1])
     mymu <- xi * mymu.pred
 }
 all.equal(mu, mymu)

 ## mupred <- ifelse(fred$repred == 0, fred$initial, c(0, mu)[fred$repred + 1])
 ## mymymyxi <- mu / mupred
 ## all.equal(xi, mymymyxi)

 ########## dxi to dmu ##########

 dmu <- transformSaturated(xi, fred, from = "xi", to = "mu",
     differential = dxi)

 mydmu <- (
     transformSaturated(xi + epsilon * dxi, fred, from = "xi", to = "mu") -
     transformSaturated(xi - epsilon * dxi, fred, from = "xi", to = "mu")
     ) / (2 * epsilon)

 all.equal(dmu, mydmu)

 ########## mu to xi ##########

 mymymymyxi <- transformSaturated(mu, fred, from = "mu", to = "xi")
 all.equal(xi, mymymymyxi)

 ########## dmu to dxi ##########

 mymymymydxi <- transformSaturated(mu, fred, from = "mu", to = "xi",
     differential = dmu)
 all.equal(dxi, mymymymydxi)

 ########## theta to mu ##########

 mu.foo <- transformSaturated(theta, fred, from = "theta", to = "mu")
 all.equal(mu, mu.foo)

 ########## dtheta to dmu ##########

 dmu.foo <- transformSaturated(theta, fred, from = "theta", to = "mu",
     differential = dtheta)
 all.equal(dmu, dmu.foo)

 ########## xi to phi ##########

 phi.foo <- transformSaturated(xi, fred, from = "xi", to = "phi")
 is.same(phi, phi.foo, fred, parm.type ="phi")

 ########## dxi to dphi ##########

 dphi.foo <- transformSaturated(xi, fred, from = "xi", to = "phi",
     differential = dxi)
 is.same(dphi, dphi.foo, fred, parm.type ="phi")

 ########## phi to xi ##########

 xi.foo <- transformSaturated(phi, fred, from = "phi", to = "xi")
 all.equal(xi, xi.foo)

 ########## dphi to dxi ##########

 dxi.foo <- transformSaturated(phi, fred, from = "phi", to = "xi",
     differential = dphi)
 all.equal(dxi, dxi.foo)

 ########## phi to mu ##########

 mu.foo <- transformSaturated(phi, fred, from = "phi", to = "mu")
 all.equal(mu, mu.foo)

 ########## dphi to dmu ##########

 dmu.foo <- transformSaturated(phi, fred, from = "phi", to = "mu",
     differential = dphi)
 all.equal(dmu, dmu.foo)

 ########## mu to theta ##########

 theta.bar <- transformSaturated(mu, fred, from = "mu", to = "theta")
 is.same(theta, theta.bar, fred, parm.type ="theta")

 ########## dmu to dtheta ##########

 dtheta.bar <- transformSaturated(mu, fred, from = "mu", to = "theta",
     differential = dmu)
 is.same(dtheta, dtheta.bar, fred, parm.type ="theta")

 ########## mu to phi ##########

 phi.bar <- transformSaturated(mu, fred, from = "mu", to = "phi")
 is.same(phi, phi.bar, fred, parm.type ="phi")

 ########## dmu to dphi ##########

 dphi.bar <- transformSaturated(mu, fred, from = "mu", to = "phi",
     differential = dmu)
 is.same(dphi, dphi.bar, fred, parm.type ="phi")

 ########## valid xi ##########

 xi <- rnorm(nrow(fred$redata))
 try(validxi(fred, xi))
 nind <- nrow(fred$redata) / nlevels(fred$redata$varb)
 xi[as.character(fred$redata$varb) == "z1"] <- 1 + rexp(nind)
 try(validxi(fred, xi))
 xi[as.character(fred$redata$varb) == "p1"] <- rexp(nind)
 try(validxi(fred, xi))
 xi[as.character(fred$redata$varb) == "b1"] <- runif(nind)
 try(validxi(fred, xi))
 xi[as.character(fred$redata$varb) == "n2"] <-
     xi[as.character(fred$redata$varb) == "n1"]^2 + rchisq(nind, df = 1)
 try(validxi(fred, xi))
 xi[as.character(fred$redata$varb) == "m1"]  <- rexp(nind)
 xi[as.character(fred$redata$varb) == "m2"]  <- rexp(nind)
 xi[as.character(fred$redata$varb) == "m3"]  <- rexp(nind)
 try(validxi(fred, xi))
 thesum <- xi[as.character(fred$redata$varb) == "m1"] +
     xi[as.character(fred$redata$varb) == "m2"] +
     xi[as.character(fred$redata$varb) == "m3"]
 xi[as.character(fred$redata$varb) == "m1"]  <-
     xi[as.character(fred$redata$varb) == "m1"] / thesum
 xi[as.character(fred$redata$varb) == "m2"]  <-
     xi[as.character(fred$redata$varb) == "m2"] / thesum
 xi[as.character(fred$redata$varb) == "m3"]  <-
     xi[as.character(fred$redata$varb) == "m3"] / thesum
 is.validxi(fred, xi)


