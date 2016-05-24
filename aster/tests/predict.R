
 library(aster)

 set.seed(42)

 nind <- 25

 vars <- c("l2", "l3", "f2", "f3", "h2", "h3")
 pred <- c(0, 1, 1, 2, 3, 4)
 fam <- c(1, 1, 1, 1, 3, 3)
 length(pred) == length(fam)
 nnode <- length(pred)

 theta <- matrix(0, nind, nnode)
 root <- matrix(1, nind, nnode)
 x <- raster(theta, pred, fam, root)
 dimnames(x) <- list(NULL, vars)

 data <- as.data.frame(x)
 site <- factor(sample(LETTERS[1:4], nind, replace = TRUE))
 foo <- rnorm(nind)
 data <- data.frame(x, site = site, foo = foo, root = 1)

 redata <- reshape(data, varying = list(vars),
     direction = "long", timevar = "varb", times = as.factor(vars),
     v.names = "resp")

 out <- aster(resp ~ foo + site + varb, pred, fam, varb, id, root,
     data = redata)
 summary(out, show.graph = TRUE)

 ##### redo with aster.default and predict.aster

 out2 <- aster(x, root, pred, fam, modmat = out$modmat)
 summary(out2)

 foo <- match(sort(unique(site)), site)
 modmat.pred <- out$modmat[foo, , ]
 origin.pred <- out$origin[foo, ]

 predict(out2, modmat = modmat.pred, parm.type = "canon")

 ##### case 1: model = "unco", obj = "unco", parm = "cano" ####

 fred <- predict(out2, modmat = modmat.pred, parm.type = "canon",
     se.fit = TRUE)

 all.equal(fred$se.fit, sqrt(diag(fred$gradient %*% solve(out2$fisher) %*%
     t(fred$gradient))))

 sally <- matrix(modmat.pred, ncol = length(out2$coef))

 all.equal(fred$gradient, sally)

 all.equal(fred$fit, as.numeric(origin.pred) + as.numeric(sally %*% out$coef))

 ##### case 1a: same but with amat

 node.names <- dimnames(out$modmat)[[2]]
 site.names <- levels(site)
 amat <- array(0, c(dim(modmat.pred)[1:2], length(site.names)))
 for (i in seq(along = site.names))
     amat[i, grep("h", node.names), i] <- 1

 alfie <- predict(out2, modmat = modmat.pred, parm.type = "canon",
     se.fit = TRUE, amat = amat)

 amatmat <- matrix(amat, ncol = dim(amat)[3])

 all.equal(alfie$fit, as.numeric(t(amatmat) %*% fred$fit))

 all.equal(alfie$gradient, t(amatmat) %*% fred$gradient)

 all.equal(alfie$se.fit, sqrt(diag(alfie$gradient %*% solve(out2$fisher) %*%
     t(alfie$gradient))))

 ##### case 2: model = "cond", obj = "cond", parm = "cano" ####
 ##### no test -- same code as case 1

 ##### case 3: model = "unco", obj = "cond", parm = "cano" ####

 out3 <- aster(x, root, pred, fam, modmat = out$modmat, type = "cond")
 summary(out3)

 fred <- predict(out3, modmat = modmat.pred, parm.type = "canon",
     se.fit = TRUE)

 nind <- dim(modmat.pred)[1]
 nnode <- dim(modmat.pred)[2]
 ncoef <- dim(modmat.pred)[3]

 aster:::setfam(fam.default())

 beta.hat <- out3$coef
 theta.hat <- as.numeric(sally %*% beta.hat)
 phi.hat <- .C("aster_theta2phi",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     theta = as.double(theta.hat),
     phi = double(nind * nnode))$phi

 all.equal(fred$fit, phi.hat)

 all.equal(fred$se.fit, sqrt(diag(fred$gradient %*% solve(out3$fisher) %*%
     t(fred$gradient))))

 my.gradient <- 0 * fred$gradient
 epsilon <- 1e-9
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     theta.epsilon <- as.numeric(sally %*% beta.epsilon)
     phi.epsilon <- .C("aster_theta2phi",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         theta = as.double(theta.epsilon),
         phi = double(nind * nnode))$phi
     my.gradient[ , k] <- (phi.epsilon - phi.hat) / epsilon
 }

 all.equal(fred$gradient, my.gradient, tolerance = sqrt(epsilon))

 alfie <- predict(out3, modmat = modmat.pred, parm.type = "canon",
     se.fit = TRUE, amat = amat)

 all.equal(alfie$fit, as.numeric(t(amatmat) %*% fred$fit))

 all.equal(alfie$gradient, t(amatmat) %*% fred$gradient)

 all.equal(alfie$se.fit, sqrt(diag(alfie$gradient %*% solve(out3$fisher) %*%
     t(alfie$gradient))))

 ##### case 4: model = "cond", obj = "unco", parm = "cano" ####

 fred <- predict(out2, modmat = modmat.pred, parm.type = "canon",
     model.type = "cond", se.fit = TRUE)

 aster:::setfam(fam.default())

 beta.hat <- out2$coef
 phi.hat <- as.numeric(origin.pred) + as.numeric(sally %*% beta.hat)
 theta.hat <- .C("aster_phi2theta",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     phi = as.double(phi.hat),
     theta = double(nind * nnode))$theta

 all.equal(fred$fit, theta.hat)

 all.equal(fred$se.fit, sqrt(diag(fred$gradient %*% solve(out2$fisher) %*%
     t(fred$gradient))))

 my.gradient <- 0 * fred$gradient
 epsilon <- 1e-9
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     phi.epsilon <- as.numeric(origin.pred) + as.numeric(sally %*% beta.epsilon)
     theta.epsilon <- .C("aster_phi2theta",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         phi = as.double(phi.epsilon),
         theta = double(nind * nnode))$theta
     my.gradient[ , k] <- (theta.epsilon - theta.hat) / epsilon
 }

 all.equal(fred$gradient, my.gradient, tolerance = sqrt(epsilon))

 alfie <- predict(out2, modmat = modmat.pred, parm.type = "canon",
     model.type = "cond", se.fit = TRUE, amat = amat)

 all.equal(alfie$fit, as.numeric(t(amatmat) %*% fred$fit))

 all.equal(alfie$gradient, t(amatmat) %*% fred$gradient)

 all.equal(alfie$se.fit, sqrt(diag(alfie$gradient %*% solve(out2$fisher) %*%
     t(alfie$gradient))))

 ##### case 5: model = "cond", obj = "cond", parm = "mean" ####

 root.pred <- matrix(1, nind, nnode)

 fred <- predict(out3, modmat = modmat.pred, parm.type = "mean",
     model.type = "cond", root = root.pred, x = root.pred)

 aster:::setfam(fam.default())

 beta.hat <- out3$coef
 theta.hat <- as.numeric(sally %*% beta.hat)
 xi.hat <- .C("aster_theta2ctau",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     theta = as.double(theta.hat),
     ctau = double(nind * nnode))$ctau

 all.equal(fred, xi.hat)

 fred <- predict(out3, modmat = modmat.pred, parm.type = "mean",
     model.type = "cond", root = root.pred, x = root.pred, se.fit = TRUE)

 all.equal(fred$fit, xi.hat)

 all.equal(fred$se.fit, sqrt(diag(fred$gradient %*% solve(out3$fisher) %*%
     t(fred$gradient))))

 aster:::setfam(fam.default())

 my.gradient <- 0 * fred$gradient
 epsilon <- 1e-9
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     theta.epsilon <- as.numeric(sally %*% beta.epsilon)
     xi.epsilon <- .C("aster_theta2ctau",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         theta = as.double(theta.epsilon),
         ctau = double(nind * nnode))$ctau
     my.gradient[ , k] <- (xi.epsilon - xi.hat) / epsilon
 }

 all.equal(fred$gradient, my.gradient, tolerance = sqrt(epsilon))

 ##### case 6: model = "unco", obj = "unco", parm = "mean" ####

 fred <- predict(out2, modmat = modmat.pred, parm.type = "mean",
     root = root.pred)

 beta.hat <- out2$coef

 beta2tau <- function(beta) {

     phi <- origin.pred + matrix(sally %*% beta, nrow = nind)

     theta <- .C("aster_phi2theta",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         phi = as.double(phi),
         theta = matrix(as.double(0), nind, nnode))$theta

     ctau <- .C("aster_theta2ctau",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         theta = as.double(theta),
         ctau = double(nind * nnode))$ctau

     tau <- .C("aster_ctau2tau",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         root = as.double(root.pred),
         ctau = as.double(ctau),
         tau = double(nind * nnode))$tau

     return(tau)
 }

 aster:::setfam(fam.default())

 tau.hat <- beta2tau(beta.hat)

 all.equal(fred, tau.hat)

 fred <- predict(out2, modmat = modmat.pred, parm.type = "mean",
     root = root.pred, se.fit = TRUE)

 all.equal(fred$fit, tau.hat)

 all.equal(fred$se.fit, sqrt(diag(fred$gradient %*% solve(out2$fisher) %*%
     t(fred$gradient))))

 aster:::setfam(fam.default())

 my.gradient <- 0 * fred$gradient
 for (k in 1:length(beta.hat)) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     tau.epsilon <- beta2tau(beta.epsilon)
     my.gradient[ , k] <- (tau.epsilon - tau.hat) / epsilon
 }

 all.equal(fred$gradient, my.gradient, tolerance = sqrt(epsilon))

 ##### case 7: model = "cond", obj = "unco", parm = "mean" ####

 fred <- predict(out2, modmat = modmat.pred, parm.type = "mean",
     model.type = "cond", root = root.pred, x = root.pred)

 beta.hat <- out2$coef

 beta2xi <- function(beta) {

     phi <- origin.pred + matrix(sally %*% beta, nrow = nind)

     theta <- .C("aster_phi2theta",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         phi = as.double(phi),
         theta = matrix(as.double(0), nind, nnode))$theta

     ctau <- .C("aster_theta2ctau",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         theta = as.double(theta),
         ctau = double(nind * nnode))$ctau

     return(ctau)
 }

 aster:::setfam(fam.default())

 xi.hat <- beta2xi(beta.hat)

 all.equal(fred, xi.hat)

 fred <- predict(out2, modmat = modmat.pred, parm.type = "mean",
     model.type = "cond", root = root.pred, x = root.pred, se.fit = TRUE)

 all.equal(fred$fit, xi.hat)

 all.equal(fred$se.fit, sqrt(diag(fred$gradient %*% solve(out2$fisher) %*%
     t(fred$gradient))))

 aster:::setfam(fam.default())

 my.gradient <- 0 * fred$gradient
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     xi.epsilon <- beta2xi(beta.epsilon)
     my.gradient[ , k] <- (xi.epsilon - xi.hat) / epsilon
 }

 all.equal(fred$gradient, my.gradient, tolerance = sqrt(epsilon))

 ##### case 8: model = "unco", obj = "cond", parm = "mean" ####

 fred <- predict(out3, modmat = modmat.pred, root = root.pred)

 beta.hat <- out3$coef

 beta2tau <- function(beta) {

     theta <- matrix(sally %*% beta, nrow = nind)

     ctau <- .C("aster_theta2ctau",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         theta = as.double(theta),
         ctau = double(nind * nnode))$ctau

     tau <- .C("aster_ctau2tau",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         root = as.double(root.pred),
         ctau = as.double(ctau),
         tau = double(nind * nnode))$tau

     return(tau)
 }

 aster:::setfam(fam.default())

 tau.hat <- beta2tau(beta.hat)

 all.equal(fred, tau.hat)

 fred <- predict(out3, modmat = modmat.pred, root = root.pred, se.fit = TRUE)

 all.equal(fred$fit, tau.hat)

 all.equal(fred$se.fit, sqrt(diag(fred$gradient %*% solve(out3$fisher) %*%
     t(fred$gradient))))

 aster:::setfam(fam.default())

 my.gradient <- 0 * fred$gradient
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     tau.epsilon <- beta2tau(beta.epsilon)
     my.gradient[ , k] <- (tau.epsilon - tau.hat) / epsilon
 }

 all.equal(fred$gradient, my.gradient, tolerance = sqrt(epsilon))

 ##### HOORAY !!!!! ##### That's it for aster.predict #####
 ##### now for aster.predict.formula #####

 ##### case 1: newdata missing

 predict(out)

 newdata <- data.frame(site = factor(LETTERS[1:4]))
 for (v in vars)
 newdata[[v]] <- 1
 newdata$root <- 1
 newdata$foo <- modmat.pred[ , "l2", "foo"]

 renewdata <- reshape(newdata, varying = list(vars),
     direction = "long", timevar = "varb", times = as.factor(vars),
     v.names = "resp")

 louise <- predict(out, newdata = renewdata, varvar = varb, idvar = id,
     root = root, se.fit = TRUE)

 all.equal(louise$modmat, modmat.pred)

 fred <- predict(out2, modmat = modmat.pred, root = root.pred, se.fit = TRUE)

 all.equal(louise$fit, fred$fit)

 all.equal(louise$se.fit, fred$se.fit)

 ##### test for global variables #####

 saves <- c("out", "renewdata", "out2", "modmat.pred", "root.pred", "louise",
     "fred")
 blurfle <- ls()
 blurfle <- ls()
 rm(list = blurfle[! is.element(blurfle, saves)])
 ls()

 louise.too <- predict(out, newdata = renewdata, varvar = varb, idvar = id,
     root = root, se.fit = TRUE)
 identical(louise, louise.too)

 fred.too <- predict(out2, modmat = modmat.pred, root = root.pred,
     se.fit = TRUE)
 identical(fred, fred.too)

 ##### test of newcoef #####

 fake <- out2
 beta.new <- fake$coefficients + rnorm(length(fake$coefficients)) * 0.1
 fake$coefficients <- beta.new
 fred.fake <- predict(fake, modmat = modmat.pred, root = root.pred,
     se.fit = TRUE)
 fred.new <- predict(out2, modmat = modmat.pred, root = root.pred,
     se.fit = TRUE, newcoef = beta.new)
 identical(fred.fake, fred.new)


