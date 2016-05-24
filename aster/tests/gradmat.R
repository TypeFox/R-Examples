
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

 out2 <- aster(x, root, pred, fam, modmat = out$modmat)
 summary(out2)

 out3 <- aster(x, root, pred, fam, modmat = out$modmat, type = "cond")
 summary(out3)

 foo <- match(sort(unique(site)), site)
 modmat.pred <- out$modmat[foo, , ]

 ##### case 1: eta = theta, zeta = phi

 beta.hat <- out3$coef

 modmat.pred.mat <- matrix(modmat.pred, ncol = length(beta.hat))

 theta.hat <- matrix(modmat.pred.mat %*% beta.hat, nrow = dim(modmat.pred)[1])

 nind <- dim(modmat.pred)[1]
 nnode <- dim(modmat.pred)[2]
 ncoef <- dim(modmat.pred)[3]

 aster:::setfam(fam.default())

 phi.hat <- .C("aster_theta2phi",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     theta = as.double(theta.hat),
     phi = matrix(as.double(0), nind, nnode))$phi

 my.phi.hat <- theta.hat
 my.phi.hat[ , 4] <- my.phi.hat[ , 4] - log(exp(exp(theta.hat[ , 6])) - 1)
 my.phi.hat[ , 3] <- my.phi.hat[ , 3] - log(exp(exp(theta.hat[ , 5])) - 1)
 my.phi.hat[ , 2] <- my.phi.hat[ , 2] - log(1 + exp(theta.hat[ , 4]))
 my.phi.hat[ , 1] <- my.phi.hat[ , 1] - log(1 + exp(theta.hat[ , 3]))
 my.phi.hat[ , 1] <- my.phi.hat[ , 1] - log(1 + exp(theta.hat[ , 2]))
 all.equal(my.phi.hat, phi.hat)
 
 gradmat <- 0 * modmat.pred
 storage.mode(gradmat) <- "double"

 gradmat <- .C("aster_D_beta2theta2phi",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     ncoef = as.integer(ncoef),
     pred = as.integer(pred),
     fam = as.integer(fam),
     theta = as.double(theta.hat),
     modmat = as.double(modmat.pred),
     gradmat = gradmat)$gradmat

 my.gradmat <- 0 * gradmat
 epsilon <- 1e-9
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     theta.epsilon <- matrix(modmat.pred.mat %*% beta.epsilon, nrow = nind)
     phi.epsilon <- .C("aster_theta2phi",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         theta = as.double(theta.epsilon),
         phi = matrix(as.double(0), nind, nnode))$phi
     my.gradmat[ , , k] <- (phi.epsilon - phi.hat) / epsilon
 }

 all.equal(gradmat, my.gradmat, tolerance = sqrt(epsilon))

 ##### case 2: eta = phi, zeta = theta

 beta.hat <- out2$coef

 phi.hat <- matrix(modmat.pred.mat %*% beta.hat, nrow = nind)

 theta.hat <- .C("aster_phi2theta",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     phi = as.double(phi.hat),
     theta = matrix(as.double(0), nind, nnode))$theta

 gradmat <- .C("aster_D_beta2phi2theta",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     ncoef = as.integer(ncoef),
     pred = as.integer(pred),
     fam = as.integer(fam),
     theta = as.double(theta.hat),
     modmat = as.double(modmat.pred),
     gradmat = gradmat)$gradmat

 my.gradmat <- 0 * gradmat
 epsilon <- 1e-9
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     phi.epsilon <- matrix(modmat.pred.mat %*% beta.epsilon, nrow = nind)
     theta.epsilon <- .C("aster_phi2theta",
         nind = as.integer(nind),
         nnode = as.integer(nnode),
         pred = as.integer(pred),
         fam = as.integer(fam),
         phi = as.double(phi.epsilon),
         theta = matrix(as.double(0), nind, nnode))$theta
     my.gradmat[ , , k] <- (theta.epsilon - theta.hat) / epsilon
 }

 all.equal(gradmat, my.gradmat, tolerance = sqrt(epsilon))

 ##### case 3: eta = phi, zeta = tau

 root.pred <- matrix(1, nind, nnode)

 beta.hat <- out2$coef

 beta2tau <- function(beta) {

     phi <- matrix(modmat.pred.mat %*% beta, nrow = nind)

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

 tau.hat <- beta2tau(beta.hat)

 gradmat <- .C("aster_D_beta2phi2tau",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     ncoef = as.integer(ncoef),
     pred = as.integer(pred),
     fam = as.integer(fam),
     beta = as.double(beta.hat),
     root = as.double(root.pred),
     origin = rep(as.double(0), nind * nnode),
     modmat = as.double(modmat.pred),
     gradmat = gradmat)$gradmat

 my.gradmat <- 0 * gradmat
 epsilon <- 1e-9
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     tau.epsilon <- beta2tau(beta.epsilon)
     my.gradmat[ , , k] <- (tau.epsilon - tau.hat) / epsilon
 }

 all.equal(gradmat, my.gradmat, tolerance = sqrt(epsilon))

 ##### case 4: eta = theta, zeta = tau

 beta.hat <- out3$coef

 beta2tau <- function(beta) {

     theta <- matrix(modmat.pred.mat %*% beta, nrow = nind)

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

 tau.hat <- beta2tau(beta.hat)

 gradmat <- .C("aster_D_beta2theta2tau",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     ncoef = as.integer(ncoef),
     pred = as.integer(pred),
     fam = as.integer(fam),
     beta = as.double(beta.hat),
     root = as.double(root.pred),
     modmat = as.double(modmat.pred),
     gradmat = gradmat)$gradmat

 my.gradmat <- 0 * gradmat
 for (k in 1:ncoef) {
     beta.epsilon <- beta.hat
     beta.epsilon[k] <- beta.hat[k] + epsilon
     tau.epsilon <- beta2tau(beta.epsilon)
     my.gradmat[ , , k] <- (tau.epsilon - tau.hat) / epsilon
 }

 all.equal(gradmat, my.gradmat, tolerance = sqrt(epsilon))

