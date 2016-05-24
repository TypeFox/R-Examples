
 library(aster)

 set.seed(42)

 pred <- c(0, 1, 1, 2)
 fam <- c(1, 1, 3, 2)
 famlist <- fam.default()

 nnode <- length(pred)
 nind <- 7

 x <- rbind(c(1, 1, 2, 4),
            c(1, 1, 3, 0),
            c(1, 1, 5, 2),
            c(1, 0, 1, 0),
            c(1, 0, 4, 0),
            c(0, 0, 0, 0),
            c(0, 0, 0, 0))

 nrow(x) == nind
 ncol(x) == nnode

 root <- matrix(1, nind, nnode)

 theta <- rnorm(nind * nnode, 0, .33)
 theta <- matrix(theta, nind, nnode)
 storage.mode(theta) <- "double"

 aster:::setfam(fam.default())

 out <- .C("aster_theta2phi",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     theta = theta,
     phi = matrix(as.double(0), nind, nnode))

 my.phi <- theta
 # poisson cumulant function
 my.phi[ , 2] <- my.phi[ , 2] - exp(theta[ , 4])
 # non-zero poisson cumulant function
 my.phi[ , 1] <- my.phi[ , 1] - log(exp(exp(theta[ , 3])) - 1)
 # bernoulli cumulant function
 my.phi[ , 1] <- my.phi[ , 1] - log(1 + exp(theta[ , 2]))

 all.equal(out$phi, my.phi)

 tout <- .C("aster_phi2theta",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     phi = out$phi,
     theta = matrix(as.double(0), nind, nnode))

 all.equal(tout$theta, theta)

 storage.mode(x) <- "double"
 storage.mode(root) <- "double"

 cout <- .C("aster_theta2ctau",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     theta = theta,
     ctau = matrix(as.double(0), nind, nnode))

 my.ctau <- theta
 # poisson mean function
 my.ctau[ , 4] <- exp(theta[ , 4])
 # non-zero poisson mean function
 foo <- exp(theta[ , 3])
 my.ctau[ , 3] <- foo / (1 - exp(- foo))
 # bernoulli cumulant function
 my.ctau[ , 2] <- 1 / (1 + exp(- theta[ , 2]))
 # bernoulli cumulant function
 my.ctau[ , 1] <- 1 * 1 / (1 + exp(- theta[ , 1]))

 all.equal(cout$ctau, my.ctau)

 xout <- .C("aster_xpred",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     x = x,
     root = root,
     xpred = matrix(as.double(0), nind, nnode))

 my.xpred <- x
 my.xpred[ , 1] <- root[ , 1]
 my.xpred[ , 2] <- x[ , 1]
 my.xpred[ , 3] <- x[ , 1]
 my.xpred[ , 4] <- x[ , 2]

 all.equal(xout$xpred, my.xpred)

 tout <- .C("aster_ctau2tau",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     root = root,
     ctau = cout$ctau,
     tau = matrix(as.double(0), nind, nnode))

 my.tau <- my.ctau
 my.tau[ , 1] <- my.xpred[ , 1] * my.ctau[ , 1]
 my.tau[ , 2] <- my.tau[ , 1] * my.ctau[ , 2]
 my.tau[ , 3] <- my.tau[ , 1] * my.ctau[ , 3]
 my.tau[ , 4] <- my.tau[ , 2] * my.ctau[ , 4]

 all.equal(tout$tau, my.tau)

 ##########

 wout <- .C("aster_theta2whatsis",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     deriv = as.integer(0),
     theta = as.double(theta),
     result = matrix(as.double(0), nind, nnode))

 my.what <- NaN * wout$result
 for (i in 1:nind)
     for (j in 1:nnode) {
         ifam <- fam[j]
         my.what[i, j] <- famfun(famlist[[ifam]], 0, theta[i, j])
     }

 all.equal(wout$result, my.what)

 ##########

 aster:::setfam(fam.default())

 wout <- .C("aster_theta2whatsis",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     deriv = as.integer(1),
     theta = as.double(theta),
     result = matrix(as.double(0), nind, nnode))

 my.what <- NaN * wout$result
 for (i in 1:nind)
     for (j in 1:nnode) {
         ifam <- fam[j]
         my.what[i, j] <- famfun(famlist[[ifam]], 1, theta[i, j])
     }
 all.equal(wout$result, my.what)

 ##########

 aster:::setfam(fam.default())

 wout <- .C("aster_theta2whatsis",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     deriv = as.integer(2),
     theta = as.double(theta),
     result = matrix(as.double(0), nind, nnode))

 my.what <- NaN * wout$result
 for (i in 1:nind)
     for (j in 1:nnode) {
         ifam <- fam[j]
         my.what[i, j] <- famfun(famlist[[ifam]], 2, theta[i, j])
     }

 all.equal(wout$result, my.what)

 ##########

 aster:::setfam(fam.default())

 vout <- .C("aster_tt2var",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     x = as.double(x),
     root = as.double(root),
     theta = as.double(theta),
     tau = as.double(tout$tau),
     var = matrix(as.double(0), nind, nnode))

 my.var <- NaN * vout$var
 for (i in 1:nind)
     for (j in 1:nnode) {
         ifam <- fam[j]
         thefam <- famlist[[ifam]]
         k <- pred[j]
         if (k > 0) {
             my.var[i, j] <- famfun(thefam, 2, theta[i, j]) * tout$tau[i, k] +
                 famfun(thefam, 1, theta[i, j])^2 * my.var[i, k]
         } else {
             my.var[i, j] <- famfun(thefam, 2, theta[i, j]) * xout$xpred[i, j]
         }
     }
 all.equal(vout$var, my.var)

