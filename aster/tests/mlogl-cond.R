
 library(aster)

 options(digits=4) # avoid rounding differences

 set.seed(42)

 nind <- 25
 nnode <- 5
 ncoef <- nnode + 1

 famlist <- fam.default()
 fam <- c(1, 1, 2, 3, 3)

 pred <- c(0, 1, 1, 2, 3)
 print(pred)

 modmat <- array(0, c(nind, nnode, ncoef))
 modmat[ , , 1] <- 1
 for (i in 2:nnode)
     modmat[ , i, i] <- 1
 modmat[ , , ncoef] <- rnorm(nind * nnode)

 beta <- rnorm(ncoef) / 10

 theta <- matrix(modmat, ncol = ncoef) %*% beta
 theta <- matrix(theta, ncol = nnode)

 root <- sample(1:3, nind * nnode, replace = TRUE)
 root <- matrix(root, nind, nnode)

 x <- raster(theta, pred, fam, root)
 
 out <- mlogl(beta, pred, fam, x, root, modmat, deriv = 2,
     type = "conditional")
 print(out)

 my.value <- 0
 for (j in 1:nnode) {
     ifam <- fam[j] 
     k <- pred[j]
     if (k > 0)
         xpred <- x[ , k]
     else
         xpred <- root[ , j]
     for (i in 1:nind)
         my.value <- my.value -
             sum(x[i, j] * theta[i, j] -
             xpred[i] * famfun(famlist[[ifam]], 0, theta[i, j]))
 }
 all.equal(out$value, my.value)

 my.grad <- NaN * out$gradient
 epsilon <- 1e-9
 for (i in 1:ncoef) {
     beta.eps <- beta
     beta.eps[i] <- beta[i] + epsilon
     out.eps <- mlogl(beta.eps, pred, fam, x, root, modmat, deriv = 0,
         type = "conditional")
     my.grad[i] <- (out.eps$value - out$value) / epsilon
 }

 all.equal(out$gradient, my.grad, tolerance = sqrt(epsilon))

 my.hess <- matrix(NaN, ncoef, ncoef)
 for (i in 1:ncoef) {
     beta.eps <- beta
     beta.eps[i] <- beta[i] + epsilon
     out.eps <- mlogl(beta.eps, pred, fam, x, root, modmat, deriv = 1,
         type = "conditional")
     my.hess[ , i] <- (out.eps$gradient - out$gradient) / epsilon
 }

 all.equal(out$hessian, my.hess, tolerance = sqrt(epsilon))

 ##########

 objfun <- function(beta) {
     out <- mlogl(beta, pred, fam, x, root, modmat, deriv = 1,
         type = "conditional")
     result <- out$value
     attr(result, "gradient") <- out$gradient
     return(result)
 }
 out <- nlm(objfun, beta, fscale = nind)
 print(out)

 ##########

 objfun <- function(beta) {
     out <- mlogl(beta, pred, fam, x, root, modmat, deriv = 2,
         type = "conditional")
     result <- out$value
     attr(result, "gradient") <- out$gradient
     attr(result, "hessian") <- out$hessian
     return(result)
 }
 out <- nlm(objfun, beta, fscale = nind)
 print(out)

 ########## expected Fisher information ##########

 aster:::setfam(fam.default())

 fout <- .C("aster_fish_cond",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     ncoef = as.integer(ncoef),
     pred = as.integer(pred),
     fam = as.integer(fam),
     beta = as.double(out$estimate),
     root = as.double(root),
     x = as.double(x),
     modmat = as.double(modmat),
     fish = matrix(as.double(0), ncoef, ncoef))

 mout <- mlogl(out$estimate, pred, fam, x, root, modmat,
     deriv = 2, type = "conditional")

 aster:::setfam(fam.default())

 theta <- .C("aster_mat_vec_mult",
    nrow = as.integer(nind * nnode),
    ncol = as.integer(ncoef),
    a = as.double(modmat),
    b = as.double(out$estimate),
    c = matrix(as.double(0), nind, nnode))$c
 ctau <- .C("aster_theta2ctau",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     theta = as.double(theta),
     ctau = matrix(as.double(0), nind, nnode))$ctau
 tau <- .C("aster_ctau2tau",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     root = as.double(root),
     ctau = as.double(ctau),
     tau = matrix(as.double(0), nind, nnode))$tau
 psiDoublePrime <- .C("aster_theta2whatsis",
     nind = as.integer(nind),
     nnode = as.integer(nnode),
     pred = as.integer(pred),
     fam = as.integer(fam),
     deriv = as.integer(2),
     theta = as.double(theta),
     result = matrix(as.double(0), nind, nnode))$result

 my.hessian <- theta * NaN
 my.fish <- theta * NaN

 for (i in 1:nind) {
     for (j in 1:nnode) {
         k <- pred[j]
         if (k > 0) {
             my.hessian[i, j] <- x[i, k] * psiDoublePrime[i, j]
             my.fish[i, j] <- tau[i, k] * psiDoublePrime[i, j]
         } else {
             my.hessian[i, j] <- root[i, j] * psiDoublePrime[i, j]
             my.fish[i, j] <- root[i, j] * psiDoublePrime[i, j]
         }
     }
 }

 modmat <- matrix(as.double(modmat), ncol = ncoef)
 my.hessian <- as.numeric(my.hessian)
 my.fish <- as.numeric(my.fish)
 my.hessian <- t(modmat) %*% diag(my.hessian) %*% modmat
 my.fish <- t(modmat) %*% diag(my.fish) %*% modmat

 all.equal(my.hessian, mout$hessian)
 all.equal(my.fish, fout$fish)

