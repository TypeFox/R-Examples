
 library(rcdd)

 set.seed(42)

 ways <- 7

 dat <- matrix(NA, nrow = 2^ways, ncol = ways)
 for (i in 1:ways)
    dat[ , i] = rep(rep(0:1, each = 2^(i - 1)), times = 2^(ways - i))

 colnames(dat) <- paste("v", 1:ways, sep = "")
 dat <- as.data.frame(dat)
 for (i in 1:ncol(dat)) dat[[i]] <- as.factor(dat[[i]])
 
 mu <- 5
 y <- rpois(nrow(dat), mu)
 dat <- cbind(dat, y = y)

 M <- model.matrix(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^3, data = dat)

 v <- M
 linearity <- y > 0

   stopifnot(is.numeric(v))
   stopifnot(all(is.finite(v)))
   stopifnot(is.matrix(v))
   if (! missing(linearity)) {
       stopifnot(is.logical(linearity))
       stopifnot(length(linearity) == nrow(v))
   } else {
       linearity <- rep(FALSE, nrow(v))
   }

   v <- d2q(v)
   lresult <- rep(TRUE, nrow(v))

       vresult <- v[lresult, , drop = FALSE]
       w <- apply(vresult, 2, qsum)
       wminus <- qmq(rep("0", length(w)), w)

       hrep <- rbind(wminus, vresult)
       fred <- c(1, rep(0, nrow(vresult)))
       hrep <- cbind(as.character(fred), hrep)
       fred <- c(0, as.numeric(linearity[lresult]))
       hrep <- cbind(as.character(fred), hrep)
       dimnames(hrep) <- NULL

 out <- lpcdd(hrep, w, minimize = FALSE)
 print(out)

 ##### check gradient of Lagrangian function zero
 b.augmented <- hrep[ , 2]
 v.augmented <- hrep[ , - c(1, 2)]
 blurfle <- out$dual.solution
 all(qpq(w, qmatmult(rbind(blurfle), v.augmented)) == "0")

 ##### check primal feasibility (trivial here since solution is zero, but)
 blurfle <- out$primal.solution
 foo <- qpq(b.augmented, qmatmult(v.augmented, cbind(blurfle)))
 all(qsign(foo) >= 0)
 
 ##### check dual feasibility
 blurfle <- qsign(out$dual.solution)
 linearity.augmented <- c(FALSE, linearity)
 length(linearity.augmented) == nrow(hrep)
 all(blurfle[! linearity.augmented] >= 0)

 ##### check complementary slackness
 all(qsign(foo) * qsign(blurfle) == 0)

 ##### now redo with ordinary computer arithmetic
 aout <- lpcdd(q2d(hrep), q2d(w), minimize = FALSE)
 names(aout)
 names(out)
 for (i in 2:length(out))
     print(all.equal(aout[[i]], q2d(out[[i]])))

