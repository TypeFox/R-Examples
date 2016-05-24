
 library(ump)
 ns <- seq(10, 40, 10)
 alphas <- ps <- seq(0, 1, 0.025)
 epsilon <- 1e-9
 for (n in ns) {

     x <- seq(0, n)

     save <- array(NA, c(length(x), length(ps), length(alphas)))

     for (p in ps) {
         for (alpha in alphas) {
             px <- dbinom(x, n, p)
             phix <- umpu.binom(x, n, p, alpha)
             save[ , ps == p, alphas == alpha] <- phix
             foo <- sum(phix * px)
             bar <- sum(x * phix * px)
             if (p == 0 & alpha > 0 & alpha < 1) {
                 fred <- phix[1:2] == alpha
                 sally <- phix[- (1:2)] == 1
                 if ((! all(fred)) | (! all(sally)))
                     cat("n =", n, ": p =", p, ": alpha =", alpha,
                         ": fred =", fred, ": sally =", sally, "\n")
             } else if (p == 1 & alpha > 0 & alpha < 1) {
                 fred <- rev(phix)[1:2] == alpha
                 sally <- rev(phix)[- (1:2)] == 1
                 if ((! all(fred)) | (! all(sally)))
                     cat("n =", n, ": p =", p, ": alpha =", alpha,
                         ": fred =", fred, ": sally =", sally, "\n")
             } else {
                 if (abs(foo - alpha) > epsilon)
                     cat("n =", n, ": p =", p, ": alpha =", alpha,
                         ": foo =", foo, "\n")
                 if (abs(bar - n * p * alpha) > epsilon)
                     cat("n =", n, ": p =", p, ": alpha =", alpha,
                         ": bar =", bar,
                         ": n * p * alpha =", n * p * alpha, "\n")
             }
         }
     }

     xs <- seq(0, n)

     for (x in xs)
         for (p in ps) {
             phix <- umpu.binom(x, n, p, alphas)
             phix.too <- save[xs == x, ps == p, ]
             if (! all.equal(phix, phix.too))
                 cat("oopsie!\n")
         }

     for (x in xs)
         for (alpha in alphas) {
             phix <- umpu.binom(x, n, ps, alpha)
             phix.too <- save[xs == x, , alphas == alpha]
             if (! all.equal(phix, phix.too))
                 cat("oopsie!\n")
         }
 }

