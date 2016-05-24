
 library(ump)

 ns <- seq(10, 40, 10)
 ps <- seq(0, 1, 0.025)
 ps <- ps[0 < ps & ps < 1]

 epsilon <- 1e-9

 for (n in ns) {
     xs <- seq(0, n)
     for (p in ps) {
         for (x in xs) {
             ## cat("n =", n, ": p =", p, ": x =", x)
             out <- arpv.binom(x, n, p, plot = FALSE)
             ## cat(": alpha =", round(out$alpha, 3),
             ##     ": phi =", round(out$phi, 3), "\n")
             aout <- approx(out$alpha, out$phi)
             bout <- umpu.binom(x, n, p, aout$x, maxiter = 100)
             if (! all(abs(aout$y - bout) < epsilon))
                 cat("oopsie!\n")
         }
     }
 }

