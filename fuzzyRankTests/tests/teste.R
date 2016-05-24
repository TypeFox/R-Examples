
 library(fuzzyRankTests)

 # follow sections 2.1.2 and 2.2 of the design doc

 myfun <- function(x, mu, alternative = c("two.sided", "less", "greater")) {
     alternative <- match.arg(alternative)
     ll <- sum(x < mu)
     tt <- sum(x == mu)
     uu <- sum(x > mu)
     n <- length(x)
     if (alternative != "two.sided") {
         if (alternative == "less") {
             foo <- ll
             ll <- uu
             uu <- foo
         }
         k <- seq(0, tt)
         values <- c(0, pbinom(k, tt, 1 / 2))
         k <- seq(0, tt + 1)
         knots <- pbinom(uu + tt - k, n, 1 / 2, lower.tail = FALSE)
     } else {
         # two sided
         if (ll > uu) {
             foo <- ll
             ll <- uu
             uu <- foo
         }
         k <- seq(0, tt + 1)
         values <- pbinom(k - 1, tt, 1 / 2)
         knots <- pbinom(uu + tt - k, n, 1 / 2, lower.tail = FALSE)
         knots <- 2 * knots
         if (max(knots) > 1) {
             kk <- seq(0, tt)
             pp <- dbinom(kk, tt, 1 / 2)
             mm.low <- 2 * (ll + kk) < n
             mm.equ <- 2 * (ll + kk) == n
             mm.hig <- 2 * (ll + kk) > n
             pp.low <- pp[mm.low]
             pp.equ <- pp[mm.equ]
             pp.hig <- pp[mm.hig]
             pp.hig <- c(pp.hig, rep(0, length(pp.low) - length(pp.hig)))
             pp.hig <- rev(pp.hig)
             pp <- c(pp.low + pp.hig, pp.equ)
             values <- c(0, cumsum(pp))
             knots <- knots[seq(along = values)]
             knots[length(knots)] <- 1
         }
     }

     return(list(knots = as.numeric(knots), values = as.numeric(values)))
 }

 x <- as.double(c(-3, -2, -2, 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 5, 6, 7))
 mu <- as.double(0)

 # less
 out <- .Call("fpvsign", sum(x > mu), sum(x == mu), sum(x < mu), 1,
     PACKAGE = "fuzzyRankTests")
 mout <- myfun(x, mu, "less")
 print(out)
 print(mout)
 all.equal(out, mout)

 # greater
 out <- .Call("fpvsign", sum(x < mu), sum(x == mu), sum(x > mu), 1,
     PACKAGE = "fuzzyRankTests")
 mout <- myfun(x, mu, "great")
 print(out)
 print(mout)
 all.equal(out, mout)

 # two-tailed
 out <- .Call("fpvsign", sum(x > mu), sum(x == mu), sum(x < mu), 2,
     PACKAGE = "fuzzyRankTests")
 mout <- myfun(x, mu, "two")
 print(out)
 print(mout)
 all.equal(out, mout)

 x2 <- as.double(c(-4, -4, x))

 # less
 out <- .Call("fpvsign", sum(x2 > mu), sum(x2 == mu), sum(x2 < mu), 1,
     PACKAGE = "fuzzyRankTests")
 mout <- myfun(x2, mu, "less")
 # print(out)
 # print(mout)
 all.equal(out, mout)

 mutoo <- as.double(2)

 # greater
 out <- .Call("fpvsign", sum(x < mutoo), sum(x == mutoo), sum(x > mutoo), 1,
     PACKAGE = "fuzzyRankTests")
 mout <- myfun(x, mutoo, "great")
 # print(out)
 # print(mout)
 all.equal(out, mout)

 # test hard case

 x3 <- c(-5, -6, - 7, x2)

 out <- .Call("fpvsign", sum(x3 < mu), sum(x3 == mu), sum(x3 > mu), 2,
     PACKAGE = "fuzzyRankTests")
 mout <- myfun(x3, mu, "two")
 # print(out)
 # print(mout)
 all.equal(out, mout)

 x3 <- c(x3, 10)

 out <- .Call("fpvsign", sum(x3 < mu), sum(x3 == mu), sum(x3 > mu), 2,
     PACKAGE = "fuzzyRankTests")
 mout <- myfun(x3, mu, "two")
 # print(out)
 # print(mout)
 all.equal(out, mout)

