
 myfun <- function(x, alternative = c("two.sided", "less", "greater")) {
     alternative <- match.arg(alternative)
     absx <- abs(x)
     ntiezero <- sum(absx == 0)
     foo <- split(x, absx)
     foo <- foo[names(foo) != "0"]
     foop <- sapply(foo, function(bar) sum(bar > 0))
     foon <- sapply(foo, function(bar) sum(bar < 0))
     mtie <- as.numeric(foon[foop > 0 & foon > 0])
     ntie <- as.numeric(foop[foop > 0 & foon > 0])
     plist <- list()
     if (ntiezero > 0) {
        i <- seq(0, ntiezero * (ntiezero + 1) / 2)
        plist[[length(plist) + 1]] <- dsignrank(i, ntiezero)
     }
     for (j in 1:length(mtie)) {
        i <- seq(0, mtie[j] * ntie[j])
        plist[[length(plist) + 1]] <- dwilcox(i, mtie[j], ntie[j])
     }
     ptie <- 1
     for (j in 1:length(plist))
        ptie <- convolve(plist[[j]], rev(ptie), type = "o")
     w <- outer(x, x, "+")
     w <- w[lower.tri(w, diag = TRUE)]
     ll <- sum(w < 0)
     tt <- sum(w == 0)
     uu <- sum(w > 0)
     if (tt != sum(sapply(plist, length) - 1))
         stop("OOPS!")
     if (alternative == "greater") {
         return(list(knots = psignrank(seq(ll - 1, ll + tt), length(x)),
             values = cumsum(c(0, ptie))))
     }
     if (alternative == "less") {
         return(list(knots = psignrank(seq(uu - 1, uu + tt), length(x)),
             values = cumsum(c(0, ptie))))
     }
     if (alternative == "two.sided") {
         wmin <- min(ll, uu)
         wmax <- wmin + tt
         w <- seq(wmin, wmax)
         pknot <- dsignrank(w, length(x))
         N <- length(x) * (length(x) + 1) / 2
         w <- pmin(w, N - w)
         foo <- split(ptie, w)
         foo <- sapply(foo, sum)
         foo <- cumsum(c(0, foo))
         bar <- pmin(1, 2 * psignrank(seq(min(w) - 1, max(w)), length(x)))
         return(list(knots = as.numeric(bar), values = as.numeric(foo)))
     }
 }

 library(fuzzyRankTests)

 x <- as.double(c(-3, -2, -2, 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 5, 6, 7))
 mu <- as.double(0)
 tol <- sqrt(.Machine$double.eps)

 out <- .Call("fpvsignrank", x, mu, "less", tol, PACKAGE = "fuzzyRankTests")
 print(out)
 all.equal(out, myfun(x, "less"))

 x2 <- as.double(c(-4, -4, x))

 out <- .Call("fpvsignrank", x2, mu, "less", tol, PACKAGE = "fuzzyRankTests")
 print(out)
 all.equal(out, myfun(x2, "less"))

 mutoo <- as.double(2)

 out <- .Call("fpvsignrank", x, mutoo, "greater", tol,
     PACKAGE = "fuzzyRankTests")
 print(out)
 all.equal(out, myfun(x - mutoo, "great"))

 ##### check for error found in devel ####
 try(.Call("fpvsignrank", x, mutoo, "two", tol, PACKAGE = "fuzzyRankTests"))

 out <- .Call("fpvsignrank", x, mutoo, "two.sided", tol,
     PACKAGE = "fuzzyRankTests")
 print(out)
 all.equal(out, myfun(x - mutoo, "two"))

