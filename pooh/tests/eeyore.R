
 library(pooh)
 set.seed(42)

 d <- 100
 nlink <- 60

 from <- sample(1:d, nlink, replace = TRUE)
 to <- sample(1:d, nlink, replace = TRUE)

 out <- weak(from, to)

 dout <- weak(from, to, domain = 1:d)

 omega <- unlist(out)
 m <- length(omega)

 fred <- matrix(0, m, m)
 sally <- match(1:d, omega)

 for (i in seq(along = from)) {
     j <- sally[from[i]]
     k <- sally[to[i]]
     fred[j, k] <- 1
     fred[k, j] <- 1
 }
 
 fran <- fred
 repeat {
    fran.old <- fran
    fran <- fran %*% fred + fred
    fran[fran > 1] <- 1
    if (identical(fran, fran.old)) break
 }
 diag(fran) <- 1

 herman <- integer(0) 
 for (i in seq(along = out))
     herman <- c(herman, rep(i, length(out[[i]])))

 blair <- outer(herman, herman, "-")
 blair <- blair == 0
 mode(blair) <- "numeric"

 identical(fran, blair)

 fawn <- lapply(out, function(x) paste(sort(x), collapse = ":"))
 dawn <- lapply(dout, function(x) paste(sort(x), collapse = ":"))
 all(fawn %in% dawn)
 setequal(1:d, unlist(dout))

