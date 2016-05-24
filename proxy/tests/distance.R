
### ceeboo 2007

## todo: special values NaN, NA, and Inf

library("proxy")

set.seed(20070630)

##

x <- matrix(runif(20),5,4)
inherits(x, "matrix")
rownames(x) <- LETTERS[1:5]
y <- x
y[2,] <- x[1,] <- 0

x
y

## dist extensions

for (i in 1:9) {
    cat("\nTesting #",i,"\n\n",sep="")
    print(.Call("R_dists", x, NULL, FALSE, i, NULL))
    print(.Call("R_dists", x, x, FALSE, i, NULL))
    print(.Call("R_dists", x, y, FALSE, i, NULL))
    print(.Call("R_dists", x, y, TRUE, i, NULL))
}

## again but via user interfaces

r <- .Call("R_minkowski_dist", x, NULL, FALSE, NULL)
all.equal(c(r), c(stats::dist(x, method = "minkowski", p = 1)))
r
.Call("R_minkowski_dist", x, x, FALSE, NULL)
.Call("R_minkowski_dist", x, y, FALSE, NULL)
.Call("R_minkowski_dist", x, y, TRUE, NULL)

dfun <- paste("R",c("euclidean", "maximum", "manhattan", "canberra", "binary", "matching", "fuzzy", "mutual"),"dist", sep = "_")

for (f in dfun) {
    cat("\nTesting ",f,"\n\n",sep="")
    r <- do.call(".Call", list(f, x, NULL, FALSE))
    s <- try(stats::dist(x, method = gsub("R_|_dist", "", f)))
    if (!inherits(s, "try-error"))
        print(all.equal(c(r), c(s)))
    print(r)
    print(do.call(".Call", list(f, x, x, FALSE)))
    print(do.call(".Call", list(f, x, y, FALSE)))
    print(do.call(".Call", list(f, x, y, TRUE)))
}

## optimized

.Call("R_ejaccard", x, NULL, FALSE)
.Call("R_ejaccard", x, x, FALSE)
.Call("R_ejaccard", x, y, FALSE)
.Call("R_ejaccard", x, y, TRUE)

.Call("R_cosine", x, NULL, FALSE)
.Call("R_cosine", x, x, FALSE)
.Call("R_cosine", x, y, FALSE)
.Call("R_cosine", x, y, TRUE)

x <- matrix(x > 0.5, 5,4)
y <- matrix(y > 0.5, 5,4)

.Call("R_bjaccard", x, NULL, FALSE)
.Call("R_bjaccard", x, x, FALSE)
.Call("R_bjaccard", x, y, FALSE)
.Call("R_bjaccard", x, y, TRUE)

###
