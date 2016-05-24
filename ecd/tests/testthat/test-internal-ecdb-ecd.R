
library(xts)

context("Test on Internal ECDB and ecdattr")

db <- ecdb("internal")

get_cnt <- function(alpha, gamma) {
    dim(read(db, alpha=alpha, gamma=gamma))[[1]]
}
rect_cnt <- function(L, step) {
    rng <- seq(-L, L, step)
    get_cnt(rng, rng)
}
                
eps <- 0.001 # default tolerance of error for real number

test_that("test ecd(0,0) is cusp",{
    a <- read(db,0,0)
  	expect_true(a$cusp > 0)
})

test_that("test ecd(2,-3) is cusp",{
    a <- read(db,2,-3)
  	expect_true(a$cusp > 0)
})

test_that("test count for 10x10 rectangle",{
    L <- 10
    step <- 0.25
    c1 <- rect_cnt(L, step)
    c2 <- (L/step*2+1)^2
  	expect_true(c1 == c2)
})

test_that("test history",{
    h <- "bootstrap-rect-internal"
  	expect_true(h %in% history(db))
})

# ---------------------------------------------
test_that("test ecdattr",{
    a <- 1
	r <- c(1, 2 ,3)
	p <- ecdattr.pairs(a, r)
	p2 <- parallel::mclapply(p, ecdattr.enrich)
	d2 <- parallel::mclapply(r, function(r) ecd(a, r))
	err_kt <- Map(function(p,d) abs(p@attr$kurtosis/ecd.kurtosis(d)-1), p2, d2)
	err_di <- Map(function(p,d) abs(p@attr$discr/discr(d)-1), p2, d2)
	err <- sum(simplify2array(err_kt)) + sum(simplify2array(err_di))
	expect_true(err < eps)
})

