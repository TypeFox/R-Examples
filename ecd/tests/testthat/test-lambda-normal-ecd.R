
context("Test on Lambda=1 Normal Distribution and BS Model")

eps <- 0.001 # default tolerance of error for real number


# ------------------------------------------------------
# tests for lambda=1
ld <- ecld(lambda=1, sigma=ecd.mp1*0.5, with.RN=TRUE)
d <- ld@ecd
d_RN <- ld@ecd_RN
stopifnot(ld@lambda==1)

x <- c(-0.1, 0, 0.2)

test_that("test CDF of normal",{
    c1 <- ecd.cdf(d_RN, x)
    c2 <- ecld.cdf(ld, x-ld@mu_D)
    expect_true( max(abs(c2/c1-1)) < eps )
})

test_that("test M(1) of normal",{
    m1 <- ecld.mgf(ld)
    m2 <- ecd.imgf(d)
    expect_true( abs(m2/m1-1) < eps )
})

test_that("test mu_D of normal",{
    s <- ld@sigma
    mu <- -s^2/4 
    mu2 <- ecld.mu_D(ld)
    expect_true( abs(mu/ld@mu_D-1) + abs(mu2/mu-1) < eps )
})

# ------------------------------------------------------

test_that("test IMGF of BS model",{
    M <- exp(-ecld.mu_D(ld))
    IMc <- ecld.imgf(ld, x) # RN=TRUE by default
    IMp <- ecld.imgf(ld, x, otype="p")
    
    M2 <- ecd.imgf(d)
    IMc2 <- ecd.imgf(d_RN, x)
    IMp2 <- 1 - IMc2
    
    expect_true( abs(M2/M-1) + max(abs(IMc2/IMc-1)) + max(abs(IMp2/IMp-1)) < eps )
})

test_that("test Lc and Lp of BS model",{
    Lc <- ecd.ogf(d_RN, x, otype="c")
    Lp <- ecd.ogf(d_RN, x, otype="p")
    
    Lc2 <- ecld.ogf(ld, x, otype="c")
    Lp2 <- ecld.ogf(ld, x, otype="p")
    
    expect_true( max(abs(Lc2/Lc-1)) + max(abs(Lp2/Lp-1)) < eps )
})

test_that("test Lc,p with op_O of BS model",{
    Lc <- ecld.op_O(ld@sigma, x, otype="c")
    Lp <- ecld.op_O(ld@sigma, x, otype="p")
    
    Lc2 <- ecld.ogf(ld, x, otype="c")
    Lp2 <- ecld.ogf(ld, x, otype="p")
    
    expect_true( max(abs(Lc2/Lc-1)) + max(abs(Lp2/Lp-1)) < eps )
})


test_that("test Lc estimate for small k,sigma of BS model",{
    ld <- ecld(lambda=1, sigma=0.1, with.RN=TRUE)
    k = seq(-0.01, 0.01, 0.001)
    Lc <- ecd.ogf(ld@ecd_RN, k)
    Lc2 <- ld@sigma/2/sqrt(pi)-k/2
    err = (Lc2-Lc)/mean(Lc)
    
    expect_true( max(abs(err)) < 0.02 )
})

# ------------------------------------------------------
ld1 <- ecld(lambda=1, sigma=ecd.mp1*0.1, mu=0.05, with.ecd=TRUE)
d1 <- ld1@ecd

test_that("test IMGF of normal, mu != 0, not BS",{
    M <- ecld.mgf(ld1)
    IMc <- ecld.imgf(ld1, x, otype="c", RN=FALSE)
    IMp <- ecld.imgf(ld1, x, otype="p", RN=FALSE)

    M2 <- ecd.imgf(d1)
    IMc2 <- ecd.imgf(d1, x)
    IMp2 <- M2 - IMc2
    expect_true( abs(M2/M-1) + max(abs(IMc2/IMc-1)) + max(abs(IMp2/IMp-1)) < eps )
})

