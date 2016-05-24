
context("Test C, m_n, S estimates on Lambda Distribution")

eps <- 0.001 # default tolerance of error for real number

lambda <- seq(2,3,by=0.125) 
L = lambda-2

betas <- if (ecd.devel()) c(-0.6, -0.3, 0.3, 0.6) else c(0.6)

for (beta in betas) {

    # ------------------------------------------------------
    # C
    
    fn_C <- function(lambda) {
        ld <- ecld(lambda=lambda, beta=beta)
        ld0 <- ecld(lambda=lambda)
        R <- ecld.const(ld)/ecld.const(ld0)
        (R-1)/beta^2
    }
    
    C <- simplify2array(parallel::mclapply(lambda, fn_C))
    D <- C/0.125
    fit <- lm(D~L)
    C2 <- 0.125*(0.98-L*0.335)
    
    test_that(paste("C estimate, beta=", beta),{
        e1 <- abs(C2/C-1)
        e2 <- abs( -coef(fit)[2]/0.335 - 1)
        expect_true(max(e1) < 0.04 & e2 < 0.01)
    })
    
    # ------------------------------------------------------
    # m1
    
    fn_m1 <- function(lambda) {
        ld <- ecld(lambda=lambda, beta=beta)
        R <- ecld.moment(ld,1)
        R/beta
    }
    
    M <- simplify2array(parallel::mclapply(lambda, fn_m1))
    M2 <- 1+L*0.235
    D <- M-1
    fit <- lm(D~L)
    
    test_that(paste("m1 estimate, beta=", beta),{
        e1 <- abs(M2/M-1)
        e2 <- abs( coef(fit)[2]/0.235 - 1)
        expect_true(max(e1) < 0.01 & e2 < 0.05)
    })
    
    # ------------------------------------------------------
    # m2
    fn_m2 <- function(lambda) {
        ld <- ecld(lambda=lambda, beta=beta)
        ld0 <- ecld(lambda=lambda)
        R <- ecld.moment(ld,2)/ecld.moment(ld0,2)
        (R-1)/beta^2
    }
    
    M <- simplify2array(parallel::mclapply(lambda, fn_m2))
    M2 <- 2*gamma(lambda/2)/gamma(lambda*3/2)
    
    test_that(paste("m2 estimate, beta=", beta),{
        e1 <- abs(M2-M)
        expect_true(max(e1) < 0.06)
    })

    # ------------------------------------------------------
    # var
    fn_var <- function(lambda) {
        ld <- ecld(lambda=lambda, beta=beta, with.ecd=TRUE)
        ld0 <- ecld(lambda=lambda, with.ecd=TRUE)
        R <- ecd.var(ld@ecd)/ecd.var(ld0@ecd)
        (R-1)/beta^2
    }
    
    M <- simplify2array(parallel::mclapply(lambda, fn_var))
    M2 <- gamma(lambda/2)/gamma(lambda*3/2)
    
    test_that(paste("m2 estimate, beta=", beta),{
        e1 <- abs(M2-M)
        expect_true(max(e1) < 0.03)
    })
    
    # ------------------------------------------------------
    # skewness
    
    fn_S <- function(lambda) {
        ld <- ecld(lambda=lambda, beta=beta, with.ecd=TRUE)
        ecd.skewness(ld@ecd)
    }
    
    M <- simplify2array(parallel::mclapply(lambda, fn_S))
    M2 <- (3/sqrt(2)) * beta*(1-3/7*L) * (1-4.5/12*(1-L)*beta^2)
    
    test_that(paste("skewness estimate, beta=", beta),{
        e1 <- abs(M2/M-1)
        expect_true(max(e1) < 0.03)
    })
}


