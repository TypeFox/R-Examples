
context("Test Incomplete Gamma and Small Sigma Limit")

eps <- 0.001 # default tolerance of error for real number
eps5 <- 0.00001 # high tolerance of error 

sigma <- 0.001

lambdas <- c(1, 2, 2.5, 3)

# ------------------------------------------------------
test_that("test Gamma(s,x) function",{
    s <- seq(0.5, 5, by=0.5)
    g1 <- ecld.gamma(s, 0)
    g2 <- gamma(s)
    expect_true(max(abs(g1-g2)) < eps5)
})

for (lambda in lambdas) {
    test_that(paste("test Gamma(s,x) hypergeo expansion, s=", lambda),{
        x <- 10
        order <- 10
        g1 <- ecld.gamma(lambda, x)
        g2 <- ecld.gamma_hgeo(lambda, x, order)
        expect_true(abs(g2/g1-1) < eps)
    })
    
    ld0 <- ecld(lambda=lambda, sigma=0.001*ecd.mp1)
    mu_D <- ecld.mu_D(ld0)
    ld <- ecld(lambda=lambda, sigma=ld0@sigma, mu=mu_D)
    
    ki <- c(2,4)*lambda
    k <- ki*ld@sigma + ld@mu
    test_that(paste("test star OGF vs full OGF, lambda=", lambda),{
        g1 <- ecld.ogf(ld, k, otype="c")
        g2 <- ecld.ogf_star(ld, ki) *ld@sigma *exp(ld@mu)
        err = max(abs(g2/g1-1)) 
        expect_true(err < 0.02)
    })
    
    test_that(paste("test star OGF btw gamma and hgeo, lambda=", lambda),{
        g1 <- ecld.ogf_star(ld, ki)
        g2 <- ecld.ogf_star_hgeo(ld, ki, order=4)
        err = max(abs(g2/g1-1)) 
        # lambda=1 has a very high error
        expect_true(err < ifelse(lambda==1, 0.15, 0.01))
    })    

    test_that(paste("test identity of star OGF btw hgeo and exp, lambda=", lambda),{
        g1 <- ecld.ogf_star_hgeo(ld, ki, order=4)
        g2 <- ecld.ogf_star_exp(ld, ki, order=3)
        err = max(abs(g2/g1-1)) 
        expect_true(err < eps)
    })    

    test_that(paste("test star OGF btw gamma and gamma_star, lambda=", lambda),{
        ki <- c(0, 0.25, 0.5, 0.75, 1)
        g1 <- ecld.ogf_star(ld, ki)
        g2 <- ecld.ogf_star_gamma_star(ld, ki)
        err = max(abs(g2/g1-1))
        expect_true(err < eps)
    })

}

