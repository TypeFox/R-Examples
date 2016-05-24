
context("Test SGED Distribution")

eps <- 0.001 # default tolerance of error for real number
eps5 <- 0.00001 # high tolerance of error 

sigma <- 0.01

lambdas <- if (ecd.devel()) seq(2, 3, by=0.5) else c(2)
betas <- if (ecd.devel()) c(0, -0.2, 0.2) else c(-0.2)

for (lambda in lambdas) {
    for (beta in betas) {

        # non-zero mu, test ecld.solve
        sged1 <- ecld(lambda=lambda, sigma=sigma*ecd.mp1, beta=beta, mu=sigma/2, is.sged=TRUE)
        ld1 <- ecld(lambda=lambda, sigma=sged1@sigma, mu=sged1@mu)  # skeleton
        d1 <- ecd(lambda=lambda, sigma=sged1@sigma, mu=sged1@mu, bare.bone=TRUE) # skeleton
        
        x <- seq(-3, 3)*sigma + sged1@mu
        
        test_that(paste("test solve of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            mu = sged1@mu
            ld0p <- ecld(lambda=lambda, sigma=sigma*(1+beta), mu=mu)
            ld0n <- ecld(lambda=lambda, sigma=sigma*(1-beta), mu=mu)
        
            y1 <- ecld.ifelse(ld1, x<mu, ecld.solve(ld0n, x), ecld.solve(ld0p, x))
            y2 <- ecld.solve(sged1, x)
            err = y1-y2
            expect_true(max(abs(err)) < eps5)
        })
        
        test_that(paste("test C of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            C1 <- ecld.const(sged1)
            C2 <- ecld.sged_const(sged1)
            expect_true( abs(C1/C2-1) < eps )
        })
        
        test_that(paste("test CDF of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            C1 <- ecld.cdf(sged1, x)
            C2 <- ecld.sged_cdf(sged1, x)
            expect_true( max(abs(C1/C2-1)) < eps )
        })
        
        test_that(paste("test y_slope of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            x <- c(-2, -1, 1, 2)*sigma + sged1@mu
            dx = 0.001*sigma
            ys1 <- ecld.y_slope(sged1, x)
            ys2 <- (ecld.solve(sged1, x+dx) - ecld.solve(sged1, x))/dx
            err = ys2/ys1-1
            expect_true(max(abs(err)) < eps)
        })

        # ----------------------------------
        n = c(1,2,3,4)
        ms1 <- ecld.moment(sged1, n)
        test_that(paste("test moment integrals of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            ms2 <- ecld.sged_moment(sged1, n)
            err2 <- ecld.ifelse(sged1, abs(ms1)<10^-10, ms1-ms2, ms2/ms1-1)
            expect_true( max(abs(err2)) < eps )
        })
        
        test_that(paste("test moment formula of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            b <- sged1@beta
            s <- sged1@sigma
            G <- gamma(lambda*(n+1)/2) /gamma(lambda/2) *s^n
            ms3 <- c(2*b, 1+3*b^2, 4*b*(1+b^2), 1+10*b^2+5*b^4) *G
            err3 <- ecld.ifelse(sged1, abs(ms1)<10^-10, ms3-ms1, ms3/ms1-1)
            #print(paste("test moment formula of lambda=", lambda, "sigma=", sigma, "beta=", beta))
            #print(err3)
            expect_true( max(abs(err3)) < eps)
        })
        
        test_that(paste("test var of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            b <- sged1@beta
            s <- sged1@sigma
            G2 <- gamma(lambda*3/2) * gamma(lambda/2)
            G3 <- gamma(lambda)^2/gamma(lambda*3/2)/gamma(lambda/2)
            var <- s^2*G2*(1+3*b^2-4*b^2*G3)
            evar <- var - ecld.var(sged1)
            expect_true( evar < eps)
        })
        
        # MGF summation tests
        if (beta==0 & lambda > 2) {
            test_that(paste("test sym mgf trunc of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
                n1 <- ecld.mgf_trunc(sged1)
                n2 <- ecld.mgf_trunc(ld1)
                expect_true( abs(n1-n2) < eps)
            })
        }
        # lambda=2 solution needs to be developed
        if (beta==0) {
            test_that(paste("test sym mgf of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
                m1 <- ecld.mgf(sged1)
                m2 <- ecld.mgf(ld1)
                expect_true( abs((m1-1)/(m2-1)-1) < eps)
            })
        }
        
        # need to use sigma=0.1 for a smaller nmax
        sged2 <- ecld(lambda=lambda, sigma=0.1*ecd.mp1, beta=beta, mu=sigma/2, is.sged=TRUE)
        if (lambda > 2) {
            test_that(paste("test mgf trunc of lambda=", lambda, "sigma=", sged2@sigma, "beta=", beta),{
                nmax <- ecd.mp2f(ecld.mgf_trunc(sged2))
                n2 <- seq(2,floor(nmax*2), by=2)
                terms <- ecld.mgf_term(sged2, n2)
                nmax2 <- n2[terms==min(terms)]
                expect_true( abs(nmax-nmax2) < 2)
            })
        }
        
        
        # ------------------------------------------------------
        sged0 <- ecld(lambda=lambda, sigma=sigma*ecd.mp1, beta=beta, is.sged=TRUE)
        # always risk netural
        mu_D <- ecld.mu_D(sged0)
        sged <- ecld(lambda=lambda, sigma=sigma*ecd.mp1, beta=beta, mu=mu_D, is.sged=TRUE)

        ki <- c(-1, 0, 1)
        k <- mu_D + ki*sigma

        test_that(paste("test CDF integrate of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            C1 <- ecld.cdf(sged, k)
            C2 <- ecld.cdf_integrate(sged0, k-mu_D)
            C3 <- ecld.sged_cdf(sged, k)
            expect_true( max(abs(C1/C2-1)) + max(abs(C1/C3-1)) < eps )
        })

        # ------------------------------------------------------
        test_that(paste("test mu_D vs mgf of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            m1 <- -log(ecld.mgf(sged0))
            m2 <- ecld.mu_D(sged0)
            expect_true( abs(m1/m2-1) < eps)
        })
        test_that(paste("test mgf vs sged_mgf of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            m1 <- -log(ecld.mgf(sged0))
            m2 <- -log(ecld.sged_mgf(sged0))
            expect_true( abs((m1-m2)/sigma) < 0.005)
        })
        test_that(paste("test sged_mgf vs imgf of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            m1 <- ecld.sged_mgf(sged0)
            m2 <- ecld.sged_imgf(sged0, 0, otype="c") + ecld.sged_imgf(sged0, -1e-8, otype="p")
            expect_true( abs((log(m1)-log(m2))/sigma) < 0.005)
        })

        mc <- ecld.sged_imgf(sged, k, otype="c") # via direct integrals of imgf and ccdf
        mp <- ecld.sged_imgf(sged, k, otype="p")
        test_that(paste("test imgf sum vs direct (c) of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            m1 <- ecld.imgf(sged0, k, otype="c")
            expect_true( max(abs(m1/mc-1)) < eps)
        })
        
        test_that(paste("test imgf sum vs direct (p) of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            m2 <- ecld.imgf(sged0, k, otype="p")
            expect_true( max(abs(m2/mp-1)) < eps)
        })

        test_that(paste("test imgf integrate vs direct (c) of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            m1 <- ecld.imgf_integrate(sged, k, otype="c")
            expect_true( max(abs(m1/mc-1)) < eps)
        })
        
        test_that(paste("test imgf integrate vs direct (p) of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            m2 <- ecld.imgf_integrate(sged, k, otype="p")
            expect_true( max(abs(m2/mp-1)) < eps)
        })
        
        test_that(paste("test imgf gamma vs direct of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            m1 <- ecld.imgf_gamma(sged, k, otype="c")
            m2 <- ecld.imgf_gamma(sged, k, otype="p")
            err1 <- max(abs(m1/mc-1))
            err2 <- max(abs(m2/mp-1))
            expect_true( max(c(err1, err2)) < eps)
        })
        
        Lc <- ecld.sged_ogf(sged, k, otype="c")
        Lp <- ecld.sged_ogf(sged, k, otype="p")
        test_that(paste("test ogf sged vs integrate of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            Lc2 <- ecld.ogf_integrate(sged0, k, otype="c")
            Lp2 <- ecld.ogf_integrate(sged0, k, otype="p")
            err1 <- max(abs(Lc2/Lc-1))
            err2 <- max(abs(Lp2/Lp-1))
            #print(paste(77, lambda, sigma, beta, ecd.mp2f(err1), ecd.mp2f(err2)))
            expect_true( max(c(err1, err2)) < eps)
        })
        
        test_that(paste("test imgf sged vs gamma of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            Lc3 <- ecld.ogf_gamma(sged0, k, otype="c")
            Lp3 <- ecld.ogf_gamma(sged0, k, otype="p")
            err1 <- max(abs(Lc3/Lc-1))
            err2 <- max(abs(Lp3/Lp-1))
            #print(paste(99, lambda, sigma, beta, ecd.mp2f(err1), ecd.mp2f(err2)))
            expect_true( max(c(err1, err2)) < eps)
        })


        # ------------------------------------------------------
        # imnt
        for (n in (0:4)) {
            #IMc_k0 <- ecld.imnt(ld, 0, n, otype="c")
            #IMp_k0 <- ecld.imnt(ld, 0, n, otype="p")
            
            test_that(paste("imnt ki=0 sum to moment, lambda=", lambda, "beta=", beta),{
                M <- 1
                e1 = 1
                expect_true(abs(e1-1) < eps)
            })

        }
        
        # ------------------------------------------------------
        sged_sml <- ecld(lambda=lambda, sigma=0.001, beta=beta, is.sged=TRUE)
        test_that(paste("test mgf_trunc and diterm for small sigma of lambda=",
                        lambda, "sigma=", sged_sml@sigma, "beta=", beta),{
            di <- ecld.mgf_diterm(sged_sml, 10000)
            nmax <- ecld.mgf_trunc(sged_sml)
            expect_true(! is.na(di) & nmax > 0)
        })

        # ------------------------------------------------------




    } # beta loop
    

}
# end of big for loop

# ------------------------------------------------------
for (beta in c(-0.2, 0, 0.2)) {
    test_that(paste("test ecld.diterm_sged_dlogF_dn of beta=", beta),{
        n <- floor(20/log(1+abs(beta)))
        f1 <- ecld.diterm_sged_dlogF_dn(n, beta)
        f2 <- log(1+abs(beta))
        expect_true(abs(f2-f1) < eps)
    })
}
