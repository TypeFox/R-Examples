
context("Test Imnt Related on Lambda Distribution")

eps <- 0.001 # default tolerance of error for real number
eps5 <- 0.00001 # high tolerance of error 

sigma <- 0.01

lambdas <- if (ecd.devel()) seq(2, 3, by=0.5) else c(2)
betas <- c(0, -0.5, 0.5)

for (lambda in lambdas) {
    for (beta in betas) {

        # always risk netural
        mu_D <- ecld.mu_D(ecld(lambda=lambda, sigma=sigma, beta=beta))
        ld <- ecld(lambda=lambda, sigma=sigma, beta=beta, mu=mu_D, with.ecd=TRUE)
        ld0 <- ecld(lambda=lambda, beta=beta, with.ecd=TRUE)

        ki <- c(-1, 0, 1)
        k <- mu_D + ki*sigma

        # ------------------------------------------------------
        # basic
        test_that(paste("test sd of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            err = ecld.sd(ld) / ecd.sd(ld@ecd) - 1
            expect_true(abs(err) < eps)
        })
        
        test_that(paste("test beta-invariant of unit dist, lambda=", lambda, "beta=", beta),{
            a <-  ecd.cdf(ld0@ecd, ld0@mu)
            b <- ecd.ccdf(ld0@ecd, ld0@mu)
            err = (b-a) * ld0@ecd@const - ld0@beta
            expect_true(abs(err) < eps)
        })

        # ------------------------------------------------------
        # imnt
        for (n in (0:4)) {
            IMc_k0 <- ecld.imnt(ld, 0, n, otype="c")
            IMp_k0 <- ecld.imnt(ld, 0, n, otype="p")
            
            test_that(paste("imnt ki=0 sum to moment, lambda=", lambda, "beta=", beta),{
                M <- if (n==0) 1 else ecld.moment(ld, n)
                e1 = IMc_k0 + IMp_k0 - M
                expect_true(abs(e1) < eps)
            })

            IMc <- ecld.imnt(ld0, ki, n, otype="c")
            IMp <- ecld.imnt(ld0, ki, n, otype="p")

            test_that(paste("imnt sigma scaling, lambda=", lambda, "beta=", beta),{
                IMc2 <- ecld.imnt(ld, ki, n, otype="c")/ld@sigma^n
                IMp2 <- ecld.imnt(ld, ki, n, otype="p")/ld@sigma^n
                e1 = IMc2/IMc-1
                e2 = IMp2/IMp-1
                expect_true(max(abs(c(e1, e2))) < eps)
            })
            
            if (n==0) {
                # CDF / CCDF
                test_that(paste("Imnt 0-th equals to ecld CDF, lambda=", lambda, "beta=", beta),{
                    e1 = IMc + IMp - 1
                    e2 = IMc - ecld.ccdf(ld0, ki)
                    e3 = IMp - ecld.cdf(ld0, ki)
                    expect_true(max(abs(c(e1, e2, e3))) < eps)
                })
                
                test_that(paste("Imnt 0-th equals to ecd CDF, lambda=", lambda, "beta=", beta),{
                    e1 = IMc + IMp - 1
                    e2 = IMc - ecd.ccdf(ld0@ecd, ki)
                    e3 = IMp - ecd.cdf(ld0@ecd, ki)
                    expect_true(max(abs(c(e1, e2, e3))) < eps)
                })
                
                test_that(paste("Imnt 0-th, ld=ld0, lambda=", lambda, "beta=", beta),{
                    IMc2 <- ecld.imnt(ld, ki, 0, otype="c")
                    IMp2 <- ecld.imnt(ld, ki, 0, otype="p")
                    e1 = IMc2/IMc - 1
                    e2 = IMp2/IMp - 1
                    e3 = IMc2 + IMp2 - 1
                    expect_true(max(abs(c(e1, e2, e3))) < eps)
                })
                
            } else {
                test_that(paste("imnt IMc+IMp equals to moment, lambda=", lambda, "beta=", beta),{
                    M <- ecld.moment(ld0, n)
                    e1 = IMc + IMp - M
                    expect_true(max(abs(e1)) < eps)
                })
            }
            
            
            if (beta==0 | lambda==2) {
                test_that(paste("Imnt", n, "-th analytic vs integrate for sym unit dist, lambda=", lambda, "beta=", beta),{
                    IMc2 <- ecld.imnt_integrate(ld0, ki, n, otype="c")
                    IMp2 <- ecld.imnt_integrate(ld0, ki, n, otype="p")
                    e1 = IMc2/IMc - 1
                    e2 = IMp2/IMp - 1
                    expect_true(max(abs(c(e1,e2))) < eps)
                })
            }
        }
        
        # ------------------------------------------------------
        # tests on IMGF
        if (beta==0 | lambda==2) {
            test_that(paste("test imgf(k) vs imnt_sum(k, Inf) of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
                Mc <- ecld.imgf(ld, k, otype="c")
                Mp <- ecld.imgf(ld, k, otype="p")
                Mc2 <- ecld.imnt_sum(ld, ki, Inf, otype="c")*exp(ld@mu)
                Mp2 <- ecld.imnt_sum(ld, ki, Inf, otype="p")*exp(ld@mu)
                
                e1 <- max(abs(Mc2/Mc-1))
                e2 <- max(abs(Mp2/Mp-1))
                expect_true( max(abs(c(e1, e2))) < eps )
            })
        }
        
        test_that(paste("test imgf(k) vs imnt(n<=2,k) apprx of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            Mc <- ecld.imgf(ld, k, otype="c")
            Mp <- ecld.imgf(ld, k, otype="p")
            Mc2 <- ecld.imnt_sum(ld, ki, 2, otype="c")*exp(ld@mu)
            Mp2 <- ecld.imnt_sum(ld, ki, 2, otype="p")*exp(ld@mu)
            
            e1 <- max(abs(Mc2/Mc-1))
            e2 <- max(abs(Mp2/Mp-1))
            e3 <- ( all(Mc >= 0.01) & all(Mp >= 0.01))
            expect_true( max(abs(c(e1, e2))) < eps & e3 )
        })
        
        # ------------------------------------------------------
        # tests on OGF
        Lc <- ecd.ogf(ld@ecd, k, otype="c")
        Lp <- ecd.ogf(ld@ecd, k, otype="p")
        
        test_that(paste("test OGF integrate of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            Lc2 <- ecld.ogf_integrate(ld, k, otype="c")
            Lp2 <- ecld.ogf_integrate(ld, k, otype="p")
            e1 = Lc2/Lc-1
            e2 = Lp2/Lp-1
            expect_true(max(abs(c(e1, e2))) < eps)
        })

        test_that(paste("test OGF imnt sum of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            # be careful to limit the order, for non-zero beta, this is a lot of integrals. 
            Lc2 <- ecld.ogf_imnt_sum(ld, k, 4, otype="c")
            Lp2 <- ecld.ogf_imnt_sum(ld, k, 4, otype="p")
            e1 = Lc2/Lc-1
            e2 = Lp2/Lp-1
            expect_true(max(abs(c(e1, e2))) < eps)
        })
        
        test_that(paste("test OGF, Lc = M(k)-e^k CCDF(k) of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            Lc2 <- ecld.imgf(ld, k, otype="c") - exp(k)*ecld.ccdf(ld, k)
            Lp2 <- -ecld.imgf(ld, k, otype="p") + exp(k)*ecld.cdf(ld, k)
            e1 = Lc2/Lc-1
            e2 = Lp2/Lp-1
            expect_true(max(abs(c(e1, e2))) < eps)
        })
        
        if (lambda==2) {
            test_that(paste("test OGF analytical of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
                Lc2 <- ecld.ogf(ld, k, otype="c")
                Lp2 <- ecld.ogf(ld, k, otype="p")
                e1 = Lc2/Lc-1
                e2 = Lp2/Lp-1
                expect_true(max(abs(c(e1, e2))) < eps)
            })
        }
        
        if (beta==0) {
            test_that(paste("test OGF gamma of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
                Lc2 <- ecld.ogf_gamma(ld, k, otype="c")
                Lp2 <- ecld.ogf_gamma(ld, k, otype="p")
                e1 = Lc2/Lc-1
                e2 = Lp2/Lp-1
                expect_true(max(abs(c(e1, e2))) < eps)
            })
            
            test_that(paste("test OGF log-slope at mu for lambda=", lambda, "sigma=", sigma, "beta=", beta),{
                Uc <- ecld.ogf_log_slope(ld, mu_D, otype="c")
                Up <- ecld.ogf_log_slope(ld, mu_D, otype="p")
                
                v <- sqrt(gamma(lambda*3/2) * gamma(lambda/2))
                f <- function(n) sigma^(n-1)/gamma(n+1)*gamma(lambda*(n+1)/2)
                qc <- sum(ecd.mpfr(sapply(1:10, f)))
                Uc2 <- -v/qc
                
                g <- function(n) (-sigma)^(n-1)/gamma(n+1)*gamma(lambda*(n+1)/2)
                qp <- sum(ecd.mpfr(sapply(1:10, g)))
                Up2 <- v/qp
                
                e1 = Uc2/Uc-1
                e2 = Up2/Up-1
                expect_true( max(abs(c(e1, e2))) < eps )
            })

        }
        
        test_that(paste("test OGF log-slope for lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            Uc <- ecld.ogf_log_slope(ld, k, otype="c")
            Up <- ecld.ogf_log_slope(ld, k, otype="p")
            
            ld0 <- ecld(lambda=ld@lambda, sigma=ld@sigma, beta=ld@beta, mu=mu_D)
            sd0 <- ecld.sd(ld0)
            
            Uc2 <- -sd0*exp(k)*ecld.ccdf(ld0,k) / ecld.ogf(ld0, k, otype="c")
            Up2 <-  sd0*exp(k)*ecld.cdf(ld0,k) / ecld.ogf(ld0, k, otype="p")
            e1 = Uc2/Uc-1
            e2 = Up2/Up-1
            expect_true( max(abs(c(e1, e2))) < eps )
        })
        
        test_that(paste("test O.V=I of lambda=", lambda, "sigma=", sigma, "beta=", beta),{
            s1 <- ecld.op_V(Lc, k, otype="c")
            s2 <- ecld.op_V(Lp, k, otype="p")
            Lc2 <- ecld.op_O(s1, k, otype="c")
            Lp2 <- ecld.op_O(s1, k, otype="p")
            e1 = Lc2/Lc-1
            e2 = Lp2/Lp-1
            e3 = s1/s2-1
            expect_true(max(abs(c(e1, e2, e3))) < eps)
        })


        # -----------------------------------------------------------------------
        # small sigma limit
        ld_sml0 = ecld(lambda=lambda, sigma=0.001, beta=beta)
        mu_D_sml = ecld.mu_D(ld_sml0)
        ld_sml = ecld(lambda=lambda, sigma=0.001, beta=beta, mu=mu_D_sml)

        test_that(paste("test OGF 1st imnt at mu_D of lambda=", lambda, "sigma=", ld_sml@sigma, "beta=", beta),{
            mu = ld_sml@mu
            Lc0 <- ecld.ogf(ld_sml, mu, otype="c")
            Lc1 <- exp(mu)*ecld.imnt(ld_sml, 0, 1, otype="c")
            Lp0 <- ecld.ogf(ld_sml, mu, otype="p")
            Lp1 <- -exp(mu)*ecld.imnt(ld_sml, 0, 1, otype="p")
            
            e1 = Lc1/Lc0-1
            e2 = Lp1/Lp0-1
            expect_true(max(abs(c(e1, e2))) < 0.004)
        })
        
        # notice k and ki is changed below this test
        test_that(paste("test OGF 0-1st imnt of lambda=", lambda, "sigma=", ld_sml@sigma, "beta=", beta),{
            mu = ld_sml@mu
            s = ld_sml@sigma
            ki <- seq(-5, 5, by=1)
            k <- mu + ki*s
            
            Lc0 <- ecld.ogf(ld_sml, k, otype="c")
            Lc1 <- exp(mu)*(ecld.imnt(ld_sml, ki, 1, otype="c") - ki*s*ecld.imnt(ld_sml, ki, 0, otype="c"))
            Lp0 <- ecld.ogf(ld_sml, k, otype="p")
            Lp1 <- -exp(mu)*(ecld.imnt(ld_sml, ki, 1, otype="p") - ki*s*ecld.imnt(ld_sml, ki, 0, otype="p"))
            
            e1 = Lc1/Lc0-1
            e2 = Lp1/Lp0-1
            expect_true(max(abs(c(e1, e2))) < 0.01)
        })
        

    }
    
    # Lc log slope for sigma->0
    test_that(paste("test OGF log-slope at mu for lambda=", lambda, "sigma->0"),{
        lds0 <- ecld(lambda=lambda, sigma=0.0001*ecd.mp1, beta=0)
        mu_D_s0 <- ecld.mu_D(lds0)
        Uc <- ecld.ogf_log_slope(lds0, mu_D_s0, otype="c")
        Up <- ecld.ogf_log_slope(lds0, mu_D_s0, otype="p")
        
        U0 <- sqrt(gamma(lambda*3/2) * gamma(lambda/2))/gamma(lambda)
        e1 = -Uc/U0-1
        e2 = Up/U0-1
        expect_true( max(abs(c(e1, e2))) < eps )
    })

}
# end of big for loop

# ------------------------------------------------------

