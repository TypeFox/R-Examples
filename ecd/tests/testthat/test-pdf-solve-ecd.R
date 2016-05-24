
context("Test on Solver, PDF, CDF")

eps <- 0.001 # default tolerance of error for real number

x <- c(-10,-1, -0.01, 0, 0.01, 1, 10)

obj_list <- list(
    "cusp_std" = ecd(), 
    "ecd(2,0)" = ecd(2,0), 
    "ecd(0,2)" = ecd(0,2), 
    "ecd(-2,0)" = ecd(-2,0), 
    "ecd(0,-2)" = ecd(0,-2), 
    "cusp(a=2)" = ecd.cusp(alpha=2), 
    "cusp(r=-2)" = ecd.cusp(gamma=-2),
    "ecd-kmax" = ecd(2.910, 0),
    "ecd(2,2)" = ecd(2,2), 
    "ecd(-2,2)" = ecd(-2,2), 
    "ecd(-2,sm)" = ecd(-2,-0.1), 
    "ecd(-sm,-2)" = ecd(-0.1,-2), 
    "ecd(2,2,s=1,b=.5)" = ecd(2,2,1,0.5),
    "ecd(2,2,s=.5,b=.5)" = ecd(2,2,0.5,0.5),
    "ecd(0,4)" = ecd(0,4), 
    "ecd(0,8)" = ecd(0,8), 
    "ecd(0,16)" = ecd(0,16)
)

for (n in names(obj_list)) {
    d <- obj_list[[n]]
    d2 <- quantilize(d)

    if (d@beta == 0) {
        test_that(paste("test half of cdf for",n),{
            c <- ecd.cdf(d, d@stats$mean)
            expect_true(abs(c-0.5) < eps)
        })

        test_that(paste("test solve and solve_sym for",n),{
            y <- solve(d,x)
            y2 <- solve_sym(d,x)
            df <- y-y2
            df <- ifelse(is.na(df), 100, df)
            c <- sum(abs(df))
            if (abs(c) >= eps) print(paste(n, "err", c, "sym_diff", paste(df,collapse=",")))
            expect_true(abs(c) < eps)
        })
    }
    test_that(paste("test solve and solve_trig for",n),{
        y <- solve(d,x)
        y2 <- solve_trig(d,x)
        df <- y-y2
        df <- ifelse(is.na(df), 100, df)
        c <- sum(abs(df))
        if (abs(c) >= eps) print(paste(n, "err", c, "trig_diff", paste(df,collapse=",")))
        expect_true(abs(c) < eps)
    })

    test_that(paste("test quantile function for",n),{
        p <- c(0.0001, 0.001, 0.01, 0.1, 0.6, 0.99, 0.999, 0.9999)
        x <- qec(p, d2)
        cdf <- pec(x, d2)
        df <- cdf/p-1
        df <- ifelse(is.na(df), 10000, df)
        c <- sum(abs(df))
        if (abs(c) >= eps) print(paste(n, "err", c, "diff", paste(df,collapse=",")))
        expect_true(abs(c) < eps)
    })

}

test_that("test relation between trig (3) and (4)",{
    z <- seq(-10,10,by=0.1)+0i
    a <- cos(acos(z^2-1)/3)
    b <- cosh(acosh(z^2-1)/3)
    # a and b are identical, and are real
    c <- sum(abs(a-b)) + sum(abs(Im(a))) + sum(abs(Im(b)))
    expect_true(c < eps)
})

# ------------------------------------------------------
# special test cases for ecd constructor

test_that("test special issue on ecd polar coordinate",{
    # the object can't even be created by older code
    # the issue is that the moment blows up in stats
	d <- ecd.polar(R=25.326, theta=96.99/180*pi, sigma=0.00106)
	d2 <- ecd(alpha=-3.082, gamma=16.218, sigma=0.00106)
    expect_true(abs(d@stats$kurtosis/d2@stats$kurtosis-1) < eps)
})

# ------------------------------------------------------
# test ecd.rational
test_that("test ecd.rational utility",{
    a <- NULL
    a[1] <- all(ecd.rational(2.5) == c(5, 2))
    a[2] <- all(ecd.rational(2.8) == c(14, 5))
    a[3] <- all(ecd.rational(3) == c(3, 1))
    a[4] <- all(ecd.rational(1.625) == c(13, 8))
    a[5] <- all(ecd.rational(0.625) == c(5, 8))
    expect_true(all(a))
})


