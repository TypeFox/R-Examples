
context("Test on High Kurtosis Region")

eps <- 0.001 # default tolerance of error for real number

d <- quantilize(ecd(2.910,0))
n <- "ecd(2.910,0)"

p_arr <- c(0.0001, 0.001, 0.01, 
           seq(0.02, 0.09, by=0.01), 
           seq(0.1, 0.9, by=0.1), 
           0.99, 0.999, 0.9999)

test_that("test half of cdf",{
    c <- ecd.cdf(d,0)
    expect_true(abs(c-0.5) < eps)
})

test_that("test half of ccdf",{
    c <- ecd.ccdf(d,0)
    expect_true(abs(c-0.5) < eps)
})

for (p in p_arr) {
    test_that(paste(n,"test quantile function at p=",p),{
        df <- pec(qec(p, d), d)/p-1
        if (abs(df) >= eps) {
            print(paste(n, "p=", p, "err", df))
        }
        expect_true(abs(df) < eps)
    })
}
    