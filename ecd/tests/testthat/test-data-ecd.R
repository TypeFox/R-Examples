
library(xts)

context("Test on Sample Data")

conf <- .ecd.data_config()
dji <- ecd.data()
djit <- ecd.data_stats("dji")
dji_attr <- xtsAttributes(djit)

eps <- 0.001 # default tolerance of error for real number

test_that("test one day in dji history explicitly",{
  	expect_true(as.numeric(dji["19350103", "Close"]) == 105.14)
})

test_that("test log return of one day in dji explicitly",{
  	r1 <- as.numeric(dji["19350104", "logr"])
  	r2 <- log(104.69/105.14)
  	expect_true(abs(log(r1/r2)) < eps)
})

for (symbol in c("dji", "chf", "gold", "r10y", "vix", "wti")) {
    test_that(paste("test one day in ", symbol, "history"),{
        c <- .ecd.data_config(symbol)
        if(!is.na(c$test_date)) {
            ts <- ecd.data(symbol)
            v1 <- as.numeric(ts[c$test_date, "Close"])
            v2 <- c$test_val
            expect_true(v1==v2)
        }
    })
}

# histogram

djih <- dji_attr$histuple
djih2 <- ecd.manage_hist_tails(djih, c(1,2))

test_that("test hx of left tail of histgram in dji",{
  	expect_true(abs( djih$hx[1] - -0.25 ) < eps)
})

test_that("test hy of left tail of histgram in dji",{
  	expect_true(djih$hy[1] == 1)
})

test_that("test hx of right tail of histgram in dji",{
  	expect_true(abs( rev(djih$hx)[1] - 0.11) < eps)
})

test_that("test hy of right tail of histgram in dji",{
  	expect_true(rev(djih$hy)[1] == 2)
})

# merged histogram

test_that("test hx of left tail of merged histgram in dji",{
  	expect_true(abs( djih2$hx[1] - -0.09) < eps)
})

test_that("test hy of left tail of merged histgram in dji",{
  	expect_true(djih2$hy[1] == 4)
})

test_that("test hx of right tail of merged histgram in dji",{
  	expect_true(abs( rev(djih2$hx)[1] - 0.07) < eps)
})

test_that("test hy of right tail of merged histgram in dji",{
  	expect_true(rev(djih2$hy)[1] == 6)
})

# lag stats

lags <- c(1,2,4,8)
lagstats <- xtsAttributes(ecd.ts_lag_stats("dji", lags))$lagstats

test_that("test scaling of var in dji lagstats",{
    fit <- lm(log(lagstats$var) ~ log(lagstats$lags))
    expect_true(abs(coef(fit)[2]-1.0) <= 0.02)
})

test_that("test kurtosis(2 days) in dji lagstats",{
    k2 <- lagstats$kurtosis[2] # within 10% of 21
    expect_true(abs(k2/21-1.0) <= 0.1)
})

# option data
spx <- ecd.read_csv_by_symbol("spxoption2")
test_that("test required dates in spx option data",{
    dates <- unique(spx$TRADE_DT)
    expect_true("20150616" %in% dates & "20150618" %in% dates)
})

test_that("test num of strikes in spx 1-day option data",{
    d1 <- subset(spx, TRADE_DT == "20150618" & EXPR_DT == "20150619")
    num <- length(d1$STRK_PRC)
    expect_true(num == 514)
})
