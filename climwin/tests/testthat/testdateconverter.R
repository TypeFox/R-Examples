# Test outcomes of convertdate #

# Test that no NAs are produced #
# Test that cintno starts at 1 #
# Test that bintno and cintno have the same range #
# Test that error is returned when cinterval is wrong #
# Test that error is returned when there are duplicate climate days #
# Test that the length of cintno and xvar are equal #
# Test for all possible combos of cinterval and cross TRUE or FALSE #
test_that("convertdate works (day, variable)", {

data(Mass, envir = environment())
data(MassClimate, envir = environment())

MassClimatedup         <- MassClimate
MassClimatedup[17533,] <- MassClimatedup[17532,]

test <- convertdate(bdate = Mass$Date, cdate = MassClimate$Date, xvar = MassClimate$Temp,
                    cinterval = "day", type = "variable")

expect_equal(length(which(is.na(test))), 0)
expect_equal(test$cintno[1], 1)
expect_equal(length(test$cintno), length(test$xvar))
expect_true((max(test$bintno) %in% test$cintno))
expect_error(convertdate(bdate = Mass$Date, cdate = MassClimate$Date,
                         cinterval = "R", type = "variable"))
expect_error(convertdate(bdate = Mass$Date, cdate = MassClimatedup$Date,
                         cinterval = "day", type = "variable"))

})

test_that("convertdate works (week, variable)", {
  
data(Mass, envir = environment())
data(MassClimate, envir = environment())

test <- convertdate(bdate = Mass$Date, cdate = MassClimate$Date, xvar = MassClimate$Temp,
                    cinterval = "week", type = "variable")

expect_equal(length(which(is.na(test))), 0)
expect_equal(test$cintno[1], 1)
expect_equal(length(test$cintno), length(test$xvar))
expect_true((max(test$bintno) %in% test$cintno))

})

test_that("convertdate works (month, variable)", {
  
data(Mass, envir = environment())
data(MassClimate, envir = environment())

test <- convertdate(bdate = Mass$Date, cdate = MassClimate$Date, xvar = MassClimate$Temp,
                    cinterval = "month", type = "variable")

expect_equal(length(which(is.na(test))), 0)
expect_equal(test$cintno[1], 1)
expect_equal(length(test$cintno), length(test$xvar))
expect_true((max(test$bintno) %in% test$cintno))

})

test_that("convertdate works (day, xvar2)", {
  
data(Mass, envir = environment())
data(MassClimate, envir = environment())

test <- convertdate(bdate = Mass$Date, cdate = MassClimate$Date, xvar = MassClimate$Temp,
                    xvar2 = MassClimate$Rain, cinterval = "day", type = "variable", cross = TRUE)

expect_equal(length(which(is.na(test))), 0)
expect_equal(test$cintno[1], 1)
expect_equal(length(test$cintno), length(test$xvar))
expect_equal(length(test$cintno), length(test$xvar2))
expect_true((max(test$bintno) %in% test$cintno))

})

test_that("convertdate works (week, xvar2)", {
  
data(Mass, envir = environment())
data(MassClimate, envir = environment())

test <- convertdate(bdate = Mass$Date, cdate = MassClimate$Date, xvar = MassClimate$Temp,
                    xvar2 = MassClimate$Rain, cinterval = "week", type = "variable", cross = TRUE)

expect_equal(length(which(is.na(test))), 0)
expect_equal(test$cintno[1], 1)
expect_equal(length(test$cintno), length(test$xvar))
expect_equal(length(test$cintno), length(test$xvar2))
expect_true((max(test$bintno) %in% test$cintno))

})

test_that("convertdate works (month, variable)", {
  
data(Mass, envir = environment())
data(MassClimate, envir = environment())

test <- convertdate(bdate = Mass$Date, cdate = MassClimate$Date, xvar = MassClimate$Temp,
                    xvar2 = MassClimate$Rain, cinterval = "month", type = "variable", cross = TRUE)

expect_equal(length(which(is.na(test))), 0)
expect_equal(test$cintno[1], 1)
expect_equal(length(test$cintno), length(test$xvar))
expect_equal(length(test$cintno), length(test$xvar2))
expect_true((max(test$bintno) %in% test$cintno))

})
