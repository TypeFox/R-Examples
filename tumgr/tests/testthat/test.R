
#------------- generate test data-------------#
#
# # empty dataset
# a <- data.frame(cbind(name = NaN, date = NaN, size = NaN))
# # dataset with non-numeric data
# b <- data.frame(cbind(name = 'd', date = 4, size= 'd'))
# # input producing not fit data only
# name <- c(rep(1,5))
# date <- c(1, 30, 60, 90, 120)
# size <- c(1, 6, .2, 7, .1)
# c <- data.frame(cbind(name, date, size))
# # input producing non-analyzed excluded cases only
# name <- c(rep(1,2),2)
# date <- c(1, 30, 1)
# size <- c(1, .9, 1)
# d <- data.frame(cbind(name, date, size))
# # input producing not fit and non-analyzed excluded cases only
# e <- data.frame(rbind(c, d))
# # input producing included cases only, check that sumstats is data frame
# data("sampleData")
# snam <- c(720001,780001,810002)
# f <- sampleData[(sampleData$name %in% snam), ]
# # input producing included, not fit, and non-analyzed excluded cases
# tall <- data.frame(rbind(c, d, f))
# # input producing included and not fit cases
# h <- data.frame(rbind(c, f))
# # input producing included and non-analyzed excluded cases
# g <- data.frame(rbind(d, f))
# # input with improper colnames
# name <- c(rep(1,5))
# day <- c(1, 30, 60, 90, 120)
# size <- c(1, 6, .2, 7, .1)
# l <- data.frame(cbind(name, day, size))
#
#
# #------------- tests-------------------------#
# test_that("expected output models is not null",{
#
#   #check that models is not null where models should exist (ie, not error)
#   expect_false(is.null(gdrate(tall, .1, FALSE)$models))
#   expect_false(is.null(gdrate(c, .1, FALSE)$models))
#   expect_false(is.null(gdrate(d, .1, FALSE)$models))
#   expect_false(is.null(gdrate(e, .1, FALSE)$models))
#   expect_false(is.null(gdrate(f, .1, FALSE)$models))
# }
# )
#
# test_that("expected output sumstats is not null",{
#
#   #check that sumstats is not null where it should exist (ie, not error)
#   expect_false(is.null(gdrate(tall, .1, FALSE)$sumstats))
#   expect_false(is.null(gdrate(c, .1, FALSE)$sumstats))
#   expect_false(is.null(gdrate(d, .1, FALSE)$sumstats))
#   expect_false(is.null(gdrate(e, .1, FALSE)$sumstats))
#   expect_false(is.null(gdrate(f, .1, FALSE)$sumstats))
# }
# )
#
#
# test_that("expected output results is not null",{
#
#   #check that results is not null where it should exist (ie, not error)
#   expect_false(is.null(gdrate(tall, .1, FALSE)$results))
#   expect_false(is.null(gdrate(c, .1, FALSE)$results))
#   expect_false(is.null(gdrate(d, .1, FALSE)$results))
#   expect_false(is.null(gdrate(e, .1, FALSE)$results))
#   expect_false(is.null(gdrate(f, .1, FALSE)$results))
# }
# )
#
# test_that("expected output allest is not null",{
#
#   #check that results is not null where it should exist (ie, not error)
#   expect_false(is.null(gdrate(tall, .1, FALSE)$allest))
#   expect_false(is.null(gdrate(c, .1, FALSE)$allest))
#   expect_false(is.null(gdrate(d, .1, FALSE)$allest))
#   expect_false(is.null(gdrate(e, .1, FALSE)$allest))
#   expect_false(is.null(gdrate(f, .1, FALSE)$allest))
# }
# )
#
#
# test_that("input colnames are correct",{
#   expect_error(gdrate(l, .1, FALSE),
#                "please rename columns as described in help page", fixed=TRUE)
# }
# )
#
# test_that("error or no output message as expected", {
#
#   expect_error(gdrate(a, .1, FALSE),
#                "input contains no non-missing data", fixed=TRUE)
#
#   expect_error(gdrate(b, .1, FALSE),
#                "all input data must be numeric", fixed=TRUE)
#
#   expect_output(gdrate(c, .1, FALSE)$allest, "no estimates when zero included cases",
#                 fixed=TRUE)
#   expect_output(gdrate(c, .1, FALSE)$sumstats, "no estimates when zero included cases",
#                 fixed=TRUE)
#
#
#
#   expect_output(gdrate(d, .1, FALSE)$allest, "no analyzable cases in input data",
#                 fixed=TRUE)
#   expect_output(gdrate(d, .1, FALSE)$sumstats, "no analyzable cases in input data",
#                 fixed=TRUE)
#
#
#
#   expect_output(gdrate(e, .1, FALSE)$allest, "no estimates when zero included cases",
#                 fixed=TRUE)
#   expect_output(gdrate(e, .1, FALSE)$sumstats, "no estimates when zero included cases",
#                 fixed=TRUE)
# }
# )
