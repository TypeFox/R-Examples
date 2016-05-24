library("papeR")
context("helper functions")

############################################################
## mySapply
############################################################
set.seed(1234)
df <- data.frame(x = rnorm(100), y = rnorm(100))

test_that("mySapply dispatches correctly on data.frames", {
              expect_equal(names(mySapply(df, mean)), c("x", "y"))
              expect_equal(mySapply(df, mean), sapply(df, mean))
          })

test_that("mySapply dispatches correctly on vectors", {
              expect_equal(length(mySapply(df$x, mean)), 1)
              expect_equal(mySapply(df$x, mean), mean(df$x))
          })


############################################################
## ADD TESTS FOR Anova.lme()
############################################################





############################################################
## get and set options
############################################################
df <- set_options(df, a = 1, b = "B", class = "test")

test_that("options and class are set correctly", {
              expect_equal(attr(df, "table.options"),
                           list(a = 1, b = "B"))
              expect_equal(class(df),
                           c("test", "data.frame"))
          })

test_that("class needs to be set", {
              expect_error(set_options(df, a = 1, b = list(x = 0, y = "y")))
          })

test_that("options are correctly obtained", {
              expect_equal(get_option(df, "a"), 1)
              expect_equal(get_option(df, "b"), "B")
              expect_equal(get_option(df, 1), 1)
              expect_equal(get_option(df, 2), "B")
              ## only one option can be obtained at once
              expect_error(get_option(df, 1:2))
          })



############################################################
## ADD TESTS FOR confint.lme
############################################################



############################################################
## ADD TESTS FOR confint.mer
############################################################


############################################################
## ADD TESTS FOR refit_model
############################################################


############################################################
## format.perc
############################################################

test_that("format.perc works", {
              expect_equal(format.perc(c(0, 0.123456, 0.1, 1), digits = 2),
                           paste0(c(0, 12, 10, 100), " %"))
              expect_equal(format.perc(c(0, 0.123456, 0.1, 1), digits = 3),
                           paste0(c("0.0", "12.3", "10.0", "100.0"), " %"))
              expect_equal(format.perc(c(0, 0.0001234, 0.1, 1), digits = 1),
                           paste0(c("0.00", "0.01", "10.00", "100.00"), " %"))

              expect_error(format.perc(c(0, 0.13456, 0.1, 1), digits = -1))
              expect_error(format.perc("A", 3))
          })


############################################################
## NAtoLvl
############################################################

test_that("NAs are correctly replaced for factors",{
              x <- xNA <- factor(rep(1:3, each = 2))
              xNA[2] <- NA

              expect_equal(NAtoLvl(x, na.lab = "missing"), x)
              expect_equal(levels(NAtoLvl(xNA, na.lab = "missing")),
                           c(1, 2, 3, "missing"))
              expect_equal(as.vector(table((NAtoLvl(xNA, na.lab = "missing")))),
                           c(1, 2, 2, 1))
          })

############################################################
## keep_levels
############################################################



############################################################
## check_which
############################################################

test_that("NAs are correctly replaced for factors",{
              expect_equal(check_which(NULL, df, "extract"), 1:2)
              expect_equal(check_which(1, df, "extract"), 1)
              expect_equal(check_which(1:2, df, "extract"), 1:2)

              expect_equal(check_which("x", df, "extract"), "x")
              expect_equal(check_which(c("x", "y"), df, "extract"), c("x", "y"))

              expect_error(check_which(-1, df, "extract"),
                           paste("One cannot extract labels for",
                                 "none-existing variables"))
              expect_error(check_which(3, df, "extract"),
                           paste("One cannot extract labels for",
                                 "none-existing variables"))
              expect_error(check_which(1.2, df, "extract"),
                           paste("One cannot extract labels for",
                                 "none-existing variables"))

              expect_error(check_which("z", df, "extract"),
                           paste("One cannot extract labels for",
                                 "none-existing variables\n  Variables",
                                 "not found in data set:\n\tz"))
              expect_error(check_which(c("a", "z"), df, "extract"),
                           paste("One cannot extract labels for",
                                 "none-existing variables\n  Variables",
                                 "not found in data set:\n\ta\n\tz"))
          })

############################################################
## get_labels
############################################################

test_that("labels are correctly extracted", {
              labels(df) <- labels(df)

              expect_equal(get_labels(df), NULL)
              expect_equal(get_labels(df$x), "x")
              expect_equal(get_labels(df$y), "y")
          })
