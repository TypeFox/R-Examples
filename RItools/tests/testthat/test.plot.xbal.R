################################################################################
# Tests for plotting functions
################################################################################

library("testthat")

context("Plotting Functions")

shouldplot <- function() {
  Sys.getenv("SHOULDPLOT")[1] == "1"
}

test_that("Basic plot", {
  # create a quick xBalance result
  set.seed(20121119)
  Z <- rep(c(0,1), 10)
  df <- data.frame(Z = Z,
                   X1 = rnorm(20, mean = Z*2),
                   X2 = rnorm(20, mean = Z*3),
                   X3 = rnorm(20, mean = Z * -1),
                   X4 = rnorm(20, mean = Z * -0.5))

  xb <- xBalance(Z ~ ., data = df)

  if (shouldplot()) {

    x11() # this might be more common than quartz
    expect_true(dev.capabilities()$capture)

    plot(xb)
    p1 <- dev.capture()
    dev.off()

    x11()
    plot(xb)
    p2 <- dev.capture()
    expect_true(class(p1) == "matrix" & class(p2) == "matrix")
    expect_identical(p1, p2) # just to prove that the same plot twice is really identical
    dev.off()

    opts <- options(warn = 2)
    # has a an argument to make only a right-sided abs difference plot
    # it will be a warning if absolute isn't a parameter
    x11()
    plot(xb, absolute = T)

    p2 <- dev.capture()
    expect_true(!identical(p1,p2))
    dev.off()

    options(opts)

    x11()
    # has an argument to order the variables (from bottom to bottom)
    plot(subset(xb, vars = c("X1", "X2", "X4", "X3")), ordered = F)
    expect_true(!identical(p1, dev.capture()))
    dev.off()

    # the order the data based on the selected variable
    x11()
    plot(xb, ordered = T)
    expect_true(!identical(p1, dev.capture()))
    dev.off()
    # note: the order should change when absolute = T, but we can't really test it as the plot is then different for two reasons
    # one way to test this would be to have a helper function that creates an array for something else to plot -- and the intermediate data could be checked

    # just a sanity check to make sure that the previous dev.capture tests worked
    x11()
    plot(xb)
    expect_identical(p1, dev.capture())
    dev.off()

    ### Error checking
    expect_error(plot(xb, statistic = "foo"), "statistic")
    expect_error(plot(xb, variable.labels = c("foo")), "labels")
    expect_error(plot(xb, strata.labels = c("foo")), "labels")
  }
})

test_that("Generic balance plots", {
  # the balanceplot function takes a matrix and plots it in the expected fashion
  # this serves as a compliment to the .xbal.plot function

  testmat <- matrix(c(4,3,2,1, 3,-2,-3,2), ncol = 2,
                    dimnames = list(c("Variable 1","Variable Two","Var 3","X4"),
                        c("Stratification 1", "Stratification 2")))
  if (shouldplot()) {

    x11() # this might be more common than quartz
    expect_true(dev.capabilities()$capture)
    dev.off()

    x11()
    balanceplot(testmat)
    p1 <- dev.capture()
    dev.off()

    x11()
    balanceplot(testmat)
    p2 <- dev.capture()
    dev.off()

    expect_identical(p1,p2)

    # no segments between points
    x11()
    balanceplot(testmat, segments = FALSE)
    p.nosegs <- dev.capture()
    dev.off()

    expect_false(identical(p1, p.nosegs))

    # using shapes in the plots
    x11()
    balanceplot(testmat, shapes = 18)
    p.oneshape <- dev.capture()
    dev.off()

    expect_false(identical(p1, p.oneshape))

    x11()
    balanceplot(testmat, shapes = c(18, 18))
    p.oneshape.vec <- dev.capture()
    dev.off()

    expect_true(identical(p.oneshape, p.oneshape.vec))

    x11()
    balanceplot(testmat, shapes = matrix(18, nrow = 4, ncol = 2))
    p.oneshape.mat <- dev.capture()
    dev.off()

    expect_true(identical(p.oneshape.mat, p.oneshape))
               
    x11()
    balanceplot(testmat, shapes = c(15, 16))
    p.twoshapes <- dev.capture()
    dev.off()

    expect_true(identical(p1, p.twoshapes))

    x11()
    balanceplot(testmat, shapes = matrix(c(rep(15, 4), rep(16, 4)), nrow = 4, ncol = 2))
    p.twoshapes.mat <- dev.capture()
    dev.off()

    expect_true(identical(p1, p.twoshapes.mat))
  }

})


test_that("Issue 21: Cairo/pango errors when running plot.xbal", {
  if (capabilities()["cairo"]) {
    set.seed(20130522)

    z <- rbinom(100, size = 1, prob = 1/2)
    x <- rnorm(100)
    y <- 3 * z * x + rnorm(100)
    data <- data.frame(z, x, y)

    xb <- xBalance(z ~ x + y, data = data)
    tmpo <- tempfile()

    # at the moment, I haven't found a way to capture the stderr output from the C level pango function
    # so, we'll just have to know that if the errors appear in the output stream during testing
    # we should come here to see the test case (not the best strategy)

    tmpf <- tempfile()
    svg(tmpf)
    plot(xb)
    dev.off()
    file.remove(tmpf)
  }
})

test_that("balanceplot can group variables", {

  testmat <- matrix(c(4,3,2,1, 3,-2,-3,2), ncol = 2,
                    dimnames = list(c("Variable 1","Variable Two","Var 3","X4"),
                        c("Stratification 1", "Stratification 2")))
  grps <- c("Group 1", "Group 2", "Group 2", "Group 1")

  if (shouldplot()) {

    x11()
    balanceplot(testmat)
    p1 <- dev.capture()
    dev.off()

    x11()
    balanceplot(testmat, groups = grps)
    p2 <- dev.capture()
    dev.off()

    expect_false(identical(p1, p2))

  }
})

test_that("preparing xbalance objects for plotting, includes groups", {

  x <- data.frame(z  = rep(c(TRUE, FALSE), 50),
                  x1 = rnorm(100),
                  x2 = sample(c("A", "B", "C"), 100, replace = T),
                  x3 = sample(c("X", "Y", "Z"), 100, replace = T),
                  x4 = sample(c(T,F), 100, replace = T),
                  x5 = sample(c("A", "B", "C"), 100, replace = T))

  xb <- xBalance(z ~ x1 * x2 * x3, data = x, strata = list(foo = ~ x4, bar = ~ x5), report = 'all')

  xbp <- prepareXbalForPlot(xb)

  grps <- attr(xbp, "groups")

  expect_true(all(grps[!is.na(grps)] %in% c("x2", "x3", "x1:x2", "x2:x3", "x1:x3", "x1:x2:x3")))

  # x1 should not have a group
  expect_equal(sum(is.na(grps)), 1)


})

test_that("Plotting using RSVGTips", {

  # this is just an existance proof: it should go through without errors
  set.seed(20140137)
  
  if(require(RSVGTipsDevice)) {
    f <- tempfile()

    x <- data.frame(z  = rep(c(TRUE, FALSE), 50),
                    x1 = rnorm(100),
                    x2 = sample(c("A", "B", "C"), 100, replace = T),
                    x3 = sample(c("X", "Y", "Z"), 100, replace = T),
                    x4 = sample(c(T,F), 100, replace = T),
                    x5 = sample(c("A", "B", "C"), 100, replace = T),
                    x6 = sample(c("X", "Y", "Z", "W"), 100, replace = T))

    xb <- xBalance(z ~ x1 * x2 * x3, data = x, strata = list(foo = ~ x4, bar = ~ x5), report = 'all')
    xb$results[, "std.diff", 2] <- xb$results[, "std.diff", 2] * 2

    devSVGTips(paste0(f, "1.svg"), height = 8, width = 8)

    plot(xb)

    dev.off()


    xb2 <- xBalance(z ~ x1 * x2 * x3, data = x, strata = list(bar = ~ x5), report = 'all')

    devSVGTips(paste0(f, "2.svg"), height = 8, width = 8)

    plot(xb2)

    dev.off()

    xb3 <- xBalance(z ~ x1 * x2 * x3, data = x, strata = list(foo = ~x4, bar = ~ x5, none = NULL), report = 'all')
    xb3$results[, "std.diff", 1] <- xb$results[, "std.diff", 1] * 2
    xb3$results[, "std.diff", 2] <- xb$results[, "std.diff", 2] * 4

    devSVGTips(paste0(f, "3.svg"), height = 8, width = 8)
    plot(xb3)
    dev.off()
  }

})

