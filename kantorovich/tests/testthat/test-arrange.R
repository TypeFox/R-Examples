context("arrange function")

test_that("test arrange in numeric class with same lengths", {
  mu <- setNames(c(1/2, 1/2), c("a", "b"))
  nu <- setNames(c(1/4, 3/4), c("a", "b"))
  expect_identical(arrange_names(mu, nu), list(mu=mu, nu=nu))
  #
  mu <- c(1/2, 1/2)
  nu <- c(1/4, 3/4)
  expect_identical(arrange_names(mu, nu),
                   structure(list(mu = structure(c(0.5, 0.5), .Names = c("1", "2")),
                                  nu = structure(c(0.25, 0.75), .Names = c("1", "2"))),
                             .Names = c("mu", "nu")))
  #
  mu <- setNames(mu, 1:2)
  expect_identical(arrange_names(mu, nu),
                   structure(list(mu = structure(c(0.5, 0.5), .Names = c("1", "2"
                   )), nu = structure(c(0.25, 0.75), .Names = c("1", "2"))), .Names = c("mu", "nu")))
  #
  mu <- setNames(mu, c("a", "b"))
  expect_error(arrange_names(mu, nu))
  #
  nu <- setNames(nu, c("b", "a"))
  expect_identical(arrange_names(mu, nu), list(mu=mu, nu=nu))
  #
  mu <- c(a=1/2, b=1/2)
  nu <- c(1/4, b=3/4)
  expect_error(arrange_names(mu, nu))
})


test_that("test arrange in bigq class with same lengths", {
    mu <- setNames(as.bigq(c(1/2, 1/2)), c("a", "b"))
    nu <- setNames(as.bigq(c(1/4, 3/4)), c("a", "b"))
    expect_identical(arrange_names(mu, nu), list(mu=mu, nu=nu))
    #
    mu <- as.bigq(c(1/2, 1/2))
    nu <- as.bigq(c(1/4, 3/4))
    expect_identical(arrange_names(mu, nu),
                     structure(list(mu = structure(as.raw(c(0x02, 0x00, 0x00, 0x00,
                                                            0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
                                                            0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00,
                                                            0x00, 0x00)), class = "bigq", denominator = as.raw(c(0x02, 0x00,
                                                                                                                 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02,
                                                                                                                 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
                                                                                                                 0x02, 0x00, 0x00, 0x00)), .Names = c("1", "2", NA, NA, NA, NA,
                                                                                                                                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                                                                                                                                      NA, NA, NA, NA, NA, NA)), nu = structure(as.raw(c(0x02, 0x00,
                                                                                                                                                                                                        0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01,
                                                                                                                                                                                                        0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
                                                                                                                                                                                                        0x03, 0x00, 0x00, 0x00)), class = "bigq", denominator = as.raw(c(0x02,
                                                                                                                                                                                                                                                                         0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
                                                                                                                                                                                                                                                                         0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
                                                                                                                                                                                                                                                                         0x00, 0x04, 0x00, 0x00, 0x00)), .Names = c("1", "2", NA, NA,
                                                                                                                                                                                                                                                                                                                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                                                                                                                                                                                                                                                                                                    NA, NA, NA, NA, NA, NA, NA, NA))), .Names = c("mu", "nu")))
    #
    mu <- setNames(mu, 1:2)
    expect_identical(arrange_names(mu, nu),
                     structure(list(mu = structure(as.raw(c(0x02, 0x00, 0x00, 0x00,
                                                            0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
                                                            0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00,
                                                            0x00, 0x00)), class = "bigq", denominator = as.raw(c(0x02, 0x00,
                                                                                                                 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02,
                                                                                                                 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
                                                                                                                 0x02, 0x00, 0x00, 0x00)), .Names = c("1", "2", NA, NA, NA, NA,
                                                                                                                                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                                                                                                                                      NA, NA, NA, NA, NA, NA)), nu = structure(as.raw(c(0x02, 0x00,
                                                                                                                                                                                                        0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01,
                                                                                                                                                                                                        0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
                                                                                                                                                                                                        0x03, 0x00, 0x00, 0x00)), class = "bigq", denominator = as.raw(c(0x02,
                                                                                                                                                                                                                                                                         0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
                                                                                                                                                                                                                                                                         0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
                                                                                                                                                                                                                                                                         0x00, 0x04, 0x00, 0x00, 0x00)), .Names = c("1", "2", NA, NA,
                                                                                                                                                                                                                                                                                                                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                                                                                                                                                                                                                                                                                                    NA, NA, NA, NA, NA, NA, NA, NA))), .Names = c("mu", "nu")))
    #
    mu <- setNames(mu, c("a", "b"))
    expect_error(arrange_names(mu, nu))
    #
    nu <- setNames(nu, c("b", "a"))
    expect_identical(arrange_names(mu, nu), list(mu=mu, nu=nu))
    #
    mu <- setNames(as.bigq(c(1/2, 1/2)), c("a", "b"))
    nu <- setNames(as.bigq(c(1/4, 3/4)), c("", "b"))
    expect_error(arrange_names(mu, nu))
})


test_that("test arrange in numeric class with different lengths", {
  mu <- setNames(c(1/3, 2/3), c("a", "b"))
  nu <- setNames(c(1/4, 2/4, 1/4), c("a", "b", "c"))
  expect_identical(arrange_names(mu, nu),
                   list(mu=setNames(c(1/3, 2/3, 0), c("a", "b", "c")), nu=nu))
  #
  mu <- setNames(c(1/3, 1/3, 1/3), c("a", "b", "c"))
  nu <- setNames(c(1/4, 3/4), c("a", "b"))
  expect_identical(arrange_names(mu, nu),
                   list(mu=mu, nu=setNames(c(1/4, 3/4, 0), c("a", "b", "c"))))
  #
  mu <- c(1/3, 1/3, 1/3)
  nu <- c(1/4, 3/4)
  expect_identical(arrange_names(mu, nu),
                   list(mu=setNames(mu,1:3), nu=setNames(c(1/4, 3/4, 0), 1:3)))
  #
  mu <- setNames(c(1/3, 1/3, 1/3), 1:3)
  nu <- c(1/4, 3/4)
  expect_identical(arrange_names(mu, nu),
                   list(mu=setNames(mu,1:3), nu=setNames(c(1/4, 3/4, 0), 1:3)))
  #
  mu <- setNames(mu, c("a", "b", "c"))
  expect_error(arrange_names(mu, nu))
  #
  nu <- setNames(nu, c("b", "a"))
  expect_identical(arrange_names(mu, nu),
                   list(mu=mu, nu=setNames(c(1/4, 3/4, 0), c("b", "a", "c"))))
  #
  mu <- c(a=1/3, b=1/3, c=1/3)
  nu <- c(1/4, b=3/4)
  expect_error(arrange_names(mu, nu))
})


test_that("test arrange in bigq class with different lengths", {
    mu <- setNames(as.bigq(c(1/4, 3/4)), c("a", "b"))
    nu <- setNames(as.bigq(c(1/8, 1/8, 3/4)), c("a", "b", "c"))
    expect_identical(arrange_names(mu, nu),
                     list(mu=setNames(as.bigq(c(1/4, 3/4, 0)), c("a", "b", "c")),
                          nu=setNames(as.bigq(c(1/8, 1/8, 3/4)), c("a", "b", "c"))))
    #
    mu <- setNames(as.bigq(c(1/8, 1/8, 3/4)), c("a", "b", "c"))
    nu <- setNames(as.bigq(c(1/4, 3/4)), c("a", "b"))
    expect_identical(arrange_names(mu, nu),
                     list(mu=setNames(as.bigq(c(1/8, 1/8, 3/4)), c("a", "b", "c")),
                          nu=setNames(as.bigq(c(1/4, 3/4, 0)), c("a", "b", "c"))))
    #
    mu <- as.bigq(c(1/4, 3/4))
    nu <- as.bigq(c(1/8, 1/8, 3/4))
    expect_identical(arrange_names(mu, nu),
                     list(mu=setNames(as.bigq(c(1/4, 3/4, 0)), 1:3),
                          nu=setNames(as.bigq(c(1/8, 1/8, 3/4)), 1:3)))
    #
    mu <- setNames(mu, 1:2)
    expect_identical(arrange_names(mu, nu),
                     list(mu=setNames(as.bigq(c(1/4, 3/4, 0)), 1:3),
                          nu=setNames(as.bigq(c(1/8, 1/8, 3/4)), 1:3)))
    #
    mu <- setNames(mu, c("a", "b"))
    expect_error(arrange_names(mu, nu))
    #
    nu <- setNames(nu, c("b", "a", "c"))
    expect_identical(arrange_names(mu, nu),
                     list(mu=setNames(as.bigq(c(1/4, 3/4, 0)), c("a", "b", "c")), nu=nu))
    #
    mu <- setNames(as.bigq(c(1/2, 1/2)), c("a", "b"))
    nu <- setNames(as.bigq(c(1/4, 0, 1/4)), c("", "b", "c"))
    expect_error(arrange_names(mu, nu))
    expect_error(arrange_names(nu, mu))
})
