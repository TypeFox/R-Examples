context("Estimate Quantiles for given Return Periods")

# importing the discharges at gauge "Wildungsmauer (Danube)"
# truncating the data to complete hydrological years
infile <- readlfdata("QMittelTag207373.txt", type="HZB", hyearstart = 4,
                     baseflow = FALSE)
wild <- subset(infile, hyear %in% 1996:2011)

am <- tyears(wild, dist = "wei", plot = F)$values

# Zeitreihe mit Lücken
wild.zeros <- wild
wild.zeros$flow[c(10, 500)] <- 0

load("reference-gl.RData")

test_that("94% quantile for Wildungsmauer is correct", {
  expect_equal(unname(Qxx(wild, 94)), 988.74, tolerance = 1e-2)
})

pars <- tyears(wild, dist="wei", plot = F)$parameters$wei

test_that("return period for 94% event is the same as in script GL", {
  rp <- 1/cdfwei(Qxx(wild, 94), pars)
  # => 1.424714 (i.e.: Return period of 1.4 years)
  # numerisch nicht ganz ident mit tolerance = 1e-5
  expect_equal(object = as.vector(rp),
               expected = 1.42471, tolerance = 1e-4)
})

test_that("flows for given return periods are the same as in script GL", {

  expect_equal(as.vector(quawei(f=1/20, pars)),
               668.4, tolerance = 1e-2)

  expect_equal(as.vector(quawei(f=1/5, pars)),
               772.3, tolerance = 1e-2)
})

# test_that("parameters estimated by evquantile() are correct", {
#   test_fit_weibull <- function(){
#     shape = runif(n = 1, min = 0.5, max = 2)
#     scale = runif(n = 1, min = 10, max = 500)
#     am <- rweibull(500, shape = shape, scale = scale)
#
#     y <- evfit(am, distribution = "wei", zeta = 0)
#
#    return(c(shape.pop = shape, scale.pop = scale,
#              shape.est = y$parameters$wei["delta"],
#              scale.est = y$parameters$wei["beta"],
#              shape.e = unname((y$parameters$wei["delta"] - shape) / shape),
#              scale.e = unname((y$parameters$wei["beta"] - scale) / scale)))
#   }
#
#   # tolerances computed with n=100 replications
#   error <- replicate(n = 10, expr = test_fit_weibull(), simplify = T)
#
#   expect_less_than(object = max(abs(error["shape.e", ])),
#                    expected = 0.1)
#
#   expect_less_than(object = max(abs(error["scale.e", ])),
#                    expected = 0.1)
# })




test_that("warnings are given", {

  # warn: "... using gevR instead."
  expect_warning(tyears(wild, dist = "gev", plot = F))

  # warn that zero flow observations are present
  expect_warning(tyears(wild.zeros, dist = "wei", plot = F))
})




test_that("gevR and wei behave identically", {

  # only for time series with strictly positive values
  # pelwei complains about "invalid L-moments"

  y <- suppressWarnings(tyears(wild, plot = F,
                               dist = c("wei", "gevR"),
                               event = c(1, 2, 7.9, 8, 8.1, 10)))

  rp <- y$T_Years_Event

  # quantiles for weibull und reversed GEV are identical
  expect_equal(rp[, "wei"], rp[, "gevR"], tolerance = 1e-10)

  # quantiles
  expect_equal(qua_ev(distribution = "gevR",
                      f = seq(0, 1, 0.01), para = y$parameters$gevR),
               qua_ev(distribution = "wei",
                      f = seq(0, 1, 0.01), para = y$parameters$wei),
               tolerance = 1e-10)

  # cdf
  expect_equal(cdf_ev(distribution = "gevR",
                      x = seq(200, 1300, 1), para = y$parameters$gevR),
               cdf_ev(distribution = "wei",
                      x = seq(200, 1300, 1), para = y$parameters$wei),
               tolerance = 1e-10)
})

test_that("quantiles are INF for return period = 1 year", {

  # doesn't work for "gev", because a Frechet type GEV is fitted, which
  # has a finite upper bound
  y <- suppressWarnings(tyears(wild.zeros, plot = F, zeta = 0,
                               dist = c("wei", "gum", "gevR"),
                               event = c(1)))
  expect_true(all(y$T_Years_Event == Inf))

  y <- suppressWarnings(tyears(wild, plot = F, zeta = 0,
                               dist = c("wei", "gum", "gevR"),
                               event = c(1)))
  expect_true(all(y$T_Years_Event == Inf))

})


test_that("quantiles for mixed distributions are plausible", {

  y <- suppressWarnings(tyears(wild.zeros, plot = F, zeta = 0,
                               dist = c("wei", "gum", "gevR"),
                               event = c(1, 2, 7.9, 8, 8.1, 10, 50)))

  rp <- y$T_Years_Event

  # return period of 1 has infinite quantiles
  expect_true(all(rp[rownames(rp) == 1, ] == Inf))


  # return periods lower than frequency of zero obersvations yield >0m³/s
  mask <- as.numeric(rownames(rp)) < 1 / y$freq.zeros
  expect_true(all(rp[mask, ] > 0))

  # return periods higher than frequency of zero obersvations yield 0m³/s
  mask <- as.numeric(rownames(rp)) >= 1 / y$freq.zeros
  expect_true(all(rp[mask, ] == 0))

})


# test_that("evquantiles() gives the same results as WMO manual", {
#
#   wmo <- list(
#     nicholson = list(data = c(8, 13, 2, 2.7, 9.1, 0, 0.7, 21.6, 35.1, 12.7, 1.3,
#                               9.1, 26.0,  22.7, 21.3, 7.9, 23, 6, 3.1, 4.1, 0,
#                               4.9, 2.3, 12.4, 4.4, 2.1, 3.1, 7.6, 9, 9, 20.4,
#                               13.6, 3.3, 10.3),
#                      lmom = matrix(c(9.765, 10.375, 4.721, 4.670, 0.273, 0.276),
#                                    ncol = 3, nrow = 2,
#                                    dimnames = list(c("values", "censored"),
#                                                    c("l_1", "l_2", "t_3"))),
#                      pmom = c(9.765, 8.68, 1.149),
#                      quantile = 1.2),
#
#     timbarra = list(data = c(27.3, 38, 72.9, 53.6, 45.7, 37, 21.7, 36.6, 43.1,
#                              12.6, 20.7, 66.7, 76.6, 53.1, 15.7, 66.6, 78, 53.1,
#                              62.6, 40.7, 51.6, 23.1, 23.4, 23, 6),
#                     lmom = matrix(c(41.977, 12.305, 0.045),
#                                   ncol = 3, nrow = 1,
#                                   dimnames = list(c("values"),
#                                                   c("l_1", "l_2", "t_3"))),
#                     pmom = c(41.977, 21.05, 0.146),
#                     quantile = 15.4))
#
#
#
#
#   # reproduce product moments (to test if input data is correct)
#   library(e1071)
#   pmom <- sapply(c(mean, sd, skewness), function(x) x(wmo$nicholson$data))
#   expect_true(abs(pmom[1] - wmo$nicholson$pmom[1]) < 1e-3)
#   expect_true(abs(pmom[2] - wmo$nicholson$pmom[2]) < 1e-2)
#   expect_true(abs(pmom[3] - wmo$nicholson$pmom[3]) < 1e-3)
#
#   pmom <- sapply(c(mean, sd, skewness), function(x) x(wmo$timbarra$data))
#   expect_true(abs(pmom[1] - wmo$timbarra$pmom[1]) < 1e-3)
#   expect_true(abs(pmom[2] - wmo$timbarra$pmom[2]) < 1e-2)
#   expect_true(abs(pmom[3] - wmo$timbarra$pmom[3]) < 1e-3)
#
#   # reproduce L-moments
#   nicholson <- evfit(x = wmo$nicholson$data,
#                      distribution = "wei", zeta = 0)
#
#   timbarra <- evfit(x = wmo$timbarra$data,
#                     distribution = "wei", zeta = 0)
#
#   expect_true(abs(nicholson$lmom[, c("l_1", "l_2", "t_3")] -
#                     wmo$nicholson$lmom[, c("l_1", "l_2", "t_3")]) <  1e-3)
#
#   expect_true(abs(timbarra$lmom["raw data", c("l_1", "l_2", "t_3")] -
#                     wmo$timbarra$lmom[, c("l_1", "l_2", "t_3")]) <  1e-3)
#
#   # reproduce parameters of distribution
#
#
#   # reproduce equantile of 10 year event
#   quant <- evquantile(nicholson, return.period = 10)
#   expect_equal(object = as.vector(quant$T_Years_Event),
#                expected = wmo$nicholson$quantile,
#                tolerance = 1e-2)
#
#   quant <- evquantile(timbarra, return.period = 10)
#   expect_equal(object = as.vector(quant$T_Years_Event),
#                expected = wmo$timbarra$quantile,
#                tolerance = 1e-2)
# })


test_that("same results as skripts from GL", {
  # all gl fits with pooling are incorrect, because of wrong pooling in old package

  # tyears_new_GL
  lfstat <- suppressWarnings(tyears(wild, dist = c("gevR", "wei"), plot = F))
  expect_equal2(object = lfstat$parameter[["wei"]],
                expected = reference.gl$tyears_new_GL$parameters[["wei"]])
  expect_equal2(object = lfstat$T_Years_Event[1, "wei"],
                expected = unname(reference.gl$tyears_new_GL$T_Years_Event["wei"]))


  expect_equal2(object = lfstat$parameter[["gevR"]],
                expected = reference.gl$tyears_new_GL$parameters[["gev"]])
  expect_equal2(object = lfstat$T_Years_Event[1, "gevR"],
                expected = unname(reference.gl$tyears_new_GL$T_Years_Event["gev"]))

  # tyears_MAX_D_hyear_GL-sum
  lfstat <- suppressWarnings(tyearsS(wild,  dist = c("gev", "wei"), variable = "d",
                                     aggr = sum, plot = F,
                                     threshold = function(x) quantile(x, probs = 0.06, na.rm = TRUE)))


  expect_equal2(object = lfstat$parameter[["gev"]],
                expected = reference.gl$"tyears_MAX_D_hyear_GL-sum"$parameters[["gev"]])
  expect_equal2(object = lfstat$T_Years_Event[1, "gev"],
                expected = unname(reference.gl$"tyears_MAX_D_hyear_GL-sum"$T_Years_Event["gev"]))

  # tyears_MAX_D_hyear_GL-max
#   lfstat <- suppressWarnings(tyearsS(wild,  dist = c("gev", "wei"),
#                                      pooling = pool_ic,
#                                      variable = "d", aggr = max,
#                                      plot = F,
#                                      threshold = function(x) quantile(x, probs = 0.06, na.rm = TRUE)))
#
#   expect_equal2(object = lfstat$parameter[["gev"]],
#                 expected = reference.gl$"tyears_MAX_D_hyear_GL-max"$parameters[["gev"]])
#   expect_equal2(object = lfstat$T_Years_Event[1, "gev"],
#                 expected = unname(reference.gl$"tyears_MAX_D_hyear_GL-max"$T_Years_Event["gev"]))

#   # tyears_MAX_V_hyear_GL-sum
  lfstat <- suppressWarnings(tyearsS(wild,  dist = c("gev", "wei"), variable = "v",
                                     threshold = function(x) quantile(x, probs = 0.06, na.rm = TRUE),
                                     aggr = sum,  plot = F))

  expect_equal2(object = lfstat$parameter[["gev"]],
                expected = reference.gl$"tyears_MAX_V_hyear_GL-sum"$parameters[["gev"]],
                tolerance = 1e-4)
  expect_equal2(object = lfstat$T_Years_Event[1, "gev"],
                expected = unname(reference.gl$"tyears_MAX_V_hyear_GL-sum"$T_Years_Event["gev"]),
                tolerance = 1e-4)

  # tyears_MAX_V_hyear_GL-max
#   lfstat <- suppressWarnings(tyearsS(wild, dist = c("gev", "wei"), variable = "v",
#                                      threshold = function(x) quantile(x, probs = 0.06, na.rm = TRUE),
#                                      aggr = max, pooling = pool_ic, plot = F))
#
#   expect_equal2(object = lfstat$parameter[["gev"]],
#                 expected = reference.gl$"tyears_MAX_V_hyear_GL-max"$parameters[["gev"]],
#                 tolerance = 1e-7)
#   expect_equal2(object = lfstat$T_Years_Event[1, "gev"],
#                 expected = unname(reference.gl$"tyears_MAX_V_hyear_GL-max"$T_Years_Event["gev"]),
#                 tolerance = 1e-6)
})
