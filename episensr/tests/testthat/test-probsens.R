context("check probabilistic sensitivity analysis")

test_that("exposure misclassification (ND): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      spca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      print = FALSE)
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      spca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 2.1756, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.7378, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 6.3929, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 2.4539, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.8722, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 13.2802, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      seexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
                      spca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      spexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
                      corr.se = .8, corr.sp = .8,
                      print = FALSE)
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      seexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
                      spca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      spexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
                      corr.se = .8, corr.sp = .8,
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 2.9099, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.6856, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 10.0550, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 3.5240, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.8107, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 47.8042, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (ND): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("uniform", c(.8, 1)),
                      spca.parms = list("uniform", c(.8, 1)),
                      print = FALSE)
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (ND): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("uniform", c(.8, 1)),
                      spca.parms = list("uniform", c(.8, 1)),
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 2.2803, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.6629, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 20.5547, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 2.4596, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.7932, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 22.1893, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("uniform", c(.8, 1)),
                      seexp.parms = list("uniform", c(.7, .95)),
                      spca.parms = list("uniform", c(.8, 1)),
                      spexp.parms = list("uniform", c(.7, .95)),
                      corr.se = .8, corr.sp = .8,
                      print = FALSE)
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("uniform", c(.8, 1)),
                      seexp.parms = list("uniform", c(.7, .95)),
                      spca.parms = list("uniform", c(.8, 1)),
                      spexp.parms = list("uniform", c(.7, .95)),
                      corr.se = .8, corr.sp = .8,
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 4.8303, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 2.1173, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 50.2278, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 5.4494, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 2.2103, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 57.4525, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (ND---logit-logistic): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("logit-logistic", c(0, .8)),
                      spca.parms = list("logit-logistic", c(0, .8)),
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 2.5231, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.6882, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 26.4179, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 2.9212, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.8739, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 31.3662, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D---logit-logistic): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("logit-logistic", c(0, 0.8)),
                      seexp.parms = list("logit-logistic", c(0, .5)),
                      spca.parms = list("logit-logistic", c(0, 0.8)),
                      spexp.parms = list("logit-logistic", c(0, .5)),
                      corr.se = .8, corr.sp = .8,
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 3.8587, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.9935, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 28.2319, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 4.7235, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.9918, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 41.8010, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND---logit-logistic): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("logit-logistic", c(0, .8)),
                      spca.parms = list("logit-logistic", c(0, .8)),
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 3.0336, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.7952, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 18.0319, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 3.6061, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.9428, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 53.8248, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D---logit-logistic): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("logit-logistic", c(0, 0.8)),
                      seexp.parms = list("logit-logistic", c(0, .5)),
                      spca.parms = list("logit-logistic", c(0, 0.8)),
                      spexp.parms = list("logit-logistic", c(0, .5)),
                      corr.se = .8, corr.sp = .8,
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 3.9263, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7675, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 39.2988, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 5.0263, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.7453, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 155.2832, tolerance = 1e-4, scale = 1)
})
####
test_that("outcome misclassification (ND---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("logit-normal", c(2.159, 0.28)),
                      spca.parms = list("logit-normal", c(2.159, .28)),
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 5.9926, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 2.6443, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 48.1013, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 6.4646, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 2.8488, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 51.8446, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("logit-normal", c(2.159, 0.28)),
                      seexp.parms = list("logit-normal", c(0, .5)),
                      spca.parms = list("logit-normal", c(2.159, 0.28)),
                      spexp.parms = list("logit-normal", c(0, .5)),
                      corr.se = .8, corr.sp = .8,
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 3.2725, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 3.2725, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 3.2725, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 3.6503, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 3.6503, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 3.6503, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("logit-normal", c(2.159, .28)),
                      spca.parms = list("logit-normal", c(2.159, .28)),
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 2.1230, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.8816, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 3.1177, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 2.3822, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 2.0578, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 3.9231, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("logit-normal", c(2.159, 0.28)),
                      seexp.parms = list("logit-normal", c(0, .5)),
                      spca.parms = list("logit-normal", c(2.159, 0.28)),
                      spexp.parms = list("logit-normal", c(0, .5)),
                      corr.se = .8, corr.sp = .8,
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 4.9560, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.2499, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 11.3726, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 7.4460, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.2844, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 69.4315, tolerance = 1e-4, scale = 1)
})

test_that("Selection bias: observed measures are correct", {
    set.seed(123)
    model <- probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      or.parms = list("triangular", c(.35, 1.1, .43)),
                      print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.7061, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.5144, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.9693, tolerance = 1e-4, scale = 1)
})

test_that("Selection bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      or.parms = list("triangular", c(.35, 1.1, .43)),
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 1.1850, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7160, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.8153, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.1793, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.6456, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 2.0608, tolerance = 1e-4, scale = 1)
})

test_that("Selection bias: adjusted measures are correct (logit)", {
    set.seed(123)
    model <- probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      or.parms = list("logit-logistic", c(0, 0.8)),
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 1.4253, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7448, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 13.8028, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.4390, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.6698, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 13.9583, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: observed measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      prev.exp = list("triangular", c(.7, .9, .8)),
                      prev.nexp = list("trapezoidal", c(.03, .04, .05, .06)),
                      risk = list("triangular", c(.6, .7, .63)),
                      corr.p = .8,
                      print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.3479, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.2757, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.4390, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 0.2180, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 0.1519, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 0.3128, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      prev.exp = list("triangular", c(.7, .9, .8)),
                      prev.nexp = list("trapezoidal", c(.03, .04, .05, .06)),
                      risk = list("triangular", c(.6, .7, .63)),
                      corr.p = .8,
                      print = FALSE)
    expect_equal(model$adj.measures[1, 1], 0.4790, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.4531, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.5082, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4785, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3769, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.6094, tolerance = 1e-4, scale = 1)
})
