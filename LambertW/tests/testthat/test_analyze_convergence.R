context("Testing bootstrap_LambertW_fit \n")
set.seed(20)
nobs <- 1e3

yy <- rt(n = nobs, df = 5)
mod.igmm <- IGMM(yy, type = "h")

sample.sizes <- round(c(0.5, 1) * nobs)
conv.analysis <- analyze_convergence(mod.igmm, R = 50, sample.size = sample.sizes)

test_that("analyze_convergence works", {
  expect_true(inherits(conv.analysis, "convergence_LambertW_fit"))
  expect_true(inherits(conv.analysis, "list"))
  expect_identical(names(conv.analysis), c("boots", "original.sample.size"))
  expect_identical(names(conv.analysis$boots), paste0(sample.sizes))
  expect_equal(length(conv.analysis$boots), length(sample.sizes))
})

sum.conv.analysis <- summary(conv.analysis, type = "norm")
kTypesCI <- c("basic", "norm", "bca", "perc")

test_that("summary on convergence results work", {
  
  expect_true(inherits(sum.conv.analysis, "summary.convergence_LambertW_fit"))
  expect_true(inherits(sum.conv.analysis, "list"),
              info = "is a list")
  expect_identical(names(sum.conv.analysis), c("boots", "cis", "estimate"),
                   info = "names are correct")
  
  expect_true(inherits(sum.conv.analysis$boots, "data.frame"),
              info = "'boots' is a data.frame")
  expect_true(inherits(sum.conv.analysis$cis, "data.frame"),
              info = "'cis' is a data.frame")
  
  expect_equal(sum.conv.analysis$estimate, mod.igmm$tau,
               info = "'estimate' element is the same as the original estimate")
  expect_error(summary(conv.analysis, type = "all"),
               info = "type = 'all' throws an error")
})

test_that("bootstrap confidence interval types lead to expected behavior in summary()", {
  for (tt in kTypesCI) {
    if (tt == "bca") {
      expect_error(summary(conv.analysis, tt),
                   info = "'bca' type for summary results in error since too little data")
    } else {
      expect_true(inherits(summary(conv.analysis, tt), 
                           "summary.convergence_LambertW_fit"),
                  info = paste("CI type", tt))
    }
  }
})

plot.tt <- plot(conv.analysis)

test_that("plotting works", {
  
  expect_true(inherits(plot.tt, "plot.convergence_LambertW_fit"))
  expect_true(inherits(plot.tt, "list"))
  
  for (pp in names(plot.tt)) {
    expect_true(inherits(plot.tt[[pp]], "ggplot"),
                info = paste("'", pp, "' plot"))
  }
})

