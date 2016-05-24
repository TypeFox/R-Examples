context("theme_optimum")
library("ggplot2")

test_that("theme_optimum - correct behaviour, defaults", {
    p <- ggplot_build(qplot(1:3, 1:3) + theme_optimum())
    expect_true(attr(p$plot$theme, "complete"))
}) 
