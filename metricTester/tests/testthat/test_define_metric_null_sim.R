library(metricTester)
context("Ensure that metrics, nulls and spatial simulations are properly defined")

metrics <- defineMetrics()

nulls <- defineNulls()

sims <- defineSimulations()

test_that("metrics are class list",
{
	expect_is(metrics, "list")
})

test_that("metrics are named",
{
	expect_true(length(names(metrics)) > 1)
})

test_that("nulls are class list",
{
	expect_is(nulls, "list")
})

test_that("null are named",
{
	expect_true(length(names(nulls)) > 1)
})

test_that("simulations are class list",
{
	expect_is(sims, "list")
})

test_that("simulations are named",
{
	expect_true(length(names(sims)) > 1)
})
