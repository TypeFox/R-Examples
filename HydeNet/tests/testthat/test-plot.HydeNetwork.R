context("plot.HydeNetwork")

data(BlackJack)

test_that("plot.HydeNetwork returns a plot under default settings",
{
  expect_that(plot(BlackJack), 
              not(throws_error()))
})

test_that("plot.HydeNetwork returns a plot with custome Node settings",
{
  expect_that(plot(BlackJack,
                   customNodes = customNode(node_id = "hit1",
                                            fillcolor = "purple", shape = "circle",
                                            fontcolor = "white", height = "2",
                                            style="filled")),
              not(throws_error()))
})

test_that("HydePlotOptions",
{
  expect_that({
    HydePlotOptions(variable=list(shape = "rect", fillcolor = "#A6DBA0"),
                  determ = list(shape = "rect", fillcolor = "#E7D4E8",
                                fontcolor = "#1B7837", linecolor = "#1B7837"),
                  decision = list(shape = "triangle", fillcolor = "#1B7837",
                                  linecolor = "white"),
                  utility = list(shape = "circle", fillcolor = "#762A83", 
                                 fontcolor = "white"))
    plot(BlackJack)},
    not(throws_error()))
})

test_that("HydePlotOptions - restoreDefaults",
{
  expect_that({
    HydePlotOptions(restorePackageDefaults = TRUE)
    plot(BlackJack)},
    not(throws_error()))
})

test_that("Remove Deterministic Nodes",
{
  expect_that({
    plot(BlackJack, removeDeterm = TRUE)},
    not(throws_error()))
})