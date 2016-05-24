##
## --- Test treePlot.R ---
##

context("check that treePlot returns warnings when providing incorrectly formatted phylo4d objects.")

test_that("phylo4d gives warning when there is no data", {
    phyd <- phylo4d(ape::rcoal(5), tip.data=data.frame())
    expect_warning(plot(phyd), "tree has no tip data to plot")
})

test_that("phylo4d gives warning when there is data but they can't be plotted", {
    phyd <- phylo4d(ape::rcoal(5), tip.data=data.frame(letters[1:5], letters[6:10]))
    expect_warning(plot(phyd), "only numeric data can be plotted at this time")
})

## test.treePlot <- function() {
## }

## test.plotOneTree <- function() {
## }

## test.phyloXXYY <- function() {
##     # function(phy, tip.order = NULL) 
## }

## test..bubLegendGrob <- function() {
## }

## test.drawDetails.bubLegend <- function() {
## }

## test.phylobubbles <- function() {
## }

## test.tip.data.plot <- function() {
## }

## test.plot.phylo4 <- function() {
##     # signature(x='phylo4', y='missing')
## }

