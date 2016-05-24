context("createPopPyramid")

test_that("createPopPyramid works", {
  data = data.frame(x=rep(1:10,2), y=rnorm(20, mean=100, sd=10), sex=rep(c('M','F'), each=10))
  
  expect_warning(p <- ggplot_build(createPopPyramid(data=data, bin='x', count='y', 
                                    divideBy='sex', values=c('M','F'))),
                 "Stacking not well defined when ymin != 0")
  
  expect_equal(nrow(p$data[[1]]), dim(data[data$sex=='M',])[[1]])
  expect_equal(nrow(p$data[[2]]), dim(data[data$sex=='F',])[[1]])
  expect_true(all(p$data[[1]][,'ymax'] == 0))
  expect_true(all(p$data[[2]][,'ymin'] == 0))
})