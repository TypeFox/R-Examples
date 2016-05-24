# Test plotting functions

context("Test that theme does not break a plot")

test_that("All themes result in a valid ggplot2 object", {
  p0 <- qplot(mpg, wt, data=mtcars) + theme_dpi()
  crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
  require(reshape) # for melt
  crimesm <- melt(crimes, id = 1)
  states_map <- map_data("state")
  p1 <- ggplot(crimes, aes(map_id = state)) + geom_map(aes(fill = Murder), map = states_map) + 
       expand_limits(x = states_map$long, y = states_map$lat)+ labs(title="USA Crime")
  p1 <- p1 + coord_map()
  p1a <- p1 + theme_dpi_map()
  p2a <- p1 + theme_dpi_map2()
  p3a <- p1 + theme_dpi_mapPNG()
  expect_identical(class(p0), class(p1))
  expect_identical(class(p1), class(p1a))
  expect_identical(class(p1), class(p2a))
  expect_identical(class(p1), class(p3a))
})

context("Test that regression can be autoplotted")

test_that("Autoplot works as expected for linear models", {
  a <- runif(1000)
  b <- 7*a+rnorm(1)
  mymod <- lm(b~a)
  test_plot_file <- "lmautoplot.png"
  png(test_plot_file)
  p1 <- autoplot(mymod)
  dev.off()
  expect_true(file.exists(test_plot_file))
  unlink(test_plot_file)
  
  # Multivariate
  data(mpg)
  mymod <- lm(cty~displ + cyl + drv, data=mpg)
  test_plot_file <- "lmautoplot.png"
  png(test_plot_file)
  p1 <- autoplot(mymod)
  dev.off()
  expect_true(file.exists(test_plot_file))
  unlink(test_plot_file)
  
  test_plot_file <- "lmautoplot.png"
  png(test_plot_file)
  p2 <- autoplot(mymod, which = 1:5, mfrow = c(1, 5))
  dev.off()
  expect_true(file.exists(test_plot_file))
  unlink(test_plot_file)

})


