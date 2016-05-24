# 
# Tests the way rollply handles bad arguments.
# 

context('Handling of bad arguments')

baddata <- cbind(runif(10,0,10),
                 runif(10,0,10),
                 runif(10,0,20))
colnames(baddata) <- c('lon','lat','temp')

badgrid <- cbind(seq(0,10,l=10),
                 seq(0,10,l=10))
colnames(badgrid) <- c('x','y')

test_that("Rollply handles missing variables", {
  # Error: bad class for .data
  expect_error(rollply(density(1:10), ~ x+y, summarise, temp.mean=mean(temp)))
  # Error: vars not found in .data
  expect_error(rollply(baddata, ~ x+y, summarise, temp.mean=mean(temp)))
  # Error: bad class for grid 
  expect_error(rollply(baddata, ~ x+y, summarise, temp.mean=mean(temp), 
               grid=density(1:10)))
  # Error: non-existent grid type
  expect_error(rollply(baddata, ~ x+y, summarise, temp.mean=mean(temp), 
               grid.type='nonexistent_grid'))
  # Error: vars not found in grid
  expect_error(rollply(baddata, ~ lon+lat, summarise, temp.mean=mean(temp), 
               grid=badgrid))
})
