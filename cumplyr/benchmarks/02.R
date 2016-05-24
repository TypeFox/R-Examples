benchmark02 <- function(N)
{
  Subject <- 1:N
  Block <- 1:2
  Run <- 1:3
  Trial <- 1:5
  
  positions <- cartesian_product(c('Subject', 'Block', 'Run', 'Trial'))
  M <- nrow(positions)
  positions <- transform(positions, x = sample(1:20, M, replace = TRUE))
  positions <- transform(positions, y = sample(1:20, M, replace = TRUE))
  
  centroid <- function(x, y)
  {
    c(mean(x), mean(y))
  }

  iddply(positions,
         equality.variables = c('Subject', 'Block', 'Run'),
         upper.bound.variables = c('Trial'),
         func = function (df)
         {
           z <- with(df, centroid(x, y))
           data.frame(cx = z[1], cy = z[2])
         })
}
