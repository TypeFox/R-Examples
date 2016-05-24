benchmark01 <- function(N)
{
  data <- data.frame(Name = 1:N,
                     Time = 1:N,
                     Value = seq(1, 2 * N - 1, by = 2))
  
  iddply(data,
         equality.variables = c('Name'),
         lower.bound.variables = c(),
         upper.bound.variables = c(),
         norm.ball.variables = list(),
         func = function (df) {with(df, mean(Value))})
}
