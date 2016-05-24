library('cumplyr')

data <- data.frame(Name = c('A', 'B', 'C', 'D', 'A'), Time = 1:5, Value = seq(1, 9, by = 2))

iddply(data,
       equality.variables = c('Name'),
       lower.bound.variables = c(),
       upper.bound.variables = c(),
       norm.ball.variables = list(),
       func = function (df) {with(df, mean(Value))})
