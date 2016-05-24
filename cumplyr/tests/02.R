library('cumplyr')

data <- data.frame(Time = 1:5, Value = seq(1, 9, by = 2))

iddply(data,
       equality.variables = c('Time'),
       lower.bound.variables = c(),
       upper.bound.variables = c(),
       norm.ball.variables = list(),
       func = function (df) {with(df, mean(Value))})

iddply(data,
       equality.variables = c(),
       lower.bound.variables = c('Time'),
       upper.bound.variables = c(),
       norm.ball.variables = list(),
       func = function (df) {with(df, mean(Value))})

iddply(data,
       equality.variables = c(),
       lower.bound.variables = c(),
       upper.bound.variables = c('Time'),
       norm.ball.variables = list(),
       func = function (df) {with(df, mean(Value))})

iddply(data,
       equality.variables = c(),
       lower.bound.variables = c(),
       upper.bound.variables = c(),
       norm.ball.variables = list('Time' = 1),
       func = function (df) {with(df, mean(Value))})

iddply(data,
       equality.variables = c(),
       lower.bound.variables = c(),
       upper.bound.variables = c(),
       norm.ball.variables = list('Time' = 2),
       func = function (df) {with(df, mean(Value))})

iddply(data,
       equality.variables = c(),
       lower.bound.variables = c(),
       upper.bound.variables = c(),
       norm.ball.variables = list('Time' = 5),
       func = function (df) {with(df, mean(Value))})
