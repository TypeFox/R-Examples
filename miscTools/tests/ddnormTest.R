library( miscTools )

eps <- 1e-7

x <- (-40:40)/10

## standard normal distribution
ddnorm( x )
all.equal( ddnorm(x), ( dnorm( x + eps ) - dnorm( x - eps ) ) / ( 2 * eps ) )

## normal distribution (non-standard)
x <- (0:100)/10
ddnorm( x, mean = 5, sd = 2 )
all.equal( ddnorm( x, mean = 5, sd = 2),
   ( dnorm( x + eps, mean = 5, sd = 2 ) - dnorm( x - eps, mean = 5, sd = 2 ) )
      / ( 2 * eps ) )

