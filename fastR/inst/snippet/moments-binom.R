x <- 0:20
sum( x * dbinom(x,20,0.25) )
sum( x^2 * dbinom(x,20,0.25) )
sum( x^2 * dbinom(x,20,0.25) ) - ( sum( x * dbinom(x,20,0.25) ) )^2
