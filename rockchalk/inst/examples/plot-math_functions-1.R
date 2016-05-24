## Paul Johnson Sept. 11, 2005  POLS 707
## Plot Some functions


x <- seq(from = 0, to = 50, length.out = 200)

y <- 3 + 4 * x

plot (x, y, main = "Linear Equation", type = "l")



y <- 3 + 4*x - 0.09 * x * x

plot (x, y, main = "Quadratic Equation", type = "l")




y <- 10 - 5 * (1/x)

plot (x,  y, main = "Reciprocal with negative coefficient", type = "l")

#Recall (1/x) = x ^(-1)


y <- 10 + 5 * (1/x)

plot (x,  y, main = "Reciprocal with positive coefficient", type = "l")



y <- exp (x)

plot (x, y, main = "Exponential of x" , type = "l")


y <- exp (-x)

plot (x, y, main = "Exponential of -x", type = "l")

# Note exp(-x) is same as 1 / exp(x)

y <- 1 / exp (x)

plot (x, y, main = "Reciprocal of Exponential of x", type = "l")






y <- exp (-(x-24)^2 )

plot (x, y, main = "Exponential of x^2", type = "l")



y <- exp ( -  (1/100) *  (x-24)^2 )
plot (x, y, main = "Exponential of x^2", type = "l")




# If y has to be constrained (say, between 0 and 100 for percents)
# there are many possibilities

z <- -10 + .4*x
y <- 100 * exp(z)/ (1 + exp(z))

# same as y <- 100 / (1 + exp(-z))
plot (x,  y, main = "S-shaped curve from the logistic", type = "l")


#Look in R's stats package for functions that start with SS.
# The SSlogis function has this formula
# y = Asym/(1+exp((xmid-input)/scale))
# you should see my example creates similar
