b0 <- 3; b1 <- 5; sigma <- 2       # set model parameters
x <- rep(1:5,each=4)               # 4 observations at each of 5 values
e <- rnorm(length(x),sd=sigma)     # error term in the model
y <- b0 + b1*x + e                 # build response according to model
###hop:3-9
model <- lm(y~x); summary(model)
confint(model)
