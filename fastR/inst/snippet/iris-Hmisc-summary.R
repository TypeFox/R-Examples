require(Hmisc) # load Hmisc package
summary(Sepal.Length~Species,iris)             # default function is mean
summary(Sepal.Length~Species,iris,fun=median)  # median instead
