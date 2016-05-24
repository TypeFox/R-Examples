require(Hmisc)
summary(Sepal.Length~Species,data=iris)
summary(Sepal.Length~Species,data=iris, fun=quantile)
