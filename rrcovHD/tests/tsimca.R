## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcovHD))

data(iris)

## New data for prediction consisting only of the first two classes
newx <- iris[iris$Species %in% c("setosa", "versicolor"), -5]

## Qda and other classification methods will keep the levels of the grouping
##  variable, even if the new data has not objects assigned to each class.
qq <- QdaClassic(Species~., data=iris)
qq
pr <- predict(qq, newdata=newx)
pr

cs <- CSimca(Species~., data=iris, k=4)
cs
pr1 <- predict(cs)
pr1
pr2 <- predict(cs, newdata=newx)
pr2

## Prediction when in the new data there are missing values
data(fish)
newfish <- na.omit(fish[fish$Species %in% c(1, 3, 7), -7])
cs <- CSimca(Species~., data=fish, k=6, kmax=6)
cs
pr1 <- predict(cs)
pr1
pr2 <- predict(cs, newdata=newfish)
pr2
