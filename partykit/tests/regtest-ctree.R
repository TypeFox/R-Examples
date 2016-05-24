## package
library("partykit")

## iris data
data("iris", package = "datasets")
irisct <- ctree(Species ~ ., data = iris)
print(irisct)
table(fit = predict(irisct), true = iris$Species)

## airquality data
data("airquality", package = "datasets")
airq <- subset(airquality, !is.na(Ozone))
airqct <- ctree(Ozone ~ ., data = airq)
print(airqct)
sum((airq$Ozone - predict(airqct))^2)
