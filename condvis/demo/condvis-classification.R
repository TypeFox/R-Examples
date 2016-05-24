library(condvis)
## demo requires 'e1071' package
library(e1071)
data(wine)
wine$Class <- as.factor(wine$Class)
## fit a model
model <- svm(Class ~ Alcohol + Malic + Ash + Magnesium +
    Phenols + Flavanoids, data = wine)
## visualise sections along 'Alcohol' and 'Phenols'
ceplot(data = wine, model = model, S = c("Alcohol", "Phenols"), sigma = 1.5)