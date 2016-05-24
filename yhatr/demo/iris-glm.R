library(yhatr)
iris$Sepal.Width_sq <- iris$Sepal.Width^2
fit <- glm(I(Species)=="virginica" ~ ., data=iris)

model.require <- function() {

}

model.transform <- function(df) {
  df$Sepal.Width_sq <- df$Sepal.Width^2
  df
}
model.predict <- function(df) {
  data.frame("prediction"=predict(fit, df, type="response"))
}
yhat.config <- c(username = "your username", apikey = "your apikey", env="http://cloud.yhathq.com")

deployment <- yhat.deploy("irisModel")

lapply(1:nrow(iris), function(i) {
  print(yhat.predict("irisModel", iris[i,]))
})
