
"superpc.decorrelate" <- function (x, competing.predictors) {
foo<- lm(t(x)~., competing.predictors)
return(foo)
}

