## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- echo=FALSE,warning=FALSE,message=FALSE,results="hide"--------------
library(pmml)
library(pmmlTransformations)
library(knitr)

## ------------------------------------------------------------------------
data(iris)
kable(head(iris,3))

## ------------------------------------------------------------------------
irisBox <- WrapData(iris)

## ------------------------------------------------------------------------
kable(head(irisBox$data,3))

## ------------------------------------------------------------------------
kable(irisBox$fieldData)

## ------------------------------------------------------------------------
irisBox <- FunctionXform(irisBox,origFieldName="Sepal.Length",
                         newFieldName="Sepal.Length.Sqrt",
                         formulaText="sqrt(Sepal.Length)")

## ------------------------------------------------------------------------
kable(head(irisBox$data,3))

## ------------------------------------------------------------------------
kable(irisBox$fieldData[6,c(1:3,14)])

## ------------------------------------------------------------------------
fit <- lm(Petal.Width ~ Sepal.Length.Sqrt, data=irisBox$data)
fit_pmml <- pmml(fit, transform=irisBox)

## ------------------------------------------------------------------------
fit_pmml[[2]] #Data Dictionary node
fit_pmml[[3]][[1]] #Mining Schema node

## ------------------------------------------------------------------------
fit_pmml[[3]][[3]]

## ------------------------------------------------------------------------
irisBox <- WrapData(iris)
irisBox <- FunctionXform(irisBox,origFieldName="Species",
                         newFieldName="Species.Setosa",
                         formulaText="if (Species == 'setosa') {1} else {0}")
kable(head(irisBox$data,3))

## ------------------------------------------------------------------------
fit <- lm(Petal.Width ~ Species.Setosa, data=irisBox$data)
fit_pmml <- pmml(fit, transform=irisBox)
fit_pmml[[3]][[3]]

## ------------------------------------------------------------------------
irisBox <- WrapData(iris)
irisBox <- FunctionXform(irisBox,origFieldName="Sepal.Length,Petal.Length",
                         newFieldName="Length.Ratio",
                         formulaText="Sepal.Length / Petal.Length")

## ------------------------------------------------------------------------
kable(head(irisBox$data,3))

## ------------------------------------------------------------------------
fit <- lm(Petal.Width ~ Length.Ratio, data=irisBox$data)
fit_pmml <- pmml(fit, transform=irisBox)

## ------------------------------------------------------------------------
fit_pmml[[2]] #Data Dictionary node
fit_pmml[[3]][[1]] #Mining Schema node

## ------------------------------------------------------------------------
fit_pmml[[3]][[3]]

## ------------------------------------------------------------------------
irisBox <- WrapData(iris)
irisBox <- FunctionXform(irisBox,origFieldName="Sepal.Length,Petal.Length",
                         newFieldName="Length.Ratio",
                         formulaText="Sepal.Length / Petal.Length")

irisBox <- FunctionXform(irisBox,origFieldName="Sepal.Length,Petal.Length,Sepal.Width",
                         newFieldName="Length.R.Times.S.Width",
                         formulaText="Length.Ratio * Sepal.Width")
kable(irisBox$fieldData[6:7,c(1:3,14)])

## ------------------------------------------------------------------------
fit <- lm(Petal.Width ~ Length.R.Times.S.Width, data=irisBox$data)
fit_pmml <- pmml(fit, transform=irisBox)


## ------------------------------------------------------------------------
fit_pmml[[2]] #Data Dictionary node
fit_pmml[[3]][[1]] #Mining Schema node

## ------------------------------------------------------------------------
fit_pmml[[3]][[3]]

## ----echo=FALSE----------------------------------------------------------

funcs <- rbind(c("+","-","/","*","^","<","<=",">",">=","&&","&","|","||","==","!=","!","ceiling","prod","log"),
c("+","-","/","*","pow","lessThan","lessOrEqual","greaterThan","greaterOrEqual","and","and","or","or","equal","notEqual","not","ceil","product","ln"))
colnames(funcs) <- funcs[1,]

kable(funcs,col.names=colnames(funcs))

## ------------------------------------------------------------------------
isIn <- function(x, ...) {
  dots <- c(...)
  if (x %in% dots) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

isIn(1,2,1,4)

## ------------------------------------------------------------------------
irisBox <- WrapData(iris)
irisBox <- FunctionXform(irisBox,origFieldName="Species",
                         newFieldName="Species.Setosa.or.Versicolor",
                         formulaText="isIn(Species,'setosa','versicolor')")

## ------------------------------------------------------------------------
kable(head(irisBox$data,3))

## ------------------------------------------------------------------------
fit <- lm(Petal.Width ~ Species.Setosa.or.Versicolor, data=irisBox$data)
fit_pmml <- pmml(fit, transform=irisBox)
fit_pmml[[3]][[3]]

## ------------------------------------------------------------------------
avg <- function(...) {
  dots <- c(...)
  return(mean(dots))
}

## ------------------------------------------------------------------------
irisBox <- WrapData(iris)
irisBox <- FunctionXform(irisBox,origFieldName="Sepal.Length,Petal.Length,Sepal.Width",
                         newFieldName="Length.Average.Ratio",
                         formulaText="avg(Sepal.Length,Petal.Length)/Sepal.Width")

## ------------------------------------------------------------------------
kable(head(irisBox$data,3))

## ------------------------------------------------------------------------
fit <- lm(Petal.Width ~ Length.Average.Ratio, data=irisBox$data)
fit_pmml <- pmml(fit, transform=irisBox)
fit_pmml[[3]][[3]]

## ------------------------------------------------------------------------
functionToPMML("1 + 2")

x <- 3
functionToPMML("foo(bar(x * y))")

## ------------------------------------------------------------------------
functionToPMML("c(1,2,3)")

## ------------------------------------------------------------------------
functionToPMML("prod(1,2,na.rm=FALSE)") #produces incorrect PMML
functionToPMML("prod(1,2)") #produces correct PMML

## ------------------------------------------------------------------------
prod(c(1,2,3))
functionToPMML("prod(c(1,2,3))")

## ------------------------------------------------------------------------
functionToPMML("pmmlT(((1+2))*(x))")

## ------------------------------------------------------------------------
functionToPMML("if(a<2) {x+3} else if (a>4) {4} else {5}")

