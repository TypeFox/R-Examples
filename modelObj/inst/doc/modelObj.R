### R code from vignette source 'modelObj.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: modelObj.Rnw:76-78
###################################################
library(modelObj)
object1 <- buildModelObj(model=~x1, solver.method='lm')


###################################################
### code chunk number 2: modelObj.Rnw:83-95
###################################################
mylm <- function(df,x){
  obj <- list()
  obj$lm <- lm(formula=x, data=df)
  obj$var <- "does something neat"
  class(obj) = "mylm"
  return(obj)
}
predict.mylm <- function(data,obj,ty){
  obj <- predict.lm(object=obj$lm, newdata=data, type=ty)
  obj <- exp(obj)
  return(obj)
}


###################################################
### code chunk number 3: modelObj.Rnw:99-105
###################################################
object2 <- buildModelObj(model=~x1, 
                         solver.method='mylm', 
                         solver.args=list('x'="formula", 'df'="data"),
                         predict.method='predict.mylm', 
                         predict.args=list('obj'="object", 'data'="newdata", 
                                           ty='response'))


###################################################
### code chunk number 4: modelObj.Rnw:171-172
###################################################
summary(pressure)


###################################################
### code chunk number 5: modelObj.Rnw:176-191
###################################################
exampleFun <- function(modelObj, data, Y){

    fitObj <- fit(modelObj, data, Y)

    ##Test that coef() is an available method
    cfs <- try(coef(fitObj), silent=TRUE)
    if(class(cfs) == 'try-error'){
      warning("Provided regression method does not have a coef method.\n")
      cfs <- NULL
    } 

    fitted <- predict(fitObj)^2

    return(list("fittedSq"=fitted, "coef"=cfs))
}


###################################################
### code chunk number 6: modelObj.Rnw:195-202
###################################################
ylog <- log(pressure$pressure)
objlm <- buildModelObj(model = ~temperature,
                       solver.method = "lm",
                       predict.method = "predict.lm",
                       predict.args = list("type"="response"))
fitObjlm <- exampleFun(objlm, pressure, ylog)
print(fitObjlm$coef)


###################################################
### code chunk number 7: modelObj.Rnw:206-213
###################################################
objnls <- buildModelObj(model = ~exp(a + b*temperature),
                        solver.method = "nls",
                        solver.args = list('start'=list(a=1, b=0.1)),
                        predict.method = "predict",
                        predict.args = list("type"="response"))
fitObjnls <- exampleFun(objnls, pressure, pressure$pressure)
print(fitObjnls$coef)


###################################################
### code chunk number 8: modelObj.Rnw:217-225
###################################################
objectnew <- buildModelObj(model = ~temperature, 
                           solver.method = 'mylm', 
                           solver.args=list('x'="formula", 'df'="data"),
                           predict.method='predict.mylm', 
                           predict.args=list('obj'="object", 'data'="newdata", 
                                             ty='response'))
fitObjnew <- exampleFun(objectnew, pressure, ylog)
print(fitObjnew$coef)


