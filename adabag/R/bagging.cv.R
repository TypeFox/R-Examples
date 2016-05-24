bagging.cv <-
function ( formula, data,v=10, mfinal=100,control) 
{
vardep<-data[,as.character(formula[[2]])]
n <- length(vardep)
#para validacion cruzada 2<v<n
if(v>n) stop(" v should be in [2, n]")
if(v<2) stop(" v should be in [2, n]")

predclass <- rep("O",n)

    for (i in 1:v) {
        test <- v * (0:floor(n/v)) + i
        test <- test[test < n + 1]
       fit <- bagging(formula, data[-test,],mfinal, control=control)
	predclass[test] <- predict.bagging(fit, data[test,])$class
    }

   # para que devuelva la matriz de confusion
tabla <- table(predclass, vardep, dnn=c("Predicted Class", "Observed Class")) 

# Para que devuelva el error en newdata
error<- 1- sum(predclass== vardep)/n

output<- list(class=predclass, confusion=tabla, error=error)

}

