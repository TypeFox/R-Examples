
predict.bagging<- function (object, newdata, newmfinal=length(object$trees), ...) 

{
    if (newmfinal > length(object$trees) | newmfinal < 1) 
        stop("newmfinal must be 1<newmfinal<mfinal")

    formula <- object$formula

 #   mfinal <- length(object$trees)
    n <- length(newdata[, 1])

#2015-07-26 lo cambio para pedict con unlabeled data 
vardep.summary<-attributes(object)$vardep.summary

   nclases <- length(vardep.summary)

    pred <-as.data.frame(sapply (object$trees[1:newmfinal], predict, newdata, type="class"))
	classfinal <- array(0, c(n, nclases))
    for (i in 1:nclases) {
        classfinal[, i] <- matrix(as.numeric(pred == names(vardep.summary)[i]), nrow = n) %*% rep(1, newmfinal)
    }
    predclass <- rep("O", n)
#2014-11-12 Se puede hacer esto usando apply para evitar el bucle? 
#Creo la funcion "select" que en caso de empate devuelve la clase mayoritaria de entre las empatadas
#2015-07-26 lo cambio para pedict con unlabeled data
predclass[]<-apply(classfinal,1,FUN=select, vardep.summary=vardep.summary)

if(sum(names(newdata)==as.character(object$formula[[2]]))==0)
{
    tabla <- NULL
    error <- NULL

}
 else{
    vardep <- newdata[, as.character(object$formula[[2]])]
    tabla <- table(predclass, vardep, dnn = c("Predicted Class", 
        "Observed Class"))
    error <- 1 - sum(predclass == vardep)/n
}

#Para que devuelva las probabilidades a posteriori
classfinal/apply(classfinal,1,sum)->votosporc

output<- list(formula=formula, votes=classfinal, prob=votosporc, class=predclass, confusion=tabla, error=error)
}

