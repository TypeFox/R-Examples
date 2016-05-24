bagging <-
function(formula, data, mfinal=100, control,...) {

formula<- as.formula(formula)
vardep <- data[,as.character(formula[[2]])]
        n <- length(data[,1])
nclases <- nlevels(vardep)
pred<- data.frame(rep(0,n)) # Dataframe para guardar las pred, al inicio esta vacio, pero luego se va agnadiendo

     arboles <- list() #Creamos una lista para guardar los arboles
replicas <- array(0, c(n,mfinal))
	#2012-05-16 nueva medida de importancia
	#sustituye a acum
	arboles[[1]] <- rpart(formula, data = data, control = rpart.control(minsplit=1, cp=-1, maxdepth=30) ) 

#Para sacar el n de variables, este luego lo sustituye en el bucle

#if( is.numeric(nrow(arboles[[1]]$splits))=="FALSE" ) stop("change rpart.control to avoid empty trees")
	nvar<-dim(varImp(arboles[[1]], surrogates = FALSE, competes = FALSE))[1]
	imp<- array(0, c(mfinal,nvar))  #Creo una matriz para guardar el "improve" de cada variable conforme evoluciona boosting/bagging



for (m in 1:mfinal) {

#boostrap<- sample(1:n,replace=TRUE)
#fit <- rpart(formula,data=data[boostrap,], control=control)

            k2 <- 1 			#Los autores agradecen su sugerencia a Ignacio Medina
            while (k2 == 1){
            
            boostrap <- sample(1:n, replace = TRUE)
            fit <- rpart(formula, data = data[boostrap, ], control = control)
            k2 <- length(fit$frame$var)
            }



arboles[[m]] <- fit #Guardamos los arboles
replicas[,m]<-boostrap

#2014-11-19 lo sustituyo por sapply
#if(m==1){pred <- predict(arboles[[m]],data,type="class")} 
#else{pred <- data.frame(pred,predict(arboles[[m]],data,type="class"))} 


		k <- varImp(arboles[[m]], surrogates = FALSE, competes = FALSE)
		imp[m,] <-k[sort(row.names(k)), ]

}



pred<-as.data.frame(sapply (arboles, predict, data=data, type="class"))




classfinal <- array(0, c(n,nlevels(vardep)))
for (i in 1:nlevels(vardep)){
 classfinal[,i] <- matrix(as.numeric(pred==levels(vardep)[i]),nrow=n)%*%rep(1,mfinal)
}

predclass <- rep("O",n)	
#2014-11-12 Se puede hacer esto usando apply para evitar el bucle? 
#Creo la funcion "select" que en caso de empate devuelva la clase mayoritaria de entre las empatadas
#2015-07-25 modifico la funcion select para poder usar predict con unlabeled data
predclass[]<-apply(classfinal,1,FUN=select, vardep.summary=summary(vardep))

#for(i in 1:n){
#predclass[i] <- as.character(levels(vardep)[(order(classfinal[i,],decreasing=TRUE)[1])])
#if(length(which(classfinal[i,]==max(classfinal[i,])))>1)
#	{predclass[i] <-names(summary(vardep)[which(classfinal[i,]==max(classfinal[i,]))])[
#order(summary(vardep)[which(classfinal[i,]==max(classfinal[i,]))],decreasing=TRUE)[1]]
#}
#else{predclass[i] <- as.character(levels(vardep)[(order(classfinal[i,],decreasing=TRUE)[1])])} 
#}
 
#normalizar la importancia de las variables, las ponderaciones son todas iguales en bagging

	pond<-rep(1,mfinal)
	imppond<-as.vector(as.vector(pond)%*%imp)
	imppond<-imppond/sum(imppond)*100
	names(imppond)<-sort(row.names(k))


#Para que devuelva las probabilidades a posteriori
classfinal/apply(classfinal,1,sum)->votosporc


ans<- list(formula=formula,trees=arboles,votes=classfinal,prob=votosporc,class=predclass, samples=replicas, importance=imppond)

#2015-07-25 pruebo a meter las clases de vardep como atributo de la salida
attr(ans, "vardep.summary") <- summary(vardep, maxsum=700)

mf <- model.frame(formula=formula, data=data) 
terms <- attr(mf, "terms") 
ans$terms <- terms 
ans$call <- match.call()




class(ans) <- "bagging"
ans

}


