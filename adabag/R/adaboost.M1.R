
boosting <-
function(formula, data,boos=TRUE, mfinal=100, coeflearn="Breiman", control, ...) {


#Exigimos que coeflearn sea uno de esos dos valores
if (!(as.character(coeflearn) %in% c("Freund","Breiman","Zhu"))){
stop("coeflearn must be 'Freund', 'Breiman' or 'Zhu' ")
}
formula<- as.formula(formula)
vardep <- data[,as.character(formula[[2]])]
        n <- length(data[,1])
nclases <- nlevels(vardep)

pesos <- rep(1/n,n)
guardarpesos <- array(0, c(n,mfinal)) #para ver los pesos de las observaciones

w <- rep(1/n,n) # desaparece el not visible binding for "<<-" que se usa en boos=F
data<-data.frame(pesos, data) #Los pesos en rpart deben ser una columna del dataframe

     arboles <- list() #Creamos una lista para guardar los arboles
pond <- rep(0,mfinal) # Un vector donde guardaremos la ponderacion de cada arbol.
pred<- data.frame(rep(0,n))

	#2012-05-16 nueva medida de importancia
	#sustituye a acum
	arboles[[1]] <- rpart(formula, data = data[,-1], control = rpart.control(minsplit=1, cp=-1, maxdepth=30) ) 
	#Para sacar el n de variables, este luego lo sustituye en el bucle
	nvar<-dim(varImp(arboles[[1]], surrogates = FALSE, competes = FALSE))[1]
	imp<- array(0, c(mfinal,nvar))  #Creo una matriz para guardar el "improve" de cada variable conforme evoluciona boosting




for (m in 1:mfinal) {
#Creamos muestras boostrap utilizando los pesos

if (boos==TRUE) {
            k <- 1		#Gracias a Ignacio Medina 2014-11-06; Evitamos arboles sin ningun corte
            while (k == 1){
            
            boostrap <- sample(1:n, replace = TRUE, prob = pesos)
            fit <- rpart(formula, data = data[boostrap, -1], control = control)
            k <- length(fit$frame$var)
            }	#Hasta aqui Gracias a Ignacio Medina 2014-11-06
		#La solucion I. Medina con boos=FALSE puede no converger


        flearn <- predict(fit,newdata=data[,-1],type="class")
ind<-as.numeric(vardep != flearn) #crear un vector indicador
err<- sum(pesos*ind)         #Calcula el error ponderado en esa iteracion

}

#limitamos el tamagno del arbol para que sean distintos?
if (boos==FALSE) {
	w<<- pesos
  	    fit <- rpart(formula=formula, data=data[,-1], weights=w, control=control) 

        flearn <- predict(fit,data=data[,-1], type="class")
ind<-as.numeric(vardep != flearn) #Crear un vector indicador
err<- sum(pesos*ind)         #Calcula el error ponderado en esa iteracion


}
# Diferenciamos entre Freund, Breiman y Zhu. 
c<- log((1-err)/err)

	if (coeflearn=="Breiman"){
	c<- (1/2)*c
	}

	if (coeflearn=="Zhu"){
	c<- c+log(nclases-1)
	}

		guardarpesos[,m]<-pesos
pesos <- pesos*exp(c*ind)
pesos<- pesos/sum(pesos)

#Si el error no es menor que la regla por defecto los pesos se inicializan
# Segun Opitz y Maclin si no 0<err<0.5 se ajustan los valores de a 3 y 0.001


maxerror<-0.5
eac<-0.001 # minimum fraction of error above 0 or under maxerror 
#maxerror<-min(1-max(summary(vardep))/sum(summary(vardep)), 0.5)
#maxerror<-1-max(summary(vardep))/sum(summary(vardep))

		if (coeflearn=="Zhu"){
		maxerror<-1-1/nclases
		}



if (err>=maxerror) {
pesos <- rep(1/n,n)
#c<-0.001
maxerror<-maxerror-eac
c<- log((1-maxerror)/maxerror)

	if (coeflearn=="Breiman"){
	c<- (1/2)*c
	}

	if (coeflearn=="Zhu"){
	c<- c+log(nclases-1)
	}

} 

if (err==0) {
pesos <- rep(1/n,n)
#c<-3
c<- log((1-eac)/eac)

	if (coeflearn=="Breiman"){
	c<- (1/2)*c
	}

	if (coeflearn=="Zhu"){
	c<- c+log(nclases-1)
	}

}

arboles[[m]] <- fit	#Guardamos los arboles
pond[m]<- c 		#Guardamos las ponderaciones

if(m==1){pred <- flearn}
else{pred <- data.frame(pred,flearn)}

if(length(fit$frame$var)>1){
		k <- varImp(fit, surrogates = FALSE, competes = FALSE)
		imp[m,] <-k[sort(row.names(k)), ]
		}
	else {imp[m,]<-rep(0,nvar)} #La solucion I. Medina da problemas con boos=FALSE


}


classfinal <- array(0, c(n,nlevels(vardep)))
for (i in 1:nlevels(vardep)){
 classfinal[,i] <- matrix(as.numeric(pred==levels(vardep)[i]),nrow=n)%*%as.vector(pond)
}

predclass <- rep("O",n)
#2014-11-12 Se puede hacer esto usando apply para evitar el bucle? 
#Creo la funcion "select" que en caso de empate devuelva la clase mayoritaria de entre las empatadas
#predclass[]<-apply(classfinal,1,FUN=select, vardep=vardep)

#2015-07-25 modifico la funcion select para poder usar predict con unlabeled data
predclass[]<-apply(classfinal,1,FUN=select, vardep.summary=summary(vardep))


#for(i in 1:n){	#Cambio esto para resolver los empates
#predclass[i] <- as.character(levels(vardep)[(order(classfinal[i,],decreasing=TRUE)[1])])
#if(length(which(classfinal[i,]==max(classfinal[i,])))>1)
#	{predclass[i] <-names(summary(vardep)[which(classfinal[i,]==max(classfinal[i,]))])[
#order(summary(vardep)[which(classfinal[i,]==max(classfinal[i,]))],decreasing=TRUE)[1]]
#}
#else{predclass[i] <- as.character(levels(vardep)[(order(classfinal[i,],decreasing=TRUE)[1])])} 
#}

#normalizar la importancia de las variables teniendo en cuenta la pond de cada arbol
	imppond<-as.vector(as.vector(pond)%*%imp)
	imppond<-imppond/sum(imppond)*100
	names(imppond)<-sort(row.names(k))





#Para que devuelva las probabilidades a posteriori
classfinal/apply(classfinal,1,sum)->votosporc




ans<- list(formula=formula, trees=arboles, weights=pond, votes=classfinal,prob=votosporc,class=predclass, importance=imppond)

#2015-07-25 pruebo a meter las clases de vardep como atributo de la salida
attr(ans, "vardep.summary") <- summary(vardep, maxsum=700)

mf <- model.frame(formula=formula, data=data) 
terms <- attr(mf, "terms") 
ans$terms <- terms 
ans$call <- match.call()


class(ans)<-"boosting"
ans
}
