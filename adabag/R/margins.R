margins <-
function(object, newdata) {


	#newdata must be the same used in object 

	#Exigimos que object contenga votes y formula
	
	if ((is.null(object$votes)==TRUE)|(is.null(object$formula)==TRUE)
	){
	stop("object must have formula and votes in the output")
}



vardep <- newdata[,as.character(object$formula[[2]])]
n <- length(newdata[,1])
nclases <- nlevels(vardep)


votos<- object$votes
votos/apply(votos,1,sum)->votosporc


margen<-rep(0,n) #Creo un vector para guardar los margenes de cada observacion

for (i in 1:n) {
 
k<-votosporc[i, as.numeric(vardep[i])]-votosporc[i,]
margen[i]<- min(k[k!=0])

}

output<- list( margins=margen)
class(output) <- "margins"
output

}

