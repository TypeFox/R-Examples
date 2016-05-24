gof.pearson <- function(object){

if(class(object)[1]!="lexpit"&class(object)[1]!="blm")
		stop("Object must be an instance of a blm or lexpit model.")

  Y <- object@y
  prediction <- predict(object)
  
  # UNIQUE COVARIATE GROUPS
  if(class(object)[1]=="blm"){
	  X <- model.matrix(object@formula, object@data)
  }
  else{
  	X <- model.matrix(object@formula.linear, object@data)
  	Z <- model.matrix(object@formula.expit, object@data)
  	X <- cbind(X, Z)
  }
  
   X <-  apply(X, 1, function(x) paste(x, sep="", collapse=""))
   X <- factor(X, labels=1:length(unique(X)))
  
  N <- tapply(Y*object@weights, X, length)
  O <- tapply(Y,X,sum)
  E <- tapply(prediction*object@weights,X,sum)
  pi <- E/N
  
  num <- (O-N*pi)^2
  denom <- N*pi*(1-pi)

  chisq = sum(num/denom)
  
  P = 1 - pchisq(chisq, length(levels(X)) - 2)

  result <- list(E = E, O = O, X2=chisq,p.value=P)
  
  
  print(cbind(E = round(result$E, 2), O = result$O))
  cat("\n","Chi-squared: ",result$chisq,"\n")
  cat("\n","P-value: ",result$p.value,"\n")

 invisible(result)
}

