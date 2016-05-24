p2distance <-
function (matriz, reference_vector = NULL, reference_vector_function = min, iterations = 20, umbral = 0.0001){

	## chequear que la matriz sea de tipo matriz
	## 
	if (!is.matrix(matriz)){
		warning("the argument 'matriz' would be a matrix object")
		matriz <- as.matrix(matriz)
	}		

	discrimination.coefficient <- function (X){
		n <- length(X) # Longitud vector
		X.ord <- sort(X) # Obtener el vector X ordenado 
		vec1 <- rep(1, length(X)) # Creamos un vector de 1 con la misma longitud que X 
		indice <- cumsum(vec1) # Obtenemos un vector con la posiciÛn 
		gini.coeff <- ((2*sum(indice*X.ord))/(n*sum(X)))-((n+1)/n) # coeficiente GINI 
		discrimination.coeff <- 2*gini.coeff*(n/(n-1))
		return(discrimination.coeff)
	} 


	################### funciones auxliares	
	# función que pre-calcula la distancia con el vector de referencia
	calcularDistancia = function (X, vector_referencia = NULL, funcion_v_referencia = min) {	
		if (is.null(vector_referencia)){
			vRef <- makeReferenceVector(X, reference_vector_function = funcion_v_referencia)
		}else{
			vRef <- reference_vector	
		}		
		m <- dim(X)[1] # Calcula el número de filas (m) de la matrix X
		n <- dim(X)[2] # Calcula el número de columnas (n) de la matrix X 
		matRef <- matrix(t(vRef), nrow=m, ncol=n, byrow=TRUE) # Crea una matrix con los valores de referencia para cada variable con el numero de filas igual a la matrix X
		mDif <- X - matRef ### Calcula la matriz de diferencias 
		mDif.abs <- abs(mDif) # matriz de valores absolutos de las diferencias. 
		desT <- as.matrix(apply(X, MARGIN=2, sd)) # Calcula le desviación tipica para cada variable (columna) de la matrix X
		desT <- desT*(sqrt((m-1)/m))
		desT.inversa <- 1/t(desT) # Calcula la inversa de las desviaciones típica
		mdesT <- matrix(t(desT.inversa), nrow=m, ncol=n, byrow=TRUE) #  Crea una matrix con las desviaciones tipicas para cada variable con el numero de filas igual a la matrix X	
		partial.Indicators <- mDif.abs * mdesT ### Crea la matrix Tipificada 
		mI <- as.matrix(partial.Indicators)
		return(partial.Indicators)
	} 
	
	# Devuelve el índice de Frechet de la matriz pasada como referencia
	indiceFrechet = function (matriz){
		iF <- matrix(rowSums(matriz), ncol=1) ### Calcula el Indice Freshet: Suma las filas (municipios) de la matrix Tipificada
		colnames(iF) <- "Frechet.Index" # Dar nombre a la variable Indice de Freshet: IF
		return(iF)
	}
	
	# Ordena las variables de matriz teniendo en cuenta el vector de referencia
	# Devuelve el orden de las variables
	ordenarVariables = function (matriz, referencia){
		colnames(referencia) <- "Reference"
		mVTF <- cbind (matriz, referencia) # Une las dos matrices y obtiene una matriz de m filas y n+1 columnas
		mCor <- cor(mVTF) # Crea la correlación entre cada variable tipificada y el IF o el DP
		vCor <- mCor[colnames(referencia),] # Creo un vector que contenga solamente la correlacion del IF (o DP) con el resto de variables
		vCorAbs <- abs(vCor) # Valores absolutos de la correlación 
		vCorSort <- sort(vCorAbs, decreasing=TRUE) # Ordeno los valores de la correlacion decrecientemente
		vCorSort <- vCorSort[!names(vCorSort) == "Reference"] # Elimino el primer valor que es la correlacion de IF (DP) con IF (DP)
		
		nombres.Ord <- names(vCorSort) # Obtengo los nombres de las variables ordenados segun la correlacion 
		return(nombres.Ord)
	}
	
	# Función que calcula la matriz con los factores de ponderacion
	calculoFactoresPonderacion = function (matriz){
		m <- dim(matriz)[1] # Calcula el número de filas de la matrix X
		n <- dim(matriz)[2] # Calcula el número de columnas de la matrix X 
		vec.results <- numeric() # Creo un  objeto vacío
		for (i in 1:(n-1)) {vec.results[i] <- summary(lm(matriz[,i+1] ~ matriz[,1:i]))$r.squared}
		vect1 <- matrix(1,m,1) # Creo una matriz de 1 columna con el valor 1 y con m filas. 
		coefs <- matrix(vec.results, m, n-1, byrow=TRUE) # Obtengo una matrix con los coeficientes y con m filas
		mR.restado <- 1 - coefs # Calculo 1 - el coeficiente Rsquared
		mFacPond <- cbind(vect1, mR.restado) # Matriz con los factores de ponderacion 
		return(mFacPond)
	}
	
	calculoDP2 = function (matriz, matrizFactores, iteracion = 1){
		mDP <- matriz * matrizFactores 
		DP2 <- t(t(apply(mDP, MARGIN=1, sum)))
		colnames(DP2) <- paste("p2distance",iteracion, sep=".") 
		return(DP2)
	}

	#################################
	#################################	
	#################################
	
	#construimos la lista de resultados
	resultados <- list()
	diff_dps <- numeric()
		
	# calculamos la matriz con las distancias, el vector de refencia el predeterminado (min)
	mDif <- calcularDistancia(matriz, vector_referencia = reference_vector, funcion_v_referencia = reference_vector_function)
	
	#devolvemos también ese valor
	resultados$partial.Indicators = mDif
	
	# calculamos el indice de frechet
	dp2_aux <- indiceFrechet(mDif)	
	dps <- dp2_aux

	iteracion <- 1 #indice de las iteraciones	
	repeat{
		
		print(paste("Iteration", iteracion))
		
		#Ordenamos las variables según la importancia con respecto al vector
		nombres.Ord <- ordenarVariables(mDif, dp2_aux)
	
		#ordenamos la matriz
		mOrdTip <- mDif[,nombres.Ord] ### Reordeno la matriz tipificada segun la correlacion 
	
		#factores de ponderacion
		mFacPond <- calculoFactoresPonderacion(mOrdTip)
	
		dp2_aux <- calculoDP2(mOrdTip, mFacPond, iteracion = iteracion)
		
		#añadimos al vector común
		dps <- cbind(dps, dp2_aux)		
		diff_dps <- cbind(diff_dps, abs(dp2_aux - dps[,iteracion])) # calculamos las diferencias
		umbral_aux <- mean(diff_dps[,iteracion]) # calculamos el umbral alcanzado
		
		#Paramos si se da la condicion de parada del método
		if((iteracion >= iterations) || (umbral >= umbral_aux)) {break}
		
		#si hemos llegado aquí, incrementamos el contador
		iteracion <- iteracion + 1
		
	}

	resultados$discrimination.coefficient <- apply(abs(matriz), 2, discrimination.coefficient)
	resultados$p2distance <- dp2_aux
	resultados$diff_p2distances <- diff_dps
	resultados$p2distances <- dps
	resultados$iteration <- iteracion
	resultados$umbral <- umbral_aux	
	resultados$variables_sort <- nombres.Ord
	resultados$correction_factors <- mFacPond
	colnames(resultados$correction_factors) <- nombres.Ord
	resultados$correction_factors <- resultados$correction_factors[1,]
	resultados$cor.coeff <- cor(matriz, dp2_aux)
	
	# devolver los resultados
	return(resultados)
}

