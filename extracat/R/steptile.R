steptileOLD <- function(x,k = 1, ... ){
	#as.data.frame(as.table(x))
	
	if(!inherits(x,"data.frame")){
			m <- length(dim(x))
			data <- subtable(x,1:m,allfactor=TRUE)
	}else{
		data <- x
		if(!("Freq")%in% names(data)){
			data$Freq <- 1
		}
		m <- ncol(data)-1
		data <- subtable(data,1:m,allfactor=TRUE)
	}


	
	ords <- optile(data[,c(1:(1+k),m+1)], return.data = FALSE, ... )[[3]]
	for(i in 1:(k+1)){
		data[,i] <- factor(data[,i],levels = levels(data[,i])[ords[[i]]])
	}
	
	s <- k+2
	
	while(s < m+1){
		kk <- min(s-1,k)
		ids <- (s-kk):s
		ord <- optile(data[,c(ids,m+1)], perm.cat = c(rep(FALSE,kk),TRUE),return.data = FALSE, ...)[[3]][[kk+1]]
		ords[[s]] <- ord
		data[,s] <- factor(data[,s],levels = levels(data[,s])[ord])
		s <- s+1
		
	}
	if(!inherits(x,"data.frame")){
		x <- xtabs(Freq~.,data=data)
		attr(x,"orders") <- ords
		return(x)
	}else{
		attr(data,"orders") <- ords
		return(data)
	}

}



steptile <- function(x,k = 1, cpcp = FALSE, ... ){
	#as.data.frame(as.table(x))
	if(!cpcp){
		return( steptileOLD(x,k, ... ) )
	}
	
	if(!inherits(x,"data.frame")){
			m <- length(dim(x))
			data <- subtable(x,1:m,allfactor=TRUE)
	}else{
		data <- x
		if(!("Freq")%in% names(data)){
			data$Freq <- 1
		}
		m <- ncol(data)-1
		data <- subtable(data,1:m,allfactor=TRUE)
	}

	
	
	#ords <- optile(data[,c(1:(1+k),m+1)], fun = fun, foreign = foreign, return.data = FALSE, iter = iter)[[3]]
	#for(i in 1:(k+1)){
	#	data[,i] <- factor(data[,i],levels = levels(data[,i])[ords[[i]]])
	#}
	ords <- optile(data[,c(1:2,m+1)], return.data = FALSE, ... )[[3]]
	for(i in 1:2){
		data[,i] <- factor(data[,i],levels = levels(data[,i])[ords[[i]]])
	}

	
	#s <- k+2
	s <- 3
	
	if(k>1){
	
		#o1 <- do.call(order,data[,(k+1):1])
		#data <- data[o1,]
	
		while(s < m+1){
			kk <- min(s-1,k)
			ids <- (s-kk):(s-1)
			
			data$pastry <- do.call(paste,lapply(data[,rev(ids)],function(z) letters[as.integer(z)]))
			o1 <- order(data$pastry)
			data <- data[o1,]
			#print("----------")
			#print(s)
			#print(cbind(data$pastry,data[,1:s]))
			#print(xtabs(Freq~data[,s]+pastry,data=data))
			ord <- optile(data[,c(m+2,s,m+1)], perm.cat = c(FALSE,TRUE),
				return.data = FALSE, ... )[[3]][[2]]
			#print(ord)
			data[,s] <- factor(data[,s],levels = levels(data[,s])[ord])
			ords[[s]] <- ord
			#print(xtabs(Freq~data[,s]+pastry,data=data))
			#print("----------")
			#os <- do.call(order,data[,s:(s-k)])
			#data <- data[os,]
			s <- s+1
		}
		data <- data[,-ncol(data)]
	}else{
		while(s < m+1){
			kk <- min(s-1,k)
			ids <- (s-1):s
			ord <- optile(data[,c(ids,m+1)], perm.cat = c(FALSE,TRUE),return.data = FALSE,  ... )[[3]][[2]]
			ords[[s]] <- ord
			data[,s] <- factor(data[,s],levels = levels(data[,s])[ord])
			s <- s+1
		}
	}

	if(!inherits(x,"data.frame")){
		x <- xtabs(Freq~.,data=data)
		attr(x,"orders") <- ords
		return(x)
	}else{
		attr(data,"orders") <- ords
		return(data)
	}

}