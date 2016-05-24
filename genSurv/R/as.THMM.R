is.THMM <- function(x) {
	return( is.data.frame(x) & inherits(x, "THMM") )
}

as.THMM <- function(x) {
	UseMethod("as.THMM")
}

as.THMM.default <- function(x) {
	stop(gettextf( "cannot coerce class '%s' into class 'THMM'", deparse( class(x) ) ), domain = NA)
}

as.THMM.THMM <- function(x) {
	if ( !is.THMM(x) ) stop("'x' must be of class 'THMM'")
	return(x)
}

as.THMM.CMM <- function(x) {
	if ( !is.CMM(x) ) stop("'x' must be of class 'CMM'")
	data <- cbind(x$id, x$start, x$stop, x$event, x$covariate, x$trans)
	data <- data.frame(data, row.names=NULL)
	names(data) <- c("id", "start", "stop", "event", "covariate", "trans")
	data2 <- matrix(ncol=4, nrow=1)
	i <- 1
	while ( i <= (nrow(data)-1) ) {
		if (data[i,1] == data[i+1,1] & data[i,4] == 0 & data[i+1,4] == 0) {
			aux1 <- c(data[i,1], 0, 1, data[i,5])
			aux2 <- c(data[i,1], data[i,3], 1, data[i,5])
			data2 <- rbind(data2, aux1, aux2)
			i <- i+2
		}
		if (i <= (nrow(data)-1) & data[i,6] == 1 & data[i,4] == 1) {
			aux1 <- c(data[i,1], 0, 1, data[i,5])
			aux2 <- c(data[i,1], data[i,3], 3, data[i,5])
			data2 <- rbind(data2, aux1, aux2)
			i <- i+2
		}
		if (i <= (nrow(data)-1) & data[i,1] == data[i+1,1] & data[i,4] == 0 & data[i+1,4] == 1) {
			if (data[i+2,4] == 0) {
				aux1 <- c(data[i,1], 0, 1, data[i,5])
				aux2 <- c(data[i,1], data[i,3], 2, data[i,5])
				aux3 <- c(data[i,1], data[i+2,3], 2, data[i,5])
				data2 <- rbind(data2, aux1, aux2, aux3)
			}
			if (data[i+2,4] == 1) {
				aux1 <- c(data[i,1], 0, 1, data[i,5])
				aux2 <- c(data[i,1], data[i,3], 2, data[i,5])
				aux3 <- c(data[i,1], data[i+2,3], 3, data[i,5])
				data2 <- rbind(data2, aux1, aux2, aux3)
			}
			i <- i+3
		}
	}
	data2 <- data.frame(data2, row.names=NULL)
	names(data2) <- c("PTNUM", "time", "state", "covariate")
	data2 <- data2[-1,]
	row.names(data2) <- as.integer( 1:nrow(data2) )
	class(data2) <- c(class(data2), "THMM")
	return(data2)
}

as.THMM.TDCM <- function(x) {
	if ( !is.TDCM(x) ) stop("'x' must be of class 'TDCM'")
	data <- cbind(x$start, x$stop, x$event, x$covariate, x$tdcov)
	data <- data.frame(data, row.names=NULL)
	names(data) <- c("start", "stop", "event", "covariate", "tdcov")
	data2 <- data.frame()
	j <- 1
	i <- 1
	doente <- 1
	ultimo <- 0
	while ( i <= nrow(data) ) {
		if ( i == nrow(data) ) ultimo <- 1
		if (ultimo == 0 & data[i,2] == data[i+1,1] & data[i+1,3] == 0) {
			data2[j,1] <- doente
			data2[j,2] <- 0
			data2[j,3] <- 1
			data2[j,4] <- data[i,4]
			j <- j+1
			data2[j,1] <- doente
			data2[j,2] <- data[i,2]
			data2[j,3] <- 2
			data2[j,4] <- data[i,4]
			j <- j+1
			data2[j,1] <- doente
			data2[j,2] <- data[i+1,2]
			data2[j,3] <- 2
			data2[j,4] <- data[i,4]
			i <- i+2
			j <- j+1
			doente <- doente+1
		} else {
			if (ultimo == 0 & data[i,2] == data[i+1,1] & data[i+1,3] == 1) {
				data2[j,1] <- doente
				data2[j,2] <- 0
				data2[j,3] <- 1
				data2[j,4] <- data[i,4]
				j <- j+1
				data2[j,1] <- doente
				data2[j,2] <- data[i,2]
				data2[j,3] <- 2
				data2[j,4] <- data[i,4]
				j <- j+1
				data2[j,1] <- doente
				data2[j,2] <- data[i+1,2]
				data2[j,3] <- 3
				data2[j,4] <- data[i,4]
				i <- i+2
				j <- j+1
				doente <- doente+1
			} else {
				if (data[i,3] == 1) {
					data2[j,1] <- doente
					data2[j,2] <- 0
					data2[j,3] <- 1
					data2[j,4] <- data[i,4]
					j <- j+1
					data2[j,1] <- doente
					data2[j,2] <- data[i,2]
					data2[j,3] <- 3
					data2[j,4] <- data[i,4]
					i <- i+1
					j <- j+1
					doente <- doente+1
				} else {
					data2[j,1] <- doente
					data2[j,2] <- 0
					data2[j,3] <- 1
					data2[j,4] <- data[i,4]
					j <- j+1
					data2[j,1] <- doente
					data2[j,2] <- data[i,2]
					data2[j,3] <- 1
					data2[j,4] <- data[i,4]
					i <- i+1
					j <- j+1
					doente <- doente+1
				}
			}
		}
	}
	names(data2) <- c("PTNUM", "time", "state", "covariate")
	row.names(data2) <- as.integer( 1:nrow(data2) )
	class(data2) <- c(class(data2), "THMM")
	return(data2)
}
