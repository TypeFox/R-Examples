is.CMM <- function(x) {
	return( is.data.frame(x) & inherits(x, "CMM") )
}

as.CMM <- function(x) {
	UseMethod("as.CMM")
}

as.CMM.default <- function(x) {
	stop(gettextf( "cannot coerce class '%s' into class 'CMM'", deparse( class(x) ) ), domain = NA)
}

as.CMM.CMM <- function(x) {
	if ( !is.CMM(x) ) stop("'x' must be of class 'CMM'")
	return(x)
}

as.CMM.TDCM <- function(x) {
	if ( !is.TDCM(x) ) stop("'x' must be of class 'TDCM'")
	data <- cbind(x$id, x$start, x$stop, x$event, x$covariate, x$tdcov)
	data <- data.frame(data, row.names=NULL)
	names(data) <- c("id", "start", "stop", "event", "covariate", "tdcov")
	data2 <- matrix(ncol=6, nrow=1)
	i <- 1
	while ( i <= (nrow(data)-1) ) {
		if (data[i,1] != data[i+1,1] & data[i,4] == 0) {
			aux1 <- c(data[i,1], 0, data[i,3], 0, data[i,5], 1)
			aux2 <- c(data[i,1], 0, data[i,3], 0, data[i,5], 2)
			data2 <- rbind(data2, aux1, aux2)
			i <- i+1
		}
		if (i <= (nrow(data)-1) & data[i,1] != data[i+1,1] & data[i,4] == 1) {
			aux1 <- c(data[i,1], 0, data[i,3], 1, data[i,5], 1)
			aux2 <- c(data[i,1], 0, data[i,3], 0, data[i,5], 2)
			data2 <- rbind(data2, aux1, aux2)
			i <- i+1
		}
		if (i <= (nrow(data)-1) & data[i,1] == data[i+1,1] & data[i,4] == 0) {
			if (data[i+1,4] == 0) {
				aux1 <- c(data[i,1], 0, data[i,3], 0, data[i,5], 1)
				aux2 <- c(data[i,1], 0, data[i,3], 1, data[i,5], 2)
				aux3 <- c(data[i,1], data[i,3], data[i+1,3], 0, data[i,5], 3)
				data2 <- rbind(data2, aux1, aux2, aux3)
			} else if (data[i+1,4] == 1) {
				aux1 <- c(data[i,1], 0, data[i,3], 0, data[i,5], 1)
				aux2 <- c(data[i,1], 0, data[i,3], 1, data[i,5], 2)
				aux3 <- c(data[i,1], data[i,3], data[i+1,3], 1, data[i,5], 3)
				data2 <- rbind(data2, aux1, aux2, aux3)
			}
			i <- i+2
		}
	}
	if ( i == nrow(data) ) {
		if (data[i,4] == 0) {
			aux1 <- c(data[i,1], 0, data[i,3], 0, data[i,5], 1)
			aux2 <- c(data[i,1], 0, data[i,3], 0, data[i,5], 2)
			data2 <- rbind(data2, aux1, aux2)
		} else if (data[i,4] == 1) {
			aux1 <- c(data[i,1], 0, data[i,3], 1, data[i,5], 1)
			aux2 <- c(data[i,1], 0, data[i,3], 0, data[i,5], 2)
			data2 <- rbind(data2, aux1, aux2)
		}
	}
	data2 <- data.frame(data2, row.names=NULL)
	names(data2) <- c("id", "start", "stop", "event", "covariate", "trans")
	data2 <- data2[-1,]
	row.names(data2) <- as.integer( 1:nrow(data2) )
	class(data2) <- c(class(data2), "CMM")
	return(data2)
}

as.CMM.THMM <- function(x) {
	if ( !is.THMM(x) ) stop("'x' must be of class 'THMM'")
	data <- cbind(x$PTNUM, x$time, x$state, x$covariate)
	data <- data.frame(data, row.names=NULL)
	names(data) <- c("PTNUM", "time", "state", "covariate")
	data2 <- matrix(ncol=6, nrow=1)
	i <- 1
	while ( i <= (nrow(data)-1) ) {
		if (data[i,2] != 0) i <- i+1
		if (data[i,2] == 0 & data[i+1,3] == 1) {
			aux1 <- c(data[i,1], 0, data[i+1,2], 0, data[i,4], 1)
			aux2 <- c(data[i,1], 0, data[i+1,2], 0, data[i,4], 2)
			data2 <- rbind(data2, aux1, aux2)
		}
		if(data[i,2] == 0 & data[i+1,3] == 3) {
			aux1 <- c(data[i,1], 0, data[i+1,2], 1, data[i,4], 1)
			aux2 <- c(data[i,1], 0, data[i+1,2], 0, data[i,4], 2)
			data2 <- rbind(data2, aux1, aux2)
		}
		if (data[i,2] == 0 & (i+2) <= nrow(data) & data[i+1,3] == 2) {
			if (data[i+2,3] == 2) {
				aux1 <- c(data[i,1], 0, data[i+1,2], 0, data[i,4], 1)
				aux2 <- c(data[i,1], 0, data[i+1,2], 1, data[i,4], 2)
				aux3 <- c(data[i,1], data[i+1,2], data[i+2,2], 0, data[i,4], 3)
				data2 <- rbind(data2, aux1, aux2, aux3)
			}
			if (data[i+2,3] == 3) {
				aux1 <- c(data[i,1], 0, data[i+1,2], 0, data[i,4], 1)
				aux2 <- c(data[i,1], 0, data[i+1,2], 1, data[i,4], 2)
				aux3 <- c(data[i,1], data[i+1,2], data[i+2,2], 1, data[i,4], 3)
				data2 <- rbind(data2, aux1, aux2, aux3)
			}
		}
		i <- i+1
	}
	data2 <- data.frame(data2, row.names=NULL)
	names(data2) <- c("id", "start", "stop", "event", "covariate", "trans")
	data2 <- data2[-1,]
	row.names(data2) <- as.integer( 1:nrow(data2) )
	class(data2) <- c(class(data2), "CMM")
	return(data2)
}
