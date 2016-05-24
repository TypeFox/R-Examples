tstartfun <- function(
	id,
	timevar,
	data)
	{	
		#save input
			tempcall <- match.call()
		#record original order of dataframe so that the output can be returned in the same order
			order.orig <- 1:nrow(data)
			order.orig <- order.orig[order(
				eval(parse(text = paste("data$", deparse(tempcall$id, width.cutoff = 500), sep = ""))),
				eval(parse(text = paste("data$", deparse(tempcall$timevar, width.cutoff = 500), sep = "")))
				)] #sort as below
		#sort dataframe on follow-up time within each individual, necessary for cumulative products below
			data <- data[order(
				eval(parse(text = paste("data$", deparse(tempcall$id, width.cutoff = 500), sep = ""))),
				eval(parse(text = paste("data$", deparse(tempcall$timevar, width.cutoff = 500), sep = "")))
				),]
		#make new dataframe for newly computed variables, to prevent variable name conflicts
			tempdat <- data.frame(
				id = data[,as.character(tempcall$id)],
				timevar = data[,as.character(tempcall$timevar)]
			)
		#compute tstart
			#time of previous record when there is a previous record
			#-1 for first record
			tempfun <- function(x) {
				tempdif <- diff(c(min(x), x))
				tempdif[tempdif == 0] <- min(x) + 1
				return(x - tempdif)
				}
			tempdat$tstart <- unsplit(lapply(split(tempdat$timevar, tempdat$id), function(x)tempfun(x)), tempdat$id)
		#return results in the same order as the original input dataframe
			return(tstart = tempdat$tstart[order(order.orig)])
}
