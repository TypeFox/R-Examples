`select.measurements` <-
function (x){
	
	m.lines <- which(x[[4]][,"sample_type"]=="measurement")
	
	x[[1]] <- x[[1]][m.lines,]
	x[[2]] <- x[[2]][m.lines,]
	x[[4]] <- x[[4]][m.lines,]
	
	return(x)
	
	}

