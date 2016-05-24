`print.raggr` <-
function(x, ...){
	k <- length(x$top.list)
	width = options("width")$width
	cat("The optimal list is: ")
	
	res.len <- cumsum(nchar(x$top.list))
	nlines <- ceiling((res.len[k]+k)/width)

	offset <- "\n       "
	offset.len <- nchar(offset)

	cutoffs <- (width-offset.len)*1:nlines	

	cat(paste(offset, paste(x$top.list[res.len < cutoffs[1]], collapse=" ")))
	if(nlines > 1)
		for(i in 2:nlines){
			cat(paste(offset, paste(x$top.list[res.len < cutoffs[i] & 
				res.len > cutoffs[i-1]], collapse=" "), collapse=" "))
		}

	cat("\n\n  Algorithm:  ", x$method)
	cat("\n  Distance:   ", x$distance)
	cat("\n  Score:      ", x$optimal.value, "\n")
   }

