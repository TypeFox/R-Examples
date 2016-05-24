randomizetimes <-
function(raw, withinvertexfrom, byvertexfrom, withreplacement=TRUE)
{
	if (byvertexfrom)
	{
		if (withinvertexfrom)
		{
			list<-by(raw,raw$VertexFrom, function(x) { newtimeorder <- sample(nrow(x), replace=withreplacement); x$TimeStart <- x$TimeStart[newtimeorder]; x$TimeStop <- x$TimeStop[newtimeorder]; return(x)})
			outdf <- do.call(rbind.fill, list)	
		}
		else
		{
			newtimeorder <- sample(nrow(raw),replace=withreplacement)
			outdf <- raw
			outdf$TimeStart <- outdf$TimeStart[newtimeorder]
			outdf$TimeStop <- outdf$TimeStop[newtimeorder]
		}
	}
	else
	{
		if (withinvertexfrom)
		{
			list<-by(raw,raw$VertexTo, function(x) { newtimeorder <- sample(nrow(x), replace=withreplacement); x$TimeStart <- x$TimeStart[newtimeorder]; x$TimeStop <- x$TimeStop[newtimeorder]; return(x)})
			outdf <- do.call(rbind.fill, list)	
		}
		else
		{
			newtimeorder <- sample(nrow(raw), replace=withreplacement)
			outdf <- raw
			outdf$TimeStart <- outdf$TimeStart[newtimeorder]
			outdf$TimeStop <- outdf$TimeStop[newtimeorder]
		}
		
	}
	return(outdf)
}

