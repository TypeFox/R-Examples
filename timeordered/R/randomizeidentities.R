randomizeidentities <-
function(raw, withinvertexfrom, byvertexfrom, withreplacement=TRUE)
{
	if (byvertexfrom)
	{
		if (withinvertexfrom)
		{
			list<-by(raw,raw$VertexFrom, function(x) { neworder <- sample(nrow(x), replace= withreplacement); x$VertexTo <- x$VertexTo[neworder]; return(x)})
			outdf <- do.call(rbind.fill, list)	
		}
		else
		{
			neworder <- sample(nrow(raw),replace=withreplacement)
			outdf <- raw
			outdf$VertexTo <- outdf$VertexTo[neworder]
		}
	}
	else
	{
		if (withinvertexfrom)
		{
			list<-by(raw,raw$VertexTo, function(x) { neworder <- sample(nrow(x), replace= withreplacement); x$VertexFrom <- x$VertexFrom[neworder]; return(x)})
			outdf <- do.call(rbind.fill, list)	
		}
		else
		{
			neworder <- sample(nrow(raw), replace= withreplacement)
			outdf <- raw
			outdf$VertexFrom <- outdf$VertexFrom[neworder]
		}		
	}
	
	# remove loops
	outdf$VertexFrom <- as.character(outdf$VertexFrom)
	outdf$VertexTo <- as.character(outdf$VertexTo)

	return(subset(outdf,outdf$VertexFrom!=outdf$VertexTo))
}

