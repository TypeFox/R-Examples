combine <- 	function(x,y){
			if(x@Number && y@Number){
				Rows <- cbind(x@RowxNumber, y@RowxNumber)
				Cols <- rbind(x@NumberxCol, y@NumberxCol)
				Number <- x@Number + y@Number
				Parameters <- list(x@Parameters, y@Parameters)
			} else if (x@Number && !y@Number){
				Rows <- x@RowxNumber
				Cols <- x@NumberxCol
				Number <- x@Number
			} else if (!x@Number && y@Number){
				Rows <- y@RowxNumber
				Cols <- y@NumberxCol
				Number <- y@Number
			} else
			{
				Rows <- x@RowxNumber
				Cols <- x@NumberxCol
				Number <- x@Number
			}
			Parameters <- list(x@Parameters, y@Parameters)
			Info <-  list(x@info,y@info)
			myBiclust = new("Biclust",Parameters=Parameters, RowxNumber=Rows, NumberxCol=Cols, 
					Number=Number, info=Info)	
	return(myBiclust)
}
		
getStats <- function(x){
	sumRows <- apply(x@RowxNumber,1,sum)
	sumCols <- apply(x@NumberxCol,2,sum)
	return(list(sumRows,sumCols))
}

