displayPlots <- function(ptList, figureRows, singlePlot, multiPage, bxpt = FALSE){
	nPlot <- length(ptList)
	if(!is.null(singlePlot)) {
		if(singlePlot > length(ptList)) stop(paste("If specified, the \"singlePlot\" argument must be an integer between 1 and ", length(ptList), " for this object."))
   		print(ptList[[singlePlot]])
   		}
   		else if(multiPage){
   			for(i in 1:length(ptList)) print(ptList[[i]])
   			}
   		else{
   			if(is.null(figureRows)){
   				if(bxpt) figureRows <- length(ptList) 
   				else{
      				figureRows <- 1					
   					if(nPlot > 2 & nPlot <= 6) figureRows <- 2
   					if(nPlot > 6) figureRows <- 3
   				}
		}

		figCol <- ceiling(nPlot/figureRows)

		if(dev.cur() == 1) dev.new()

		curCol <- curRow <- 1

		for(i in 1:(nPlot-1)){
			print(ptList[[i]], split = c(curCol,curRow,nx = figCol,ny = figureRows), more = TRUE)
			if(curCol < figCol){
				curCol <- curCol + 1
			}
			else {
				curCol <- 1
				curRow <- curRow + 1
			}
		}

		print(ptList[[nPlot]], split = c(curCol,curRow,nx = figCol,ny = figureRows), more = FALSE)
   		
   		}
	
}