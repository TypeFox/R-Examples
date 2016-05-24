"getConToPlot" <- function (C, cohspec, cohcol) 
{
    cohmax <- vector()
	  if(!identical(cohcol, 0)) {	
		  for (i in 1:length(cohcol)) { 

		        cohmax <- append(max(abs( C[,ncol(C)-i+1])), cohmax)
		  	cvec <- C[, ncol(C)-i+1 ] / 
					 max(abs( C[, ncol(C)-i+1 ] ))
			cvec[is.nan(cvec)]<-0 
			C[, ncol(C)-i+1 ] <- cvec
			
	          }       
	  }
	  attr(C, "max") <- cohmax

	  C
}

