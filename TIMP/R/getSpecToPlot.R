"getSpecToPlot" <-
function (E, cohmax, cohcol, plotcohcolspec = TRUE) 
{
	if(plotcohcolspec) {
	 if(!identical(cohcol,0)){
	   E[, cohcol] <- E[, cohcol] * cohmax 
	 }
	}
	else {
	  if(!identical(cohcol,0)){
	     E <- E[, -cohcol]
	  }
	}
	E
}

