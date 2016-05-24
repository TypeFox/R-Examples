library(testthat);

#Checks that the two objects in the plot are the same with the exception of the name
#Reports the number of errors along with the types of the differing fields
is_identical_without_name <- function(y,maxLength=5){
	function(x){
		list.x <- unlist(x);
		list.y <- unlist(y);
		raw.x <- list.x[!names(list.x) %in% c("name")];
		raw.y <- list.y[!names(list.y) %in% c("name")];
		
		ret <- isTRUE(all.equal(raw.x,raw.y));
		retStr <- "";#Initialize the return string for later processing if ret isn't true
		
		if(!ret)#If there are differences between them, then print them out
		{
			diffInd <- which(raw.x!=raw.y);
			diffNames <- names(raw.x)[diffInd];#Get the name of the differences
			diffValuesX <- raw.x[diffInd];
			diffValuesY <- raw.y[diffInd];
			
			totalDiff <- length(diffValuesY);
			numericDiff <- length(which(!is.na(as.numeric(diffValuesX))));
			characterDiff <- length(which(is.na(as.numeric(diffValuesX))));
			
			#If there are more than maxLength values to print, only print the first maxLength differences
			if (length(diffInd) > maxLength){
				diffNameStr <- paste0(toString(diffNames[1:maxLength]),"...");
				diffStrX <- paste0(toString(diffValuesX[1:maxLength]),"...");
				diffStrY <- paste0(toString(diffValuesY[1:maxLength]),"...");
			}
			else{
				diffNameStr <- toString(diffNames);
				diffStrX <- toString(diffValuesX);
				diffStrY <- toString(diffValuesY);
			}
			
			retStr <- paste("has different (",diffNameStr,") in",x);
			retStr <- paste(retStr,"\n\tTotal:",totalDiff,"| Numeric:",numericDiff,"| Character:",characterDiff);
			retStr <- paste(retStr,"\n\tThe values are (",diffStrX,") compared to (",diffStrY,")");
		}
		
		#print(retStr);
		
		expectation(
			ret,
			retStr
		)
	}
}
