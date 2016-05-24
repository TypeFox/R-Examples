as.scalar <- function(obj){
	return(structure(obj, class=c("scalar",class(obj))));
}


recordlist <- function(mydataframe){
	L <- list();
	if(nrow(mydataframe) == 0) return(L);
	for(i in 1:nrow(mydataframe)){
		L[[i]] <- lapply(as.list(mydataframe[i,]), as.scalar);
	}
	return(L);
}


#' Not a plot, but generates some data with geotags for the demo front-end
#' @param prompt_id an optional prompt id to be added to the output
#' @param ... arguments passed on to oh.survey_response.read
#' @return a list with geotags and optionally pictures
#' @export
gmapdata <- function(prompt_id=NULL, ...){
	
  myData <- oh.survey_response.read(column_list="urn:ohmage:prompt:response,urn:ohmage:context:location:latitude,urn:ohmage:context:location:longitude,urn:ohmage:user:id,urn:ohmage:context:timestamp,urn:ohmage:survey:id", ...);

  #specify output data
  outputcolumns <- c("context.location.longitude","context.location.latitude","context.timestamp","user.id"); 
  newnames <- c("lng", "lat", "timestamp", "user_id");
  
  #check for a photoprompt
  photoprompts <- which(lapply(myData, attr, "prompt_type") == "photo");
  if(length(photoprompts) > 0){
	mostpictures <- which.min(sapply(lapply(myData[photoprompts], is.na),sum));
	promptname <- names(photoprompts)[mostpictures];

    outputcolumns <- c(outputcolumns, promptname);
    newnames <- c(newnames, "photo");
  }
 
  #add a custom column 
  if(!is.null(prompt_id)){
	  outputcolumns <- c(outputcolumns, paste("prompt.id.",prompt_id,sep=""));   
	  newnames <- c(newnames, prompt_id);
 	}
  
  #select rows without missing data
  myData <- myData[outputcolumns];
  #myData <- na.omit(myData);

  #get it in the right format:
  myData[["context.timestamp"]] <- as.character(myData[["context.timestamp"]]);  
  names(myData) <- newnames;

  return(recordlist(myData)); 	
	
}