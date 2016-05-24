getweight <-
function(item,option,scale=NULL,key=NULL,weights=NULL,NAweight){


	if(option == -1 | is.na(option)){return(NAweight)}
	else if(is.na(option)){return(0)}
	else if(!is.null(weights)){ return(weights[2,which(weights[1,]==option)])}
	else if(scale==1 & option==key){return(1)}
	else if(scale==1 & option!=key){return(0)}
	else if(scale==0){return(option)}
	else if(scale==3){return(0)}

}

