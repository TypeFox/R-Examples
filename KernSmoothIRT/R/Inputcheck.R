
Inputcheck<-function(responses,key,format,itemlabels,weights,miss,evalpoints,bandwidth,nitem,nsubj,kernel,NAweight,nevalpoints,thetadist,groups, SubRank){

	
	if(is.na(min(apply(responses,2,sd, na.rm=TRUE)))){
		return("One or more item(s) has no variation in responses (ie. All subjects chose the same option or all responses are missing). Please exclude this item from the responses matrix")
		
	}
	if(min(apply(responses,2,sd, na.rm=TRUE))==0){
		return("One or more item(s) has no variation in responses (ie. All subjects chose the same option). Please exclude this item from the responses matrix")
		
	}
	if(!is.null(weights) & !is.null(format)){
		return('Weights and a format/key were input, use one or the other')
		
	
	}
	if(is.null(weights)){
		if(length(format) != 1){
			if(length(format) != nitem){
			return('format must be a single number or a vector equal to the number of items specifying the format for each.')
			
			}
			if(length(key) != nitem){
			return('format must either be 1 for all Multiple-choice items (with a key specified for each item), 2 for all Rating scale/Partial credit items (with a key specified); 3 for all nominal items (key omitted); a numeric vector of length equal to the number of items specifing the format for each item (key specified); or omitted with the weights argument specifying the weight for each option and item.')
			
		
			}
		}
		else if((format[1] == 1 | format[1] == 2) & length(key) != nitem){
			return('format must either be 1 for all Multiple-choice items (with a key specified for each item), 2 for all Rating scale/Partial credit items (with a key specified); 3 for all nominal items (key omitted); a numeric vector of length equal to the number of items specifing the format for each item (key specified); or omitted with the weights argument specifying the weight for each option and item.')
			
		}
		else if(format[1] == 3 & is.null(SubRank)){
			return('If format = 3 then SubRank must be specified')
			
		}
	}
	if(!is.null(weights)){
		if(class(weights)!="list" | length(weights) != nitem | nrow(weights[[2]])!=2){
		
			return('weights must be a list of length n items. Each list element must have two rows. The first row is the response, the second row is the weight')
	
			}
		
		}

	if(!is.null(key) & length(key)!=nitem ){
		return('key must be of length n items with the correct response for each nominal item or the highest level for each ordinal item')	
		
	}
	if(!is.null(itemlabels) & length(itemlabels)!=nitem){
		return('itemlabels must be of length items or left blank')
		
	}
	if(bandwidth[1] != "Silverman" & bandwidth[1] !="CV" & length(bandwidth) != nitem){
		return('bandwidth much either be: "CV", "default", a vector of length equal to then number of items or left blank')
		
	}
	if(kernel != "gaussian" & kernel != "quadratic" & kernel != "uniform"){
		return('kernel must be either: "gaussian", "quadratic", or "uniform"')
		
	}
	if(NAweight <0 | !is.numeric(NAweight)){
		return('NAweight must be numeric greater than 0');
		
	}
	if(nevalpoints <0 | !is.numeric(nevalpoints)){
		return('nevalpoints must be numeric greater than 0');
		
	}
	if(!is.list(thetadist)){
		return('thetadist must be a list with the first element containing the type of distribution ("norm", "beta", "unif", "cauchy", etc.). The other arguments are the parameters of the distribution ');
		
	}
	if(!is.null(SubRank) & length(SubRank) != nsubj){
		return("SubRank must either be a vector of length equal to the number of subjects with each element specifying the relative rank of each subject or unspecified in which case subjects will be ranked according to the function specified with the RankFun argument")
		
	}



	return("NoErr")


}

