parse.prompt <- function(obj){
	varname <- names(obj);
	prompt_type <- obj[[1]]$context$prompt_type;
	values <- lapply(obj[[1]]$values, parsevector);
	if(prompt_type=="multi_choice" || prompt_type=="multi_choice_custom"){
		values <- identity(values);
	} else {
		values <- unlist(values);
	}
	choice_glossary <- obj[[1]]$context$choice_glossary;
	
	#check if there is a glossary:
	if(!is.null(choice_glossary) && !is.na(choice_glossary)){
		#if there is an explicit 'value' property, that means it is ordered factor
		if(!is.null(choice_glossary[[1]]$value)){
			#if there are explicit "values", then use these
			#levels <- unname(sapply(choice_glossary,unlist)["value",]);
			#labels <- unname(sapply(choice_glossary,unlist)["label",]);
			
			#2.5 update: no longer values as keys. Maybe remove if-else completely.
			levels <- sort(names(choice_glossary));
			labels <- unname(unlist(choice_glossary)[paste(levels,"label",sep=".")]);			
			
		} else {
			#else use the json keys as values
			levels <- sort(names(choice_glossary));
			labels <- unname(unlist(choice_glossary)[paste(levels,"label",sep=".")]);
		}
		
		#in case of ordinal variables, try to sort
		theorder <- order(as.numeric(levels));
		levels <- levels[theorder];
		labels <- labels[theorder];		
		
	}
	
	newvar <- switch(prompt_type,
		single_choice = factor(values, levels, labels, ordered=TRUE),
		#single_choice_custom = factor(values, levels, labels, ordered=TRUE),
		single_choice_custom = factor(values, levels=unique(c(labels, values)), ordered=TRUE),
		multi_choice = multifactor(values, levels, labels),
		#multi_choice_custom = multifactor(values, levels, labels),
		multi_choice_custom = multifactor(values, unique(c(labels, unlist(values)))),
		number = as.numeric(values), #strings are converted to NA without warning
		remote_activity = as.numeric(values),
		timestamp = as.POSIXct(strptime(values, format="%Y-%m-%dT%H:%M:%S")),
		hours_before_now = structure(as.numeric(values), class=c("hours_before_now", "numeric")),
		photo = values,
		text = values,
		#remote_activity = suppressWarnings(sapply(values, mean)), #NOTE: this is a temp solution for the sleep study!!
		stop("Don't know how to parse item: ", varname, "of prompt_type: ", prompt_type, "\n")
	);
	
	attr(newvar, "prompt_type") <- prompt_type;
	return(newvar);
	
}
