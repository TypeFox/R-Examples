separate <-
function(x){
	species <- list()
	p = 1
	i = 1
	final = FALSE
	#A loop to check each specimen
	while(i < (nrow(x))){
		#compares name of specimen to the specimen after it
		if(x[i,2] == x[i+1,2])
		{
			start = i #Records what specimen the matches started at
			count = 2 #A counter for how many specimens match
			if((i+1) == nrow(x)){
				final = TRUE
				last = i + 1
			}
			i = i + 1
			Match = TRUE
			#While the specimen names match check the next	specimen
			while(final == FALSE && Match == TRUE)
			{
				if(x[i,2] == x[i+1,2]){
					count = count + 1 #Counts how many specimen match species names
					i = i + 1
				}
				else{
					Match = FALSE #When specimen names no longer match, sets Match to False
					last = i #Records what specimen matches stopped at
				}
				if((i) == nrow(x)){
				final = TRUE
				last = i
				}
			}
			#Creates a matrix of specimen of the same species
			species[[p]] <- as.matrix(data.frame(x[start:last,]))
			i = i + 1
			p = p + 1
		}		
		else{
			if(final == FALSE){
				#If specimen has no matches, records in list
				species[[p]] = t(as.matrix(x[i,]))
				p = p + 1
				i = i + 1
			}
		}
	}
	if(final == FALSE){
		#Records last specimen in the list
		species[[p]] = t(as.matrix(x[nrow(x),]))
	}
	#Returns list of specimens separated by species
	return(species)
}
