find.singles <-
function(species,seqlength){
	single <- rep(0,seqlength) # creates a vector of seqlength 0's
	shared <- rep(0, seqlength) # creates a vector of seqlength 0's
	slength<-seqlength+2
	#Checks how many specimen are for one species
	for(i in 1:length(species)){
		if(nrow(species[[i]]) == 1){
			#Counts the number of singleton VLFs
			for(l in 3:slength){
				if(is.na(species[[i]][l]) == FALSE){
					if(species[[i]][l] == "A"){
						single[l-2] = single[l-2] + 1
					}
					else{if(species[[i]][l] == "C"){
						single[l-2] = single[l-2] + 1
					}
					else{if(species[[i]][l] == "G"){
						single[l-2] = single[l-2] + 1
					}
					else{if(species[[i]][l] == "T"){
						single[l-2] = single[l-2] + 1
				}}}}}
			}
		}
		else{
			comp <- compare(species[[i]], seqlength)
			single <- single + comp[1,]
			shared <- shared + comp[2,]
		}
	}
	return(rbind(single, shared))
}
