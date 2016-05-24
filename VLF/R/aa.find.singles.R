aa.find.singles <-
function(aaSpecies, seqlength){
	aa.single <- rep(0,seqlength) # creates a vector of 216 0's
	aa.shared <- rep(0, seqlength) # creates a vector of 216 0's

	#Checks how many specimen are for one species
	for(i in 1:length(aaSpecies)){
		if(nrow(aaSpecies[[i]]) == 1){
			#Counts the number of singleton VLFs
			for(l in 1:seqlength){
				if(is.na(aaSpecies[[i]][l+2]) == FALSE){
					aa.single[l] = aa.single[l] + 1
				}
			}
		}
		else{
			comp <- aa.compare(aaSpecies[[i]], seqlength)
			aa.single <- aa.single + comp[1,]
			aa.shared <- aa.shared + comp[2,]
		}
	}
	return(rbind(aa.single, aa.shared))
}
