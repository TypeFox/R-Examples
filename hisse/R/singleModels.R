CheckFull <- function(x) {
	full=FALSE
	if(min(x)==0) {
		if(length(unique(x))==1+length(sequence(max(x)))) {
			full=TRUE
		}
	} else {
		if(length(unique(x))==length(sequence(max(x)))) {
			full=TRUE
		}	
	}
	return(full)
}



GetAllModels.anc <- function(discrete.states=2, hidden.states=2) {
	combo.length <- discrete.states * hidden.states
	sequence.possibilities <- list(c(0,1)) 
	for (i in 2:combo.length) {
		sequence.possibilities[[i]] <- -1+sequence(i+1)	
	}
	all.combos <- expand.grid(sequence.possibilities)
	good.combos <- apply(all.combos, 1, CheckFull)
	return(all.combos[good.combos,])
}
