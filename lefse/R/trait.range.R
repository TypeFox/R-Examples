trait.range <-
function(my.sample, traits){

	range.sub = function(x){

	com.names = names(x[x > 0])
	
	apply(traits[com.names, ], MARGIN = 2, max,na.rm=T) - apply(traits[com.names, ], MARGIN = 2, min,na.rm=T)
	}

apply(my.sample, MARGIN = 1, range.sub)

}
