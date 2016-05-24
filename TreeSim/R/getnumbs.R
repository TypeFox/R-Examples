getnumbs <-
function(tree){
		tipdepth <- age.tips(tree)		
		# entry i is the age of tip i
		branching <- branching.times.complete(tree)
		numbs <- numbspecies(tipdepth,branching)
		 
		#numbs - left column node name, middle column time of node, right column number of species after event
		event <- c(0, 0, 0)   #for last interval in next for loop being ok
		numbs <- rbind(numbs,event)
		numbs
		}

