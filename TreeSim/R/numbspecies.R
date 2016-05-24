numbspecies <-
function(tipdepth,branching){
	branchingorder <- order(branching, decreasing=TRUE)
	branchingdepth<-sort(branching,decreasing=TRUE)
	tiporder<-order(tipdepth[,2],decreasing=TRUE)
	tipdepth<-tipdepth[,2][order(tipdepth[,2],decreasing=TRUE)]
	br<-1
	ext<-1
	numbscur <- 1
	event <- vector()
	numbs <- vector()		#numbs - left column node name, middle column time of node, right column number of species after event
	while (tipdepth[ext] > 0 || br <= length(branching)){  # (ext <= length(tipdepth) || br <= length(branching)){        #
		if (br <= length(branching) && branchingdepth[br] > tipdepth[ext]) {
			numbscur <- numbscur + 1
			event <- c(branchingorder[br]+length(tipdepth), branchingdepth[br], numbscur)
			br <- br+1
		} else {
			numbscur <- numbscur - 1
			event <- c(tiporder[ext], tipdepth[ext], numbscur)
			ext <- ext + 1
		}
		numbs <- rbind(numbs,event)
	}
	numbs
	}

