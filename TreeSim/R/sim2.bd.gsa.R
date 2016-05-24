sim2.bd.gsa <-
function(n,m,numbtrees,lambda,mu){
	curnumbtrees<-0
	treearray <- list()
	timemax <- 10/(lambda+mu)
	while (curnumbtrees<numbtrees) {
	tree=0
	timeperiods <-vector()
	while(class(tree) != "phylo" || length(timeperiods) == 0){
	tree<-sim2.bd.ssa(m,1,lambda,mu)[[1]]
	if (class(tree)=="phylo") {
		tipdepth <- age.tips(tree)		
		# entry i is the age of tip i
		branching <- branching.times.complete(tree)
		numbs <- numbspecies(tipdepth,branching)
		 
		#numbs - left column node name, middle column time of node, right column number of species after event
		event <- c(0, 0, 0)   #for last interval in next for loop being ok
		numbs <- rbind(numbs,event)
		timeperiods <- which(numbs[,3]==n)
		}
	}
	#class(tree)== "phylo" and length(timeperiods) == 0
		timelength<-vector()
		if (length(timeperiods)>0){
		for (i in 1:length(timeperiods)){
			timelength <- c(timelength, numbs[timeperiods[i],2]-numbs[timeperiods[i]+1,2])
		}	
		}
		timeoverall <- sum(timelength)	
		
		rsamp <- runif(1,min=0,max=1)
		if (rsamp<=(timeoverall/timemax)){
			r <- runif(1,min=0,max=timeoverall)
			timepassed<-0
			j<-0
			while (r>timepassed) {
				j <-j+1
				timepassed<-timepassed+timelength[j]
			}
			cutinterval <- timeperiods[j]
			cuttime <- numbs[cutinterval,2]-runif(1,min=0,max=timelength[j])
	
			if (cuttime>max(branching.times.complete(tree))){
				stop("error")
			}
			treecut <- cuttree(tree,cuttime)
			treearray <- c(treearray,list(reorder(treecut)))
			#treecut<-reorder(treecut)
			curnumbtrees<-curnumbtrees+1
		} else {
				tree=0
		}
	}	
	treearray
	}

