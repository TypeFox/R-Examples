LTT.average.root <- function (trees) { 
furcation <- vector()
branchingtree <- vector()
notultra <- vector()
for (j in 1:length(trees)){
	if (is.ultrametric(trees[[j]]) == TRUE) {
		branching <- LTT.general(trees[[j]])
		branchingdiff <- branching[1,2]
		for (k in 2:(length(branching[,1])-1)){
			branchingdiff <- c(branchingdiff,(branching[k,2]-branching[(k-1),2]))
			}
		furcation <- c(furcation, branchingdiff)
		branchingtemp <- branching[,1]
		branchingtemp <- branchingtemp[-length(branchingtemp)]
		branchingtree <- c(branchingtree, branchingtemp)

	} else {
			notultra <- c(notultra, j)
			}
}
furcation <- furcation[order(branchingtree)]
current <- furcation[1]
linnumber <- current
for (j in 2:(length(furcation))){
	current <- current + furcation[j]
	linnumber <- c(linnumber, current)
	}
branchingtree <- sort(branchingtree)
branchingtree <- c(branchingtree,0)
linnumber <- c(linnumber,linnumber[length(linnumber)])
obj<-cbind(branchingtree, linnumber)
numbevents <- length(obj[,1])
# for (j in 1: (numbevents-2)) {
	# if (obj[numbevents-j,1] == obj[numbevents-(j+1),1]){
		# obj <- obj[- (numbevents-(j+1)),]
		# }
	# }
del<-vector()
for (j in 1: (numbevents-1)) {
	if (abs(obj[j,1] - obj[(j+1),1]) < (10^(-10))){
		del<-c(del,j)
		}
	}
	obj<-obj[-del,]	
for (j in 1:length(obj[,2])){
	obj[j,2]<- obj[j,2]/length(trees)
	}
obj
}
