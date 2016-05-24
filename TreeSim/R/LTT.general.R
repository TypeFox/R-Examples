LTT.general <- function (t) {
branchingtree <- branching.times(t)
furcation <- table(t[[1]][,1])-1
furcation <- furcation[order(branchingtree)]
current <- furcation[length(furcation)]+1
linnumber <- current
for (j in 1:(length(furcation)-1)){
	current <- current + furcation[length(furcation)-j]
	linnumber <- c(linnumber, current)
	}
branchingtree <- -sort(branchingtree,decreasing=TRUE)
branchingtree <- c(branchingtree,0)
linnumber <- c(linnumber,linnumber[length(linnumber)])
obj<-cbind(branchingtree, linnumber)
obj
}
