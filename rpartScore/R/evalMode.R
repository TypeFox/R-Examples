evalMode <-
function(y, wt, parms)
{
	class.labs <- as.numeric(names(table(y)))
	node<-xtabs(wt~y)
	id <- class.labs[which(node == max(node))]
    	newid <- ifelse(length(id) > 1, id[sample(1:length(id), size = 1)],id)
	modal.class<-newid
 
	sum.err <- sum((class.labs!=modal.class)*node)
	list(label = modal.class, deviance = sum.err)	
}

