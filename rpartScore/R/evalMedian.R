evalMedian <-
function(y, wt, parms)
{
	class.labs <- as.numeric(names(table(y)))
	node<-xtabs(wt~y)
	cum.node<-cumsum(node)
	{
	if(max(cum.node) %% 2 == "0")
		{
 		median <- c(max(cum.node)/2, (max(cum.node)/2)+1)
		id <- c(min(class.labs[which(median[1] <= cum.node)]), min(class.labs[which(median[2] <= cum.node)]))
		id.freq <- c(node[class.labs == id[1]], node[class.labs == id[2]])
			{
			if(id.freq[1] != id.freq[2])
				{
				median.class <- id[which(id.freq == max(id.freq))]
				}			
			else
				{
				median.class <- id[sample(1:length(id), size = 1)]
				}
			}		
		}
	else
		{
		median <- (max(cum.node)+1)/2
		id <- min(class.labs[which(median <= cum.node)])
		median.class <- id
		}
	}
 
	sum.diff <- sum(abs(class.labs - median.class)*node)
	list(label = median.class, deviance = sum.diff)	
}

