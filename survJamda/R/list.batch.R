list.batch <-
function(saminfo){
	tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
	batches <- NULL
	for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
	batches
	}

