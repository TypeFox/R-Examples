identifyVisits <-
function(directory=directory,name.of.log.subjects)
{
	list.subs <-read.csv(paste(directory,name.of.log.subjects,".csv",sep=""))
	visit <- unique(as.character(list.subs$visit))
	return(visit)
	}

