identifySubjects <-
function(directory=directory,name.of.log.subjects)
{
	list.subs <-read.csv(paste(directory,name.of.log.subjects,".csv",sep=""))
	subs <- unique(as.numeric(list.subs$id))
	return(subs)
}

