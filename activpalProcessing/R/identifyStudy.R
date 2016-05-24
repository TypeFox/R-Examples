identifyStudy <-
function(directory=directory,name.of.log.subjects)
{
	list.subs <-read.csv(paste(directory,name.of.log.subjects,".csv",sep=""))
	study <- unique(as.character(list.subs$study))
	return(study)
}

