generatevertexedgelist <-
function(raw, allindivs)
{	
	if (length(allindivs) > length(unique(allindivs)))
	{
		stop("Allindivs is not a list of unique objects - check for duplicates.")	
	}
	
	noninteractors <- setdiff(as.character(allindivs), union(as.character(raw$VertexFrom), as.character(raw$VertexTo)))
	
	DatasetStartTime <- min(c(raw$TimeStart, raw$TimeStop),na.rm=TRUE)
	DatasetStopTime <- max(c(raw$TimeStart, raw$TimeStop),na.rm=TRUE)
	
	startingevents <- data.frame(TimeVertexFrom=paste(raw$VertexFrom, raw$TimeStart), TimeVertexTo=paste(raw$VertexTo, raw$TimeStart),VertexFrom=raw$VertexFrom, VertexTo=raw$VertexTo, TimeStart=raw$TimeStart, TimeStop=raw$TimeStart, TimeCost=0, HopCost=1)
	
	stoppingevents <- data.frame(TimeVertexFrom=paste(raw$VertexFrom, raw$TimeStop), TimeVertexTo=paste(raw$VertexTo, raw$TimeStop), VertexFrom=raw$VertexFrom, VertexTo=raw$VertexTo, TimeStart=raw$TimeStop, TimeStop=raw$TimeStop, TimeCost=0, HopCost=1)
	
	nodepoints <- rbind(data.frame(Vertex=raw$VertexFrom, Time=raw$TimeStart), data.frame(Vertex=raw$VertexFrom, Time=raw$TimeStop), data.frame(Vertex=raw$VertexTo, Time=raw$TimeStart), data.frame(Vertex=raw$VertexTo, Time=raw$TimeStop))
	
	selflist <- by(nodepoints, nodepoints$Vertex, function(x) { timelist <- unique(c(DatasetStartTime, sort(x$Time), DatasetStopTime)); return(data.frame(TimeVertexFrom=paste(x$Vertex[1], head(timelist,-1)), TimeVertexTo=paste(x$Vertex[1], tail(timelist,-1)), VertexFrom=x$Vertex[1], VertexTo=x$Vertex[1], TimeStart=head(timelist,-1), TimeStop=tail(timelist,-1), TimeCost=diff(timelist), HopCost=0   )) })
	
	finaleventlist <- rbind(startingevents, stoppingevents)
	#finalnames <- names(finaleventlist)
	
	selflistall <- do.call("rbind.fill", selflist)
	finaleventlist <- rbind(finaleventlist, selflistall)
	
	
	# also add the non-interactors
	if (length(noninteractors) > 0)
	{
		nonlist <- data.frame(TimeVertexFrom=paste(noninteractors, DatasetStartTime),TimeVertexTo=paste(noninteractors, DatasetStopTime), VertexFrom=noninteractors, VertexTo=noninteractors, TimeStart=DatasetStartTime, TimeStop=DatasetStopTime, TimeCost=(DatasetStopTime-DatasetStartTime), HopCost=0)
	
		finaleventlist <- rbind(finaleventlist, nonlist)
	}
	
	finaleventlist$TimeVertexFrom <- as.vector(finaleventlist$TimeVertexFrom)
	finaleventlist$TimeVertexTo <- as.vector(finaleventlist$TimeVertexTo)
	finaleventlist$VertexFrom <- as.vector(finaleventlist$VertexFrom)
	finaleventlist$VertexTo <- as.vector(finaleventlist$VertexTo)
	
	
	factorlist <- factor(unique(c(finaleventlist$TimeVertexFrom, finaleventlist$TimeVertexTo)))
	finaleventlist$TimeVertexFrom <- factor(as.vector(finaleventlist$TimeVertexFrom), levels=factorlist)
	finaleventlist$TimeVertexTo <- factor(as.vector(finaleventlist$TimeVertexTo), levels=factorlist)

	knownvertexlist <- data.frame(Identity=paste(nodepoints$Vertex, nodepoints$Time),Name=nodepoints$Vertex, Time=nodepoints$Time)
	startvertexlist <- data.frame(Identity=paste(unique(nodepoints$Vertex),DatasetStartTime),Name=unique(nodepoints$Vertex),Time=DatasetStartTime)
	stopvertexlist <- data.frame(Identity=paste(unique(nodepoints$Vertex),DatasetStopTime),Name=unique(nodepoints$Vertex),Time=DatasetStopTime)
	
	tempoutputvertexlist <- rbind(knownvertexlist, startvertexlist, stopvertexlist)
	
	if (length(noninteractors) > 0)
	{
		nonvertexstartlist <- data.frame(Identity=paste(noninteractors, DatasetStartTime),Name=noninteractors,Time=DatasetStartTime)
		nonvertexstoplist <- data.frame(Identity=paste(noninteractors, DatasetStopTime),Name=noninteractors,Time=DatasetStopTime)	
		tempoutputvertexlist <- rbind(tempoutputvertexlist, nonvertexstartlist, nonvertexstoplist)
	}
	
	
	finalvertexlist <- unique(tempoutputvertexlist)

	return(list(edgelist=finaleventlist,vertexlist=finalvertexlist))
}

