generatelatencies <-
function(raw, allindivs=union(raw$VertexFrom, raw$VertexTo))
{
	maxTime <- max(c(raw$TimeStart,raw$TimeStop))	
	
	vectorClocks = array(-10*maxTime,dim=c(length(allindivs),length(allindivs),maxTime+1),dimnames=list(allindivs,allindivs,0:maxTime))

	timeOrderedEvents <- data.frame(Actor=raw$VertexFrom, Target=raw$VertexTo, Time=raw$TimeStart)
	timeOrderedEvents <- timeOrderedEvents[order(timeOrderedEvents$Time),]

	for (t in 1:(maxTime+1))
	{
		if (t > 1)
		{
			vectorClocks[,,t] <- vectorClocks[,,t-1]
		}
	
		# do the diagonal
		for (i in 1:length(allNames))
		{
			vectorClocks[i,i,t] <- (t-1)	# t=0 at index one
		}
	
		theseEvents <- timeOrderedEvents[timeOrderedEvents$Time == (t-1),]
		if (nrow(theseEvents) > 0)
		{
			for (i in 1:nrow(theseEvents))
			{
				# put in the observed event
				thisEvent <- theseEvents[i,]
				actorID <- as.vector(thisEvent$Actor)
				targetID <- as.vector(thisEvent$Target)
			
				# put in the interaction
				vectorClocks[actorID, targetID, t] <- (t-1)
				# update the clocks
				actorClock <- vectorClocks[actorID, , t]
				targetClock <- vectorClocks[targetID, , t]
				vectorClocks[actorID, , t] <- pmax(actorClock, targetClock)
			

				# the reciprocal interaction
				vectorClocks[targetID, actorID, t] <- (t-1)
				actorClock <- vectorClocks[actorID, , t]
				targetClock <- vectorClocks[targetID, , t]
				vectorClocks[targetID, , t] <- pmax(actorClock, targetClock)
			
			}
		}
	}
	
	vectorClocks[which(vectorClocks < 0)] <- NA
	latencies <- vectorClocks
	for (t in 0:maxTime)
	{
		latencies[,,(t+1)] <- t - vectorClocks[,,(t+1)]
	}
	
	return(latencies)
}

