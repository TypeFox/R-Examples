generatetonetwork <-
function(raw, allindivs=union(raw$VertexFrom, raw$VertexTo))
{
	vel <- generatetonetworkfromvel(generatevertexedgelist(raw, allindivs))
	return(vel)
}

