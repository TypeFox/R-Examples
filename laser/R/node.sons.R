`node.sons` <-
function(phy, node)
{
	mat <- phy$edge;
	sons <- mat[,2][which(mat[,1] == node)];
	sons;
}

