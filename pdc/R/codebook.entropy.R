codebook.entropy <-
function(data, m, t)
{

	ent <- entropy(codebook(data,m, t))
	
	return( ent);
}
