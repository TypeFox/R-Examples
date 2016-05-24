`grouping.levels` <-
function(model, group)
{
	grouping.original <- model@flist[[group]]
	grouping.names <- levels(grouping.original)
	return(grouping.names)
}

