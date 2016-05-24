`shading` <-
function(breaks,cols=brewer.pal(length(breaks),'Reds')) 
{	res=list(breaks=breaks,cols=cols)
	class(res) = 'shading'
	res }

