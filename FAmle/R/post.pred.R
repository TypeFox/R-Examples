post.pred <-
function(z,fun=NULL)
{
	x.n.new <- sapply(as.list(1:z$iter),function(g) distr(z$input$n,z$input$dist,z$sims[g,],'r'))
	if(!is.null(fun)) return(apply(x.n.new,2,fun))
	else return(x.n.new)
}

