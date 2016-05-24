msBP.pdf <-
function(weights, n.points, y=NULL)
{
if(is.null(y))
{
if(n.points<2) stop("n.points must be an integer greater than 2")
grid <- seq(0.001, 0.999, length=n.points)
}
else
{
	grid = y
	n.points = length(y)
}
	w <- tree2vec(weights)
res <- .C("dmsBP_C", as.double(w), as.double(grid), as.integer(n.points), as.integer(weights$max.s), 
	ans = as.double(rep(0, n.points)), PACKAGE = "msBP")
list(y=grid,dens=res$ans)
}

msBP.cdf <-
function(weights, n.points, log=FALSE, y=NULL)
{
	if(is.null(y))
	{
	if(n.points<2) stop("n.points must be an integer greater than 2")
	grid <- seq(0.001, 0.999, length=n.points)
	}
	else
	{
		grid = y
		n.points = length(y)
	}
	w <- tree2vec(weights)
	res <- .C("pmsBP_C", as.double(w), as.double(grid), as.integer(n.points), as.integer(weights$max.s), 
	ans = as.double(rep(0, n.points)), as.integer(log), PACKAGE = "msBP")
	list(y=grid,prob=res$ans)
}
