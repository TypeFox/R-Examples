ternary.apply <-
function(grid, f, ...)
{
	f = match.fun(f)
	abc = grid[,letters[1:3]]
	g = function(i) f(abc[i,], ...)
	unlist(lapply(1:nrow(grid), g))
}

