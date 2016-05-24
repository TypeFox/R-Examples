prepzscore <-
function (i, j)
{
	i = scale(t(scale(t(i))))
	j = scale(t(scale(t(j))))
	mat = rbind(i,j)
	return(mat)
}

