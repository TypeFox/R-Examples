znorm <-
function(m)
{
	m = scale(t(scale(t(m))))
	return(m)
}

