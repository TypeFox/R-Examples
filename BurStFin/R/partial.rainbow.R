"partial.rainbow" <-
function (start=0, end=.35) 
{
	fun.copyright <- "Placed in the public domain 2003-2012 by Burns Statistics Ltd."
	fun.version <- "partial.rainbow 002"

	rainarg <- formals(rainbow)
	rainarg$start <- start
	rainarg$end <- end
	ans <- rainbow
	formals(ans) <- rainarg
	ans
}

