checkphengen <- function(data)
{
	if (!is(data,"gwaa.data"))
		stop("wrong class of argument")
	if (any(data@phdata$id != data@gtdata@idnames))
		stop("no correspondence @gtdata$idnames <=> @phdata$id")
}
