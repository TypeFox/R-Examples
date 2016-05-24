### Try to source a script file and returns a structured error if it fails
### Romain Francois <francoisromain@free.fr>
trySource <- function (file)
{
	out <- try(source(file), silent = TRUE)
	if (inherits(out, "try-error"))
		out <- parseError(out)
	return(invisible(out))  
}
