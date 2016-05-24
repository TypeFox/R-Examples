Source <- function (...) {
	.Deprecated("sourceFormat")
	sourceFormat(...)
}

sourceFormat <- function (file, out.form = getOption("R.output.format"), local = FALSE,
echo = FALSE, print.eval = TRUE, verbose = getOption("verbose"),
prompt.echo = getOption("prompt"), max.deparse.length = 150,
chdir = FALSE, prompt = FALSE)
{

    ## This is a reworked version of .Rsource from RpadUtils (Tom Short)
    ## but this version uses source() itself

    if (is.null(out.form)) out.form <- "text"
    ## capture.all() is inspired from capture.output(), but it captures
    ## both the output and the message streams and it evaluates in .GlobalEnv
    capture.all <- function (...) {
		args <- substitute(list(...))[-1]
    	file <- textConnection("rval", "w", local = TRUE)
        sink(file, type = "output")
        sink(file, type = "message")
        on.exit({
            sink(type = "output")
            sink(type = "message")
            close(file)
        })

    	for (i in seq(length = length(args))) {
            expr <- args[[i]]
            if (mode(expr) == "expression")
				tmp <- lapply(expr, withVisible) #tmp <- lapply(expr, evalVis)
            else if (mode(expr) == "call")
				tmp <- list(withVisible(expr))	 #tmp <- list(evalVis(expr))
			else if (mode(expr) == "name")
            	tmp <- list(withVisible(expr))   #tmp <- list(evalVis(expr))
            else stop("bad argument")
            for (item in tmp) {
				if (item$visible)
				print(item$value)
            }
    	}
    	sink(type = "output")
        sink(type = "message")
		cat("====\n")
    	print(file)
    	cat("====\n")
    	return(file)
    }

    ## We capture output from source() with default args slightly modified
### TODO: get rid of source() and use something like:
    ## (try(parse(textConnection("ls()")), silent = TRUE))
    ## with detection of incomplete lines and other error messages!
    res <- capture.all(source(file = file, local = FALSE, echo = echo,
		print.eval = print.eval, verbose = verbose, prompt.echo = prompt.echo,
		max.deparse.length = max.deparse.length, chdir = chdir))
    if (inherits(res, "list"))
    	res <- paste(res, collapse = "\n")
    if (!out.form %in% c("none", "html"))
        res <- paste(paste(res, collapse="\n"), "\n", sep = "")
    ## Note for out.form == "html", we want to use something like:
    ##require(R2HTML) || stop("Package 'R2HTML' is required!")
    ##res <- HTML(res, file = "")
    ## But since we do not want a dependency to R2HTML here,
    ## we should better put this in the SciViews-R manual
    if (prompt)
    	res <- paste(res, options()$prompt, sep = "")
### TODO: possibly use a continue prompt!
    invisible(res)
}
