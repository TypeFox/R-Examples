##' A driver to parse ascii noweb files with Sweave tool - cacheSweave based
##'
##' @keywords internal
##' @author David Hajage
cacheSweaveAsciiSetup <- function(..., cache = FALSE, trace = FALSE, dependson = NULL) {
  out <- RweaveAsciiSetup(...)

  ## Add the (non-standard) options for code chunks with caching
  out$options[["cache"]] <- cache
	out$options[["dependson"]] <- dependson
	out$options[["trace"]] <- trace
  out
}

##' A driver to parse ascii noweb files with Sweave tool - cacheSweave based
##'
##' @keywords internal
##' @author David Hajage
makeCacheSweaveAsciiCodeRunner <- function(evalFunc = cacheSweave:::cacheSweaveEvalWithOpt) {
  runner <- makeRweaveAsciiCodeRunner(evalFunc)
	function(object, chunk, options) {
		updatedChunk <- FALSE
		e <- runner(object, chunk, options)
		flag <- 'L'
		if(updatedChunk) {
			chunkName <- cacheSweave:::metaChunkName(options)
			cacheSweave:::metaSetCreationTime(chunkName)
			flag <- 'S'
		}
		n <- nchar(as.character(options$chunknr))
    # overwrites the : on preceding row with flag using ANSI
		cat("[F[",n+2,"C",flag,"\n",sep='')
		e
	}
}

cacheSweaveAsciiOptions <- function(options) {
	moreoptions <- c('dependson','cache')
	oldoptions <- options[setdiff(names(options),moreoptions)]
	newoptions <- options[intersect(names(options),moreoptions)]
	Rweaveoptions <- RweaveAsciiOptions(oldoptions)
	options <- unlist(list(Rweaveoptions,newoptions),recursive=FALSE)
  options
}

###############################################

##' A driver to parse asciidoc noweb files with Sweave tool - cacheSweave based
##'
##' @export
##' @rdname RweaveAsciidoc
##' @author David Hajage
cacheSweaveAsciidoc <- function()
{
    require(cacheSweave)
    list(setup = cacheSweaveAsciiSetup,
         runcode = makeCacheSweaveAsciiCodeRunner(),
         writedoc = RweaveAsciiWritedoc,
         finish = RweaveAsciiFinish,
         checkopts = cacheSweaveAsciiOptions)
}

##' A driver to parse txt2tags noweb files with Sweave tool - cacheSweave based
##'
##' @export
##' @rdname RweaveT2t
##' @author David Hajage
cacheSweaveT2t <- function()
{
    list(setup = cacheSweaveT2tSetup,
         runcode = makeCacheSweaveAsciiCodeRunner(),
         writedoc = RweaveAsciiWritedoc,
         finish = RweaveAsciiFinish,
         checkopts = cacheSweaveAsciiOptions)
}

##' A driver to parse txt2tags noweb files with Sweave tool - cacheSweave based
##'
##' @param trace trace
##' @param dependson dependson
##' @keywords internal
##' @rdname RweaveT2t
##' @author David Hajage
cacheSweaveT2tSetup <- function(..., cache = FALSE, trace = FALSE, dependson = NULL) {
  out <- RweaveT2tSetup(...)

  ## Add the (non-standard) options for code chunks with caching
  out$options[["cache"]] <- cache
	out$options[["dependson"]] <- dependson
	out$options[["trace"]] <- trace
  out
}

##' A driver to parse org noweb files with Sweave tool - cacheSweave based
##'
##' @export
##' @rdname RweaveOrg
##' @author David Hajage
cacheSweaveOrg <- function()
{
    list(setup = cacheSweaveOrgSetup,
         runcode = makeCacheSweaveAsciiCodeRunner(),
         writedoc = RweaveAsciiWritedoc,
         finish = RweaveAsciiFinish,
         checkopts = cacheSweaveAsciiOptions)
}

##' A driver to parse org noweb files with Sweave tool - cacheSweave based
##'
##' @param trace trace
##' @param dependson dependson
##' @keywords internal
##' @rdname RweaveOrg
##' @author David Hajage
cacheSweaveOrgSetup <- function(..., cache = FALSE, trace = FALSE, dependson = NULL) {
  out <- RweaveOrgSetup(...)

  ## Add the (non-standard) options for code chunks with caching
  out$options[["cache"]] <- cache
	out$options[["dependson"]] <- dependson
	out$options[["trace"]] <- trace
  out
}

##' A driver to parse pandoc noweb files with Sweave tool - cacheSweave based
##'
##' @export
##' @rdname RweavePandoc
##' @author David Hajage
cacheSweavePandoc <- function()
{
    list(setup = cacheSweavePandocSetup,
         runcode = makeCacheSweaveAsciiCodeRunner(),
         writedoc = RweaveAsciiWritedoc,
         finish = RweaveAsciiFinish,
         checkopts = cacheSweaveAsciiOptions)
}

##' A driver to parse pandoc noweb files with Sweave tool - cacheSweave based
##'
##' @param trace trace
##' @param dependson dependson
##' @keywords internal
##' @rdname RweavePandoc
##' @author David Hajage
cacheSweavePandocSetup <- function(..., cache = FALSE, trace = FALSE, dependson = NULL) {
  out <- RweavePandocSetup(...)

  ## Add the (non-standard) options for code chunks with caching
  out$options[["cache"]] <- cache
	out$options[["dependson"]] <- dependson
	out$options[["trace"]] <- trace
  out
}

##' A driver to parse textile noweb files with Sweave tool - cacheSweave based
##'
##' @export
##' @rdname RweaveTextile
##' @author David Hajage
cacheSweaveTextile <- function()
{
    list(setup = cacheSweaveTextileSetup,
         runcode = makeCacheSweaveAsciiCodeRunner(),
         writedoc = RweaveAsciiWritedoc,
         finish = RweaveAsciiFinish,
         checkopts = cacheSweaveAsciiOptions)
}

##' A driver to parse textile noweb files with Sweave tool - cacheSweave based
##'
##' @param trace trace
##' @param dependson dependson
##' @keywords internal
##' @rdname RweaveTextile
##' @author David Hajage
cacheSweaveTextileSetup <- function(..., cache = FALSE, trace = FALSE, dependson = NULL) {
  out <- RweaveTextileSetup(...)

  ## Add the (non-standard) options for code chunks with caching
  out$options[["cache"]] <- cache
	out$options[["dependson"]] <- dependson
	out$options[["trace"]] <- trace
  out
}

##' A driver to parse rest noweb files with Sweave tool - cacheSweave based
##'
##' @export
##' @rdname RweaveReST
##' @author David Hajage
cacheSweaveReST <- function()
{
    list(setup = cacheSweaveReSTSetup,
         runcode = makeCacheSweaveAsciiCodeRunner(),
         writedoc = RweaveAsciiWritedoc,
         finish = RweaveAsciiFinish,
         checkopts = cacheSweaveAsciiOptions)
}

##' A driver to parse rest noweb files with Sweave tool - cacheSweave based
##'
##' @keywords internal
##' @rdname RweaveReST
##' @author David Hajage
##' @param trace trace
##' @param dependson dependson
cacheSweaveReSTSetup <- function(..., cache = FALSE, trace = FALSE, dependson = NULL) {
  out <- RweaveReSTSetup(...)

  ## Add the (non-standard) options for code chunks with caching
  out$options[["cache"]] <- cache
	out$options[["dependson"]] <- dependson
	out$options[["trace"]] <- trace
  out
}
