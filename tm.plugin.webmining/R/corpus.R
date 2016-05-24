#' @title WebCorpus constructor function.
#' @description \code{WebCorpus} adds further methods and meta data to \code{\link[tm]{Corpus}} and therefore
#' constructs a derived class of \code{\link[tm]{Corpus}}. Most importantly, \code{WebCorpus}
#' calls \code{$PostFUN} on the generated \code{WebCorpus}, which retrieves the main content
#' for most implemented \code{WebSource}s. Thus it enables an efficient retrieval of new feed items
#' (\code{\link{corpus.update}}). All additional WebCorpus fields are added to \code{tm$meta}
#' like \code{$source}, \code{$readerControl} and \code{$postFUN}.
#' @param x object of type Source, see also \code{\link{Corpus}}
#' @param readerControl specifies reader to be used for \code{Source}, defaults to
#' list(reader = x$DefaultReader, language = "en"
#' @param postFUN function to be applied to WebCorpus after web retrieval has been completed,
#' defaults to x$PostFUN
#' @param retryEmpty specifies if retrieval for empty content elements should be repeated, 
#' defaults to TRUE
#' @param ... additional parameters for Corpus function (actually Corpus reader)
#' @importFrom tm Corpus reader getElem stepNext eoi SimpleSource
#' @export
WebCorpus <- function(x, readerControl = list(reader = reader(x), language = "en"),
    postFUN = x$postFUN, retryEmpty = TRUE, ...)
{
  stopifnot(inherits(x, "WebSource"))
  
  readerControl <- prepareReader(readerControl, reader(x))
  
  if (is.function(readerControl$init))
    readerControl$init()
  
  if (is.function(readerControl$exit))
    on.exit(readerControl$exit())
  
  tdl <- vector("list", length(x))
  counter <- 1
  while (!eoi(x)) {
    x <- stepNext(x)
    elem <- getElem(x)
    doc <- readerControl$reader(elem,
        readerControl$language,
        as.character(counter))
    tdl[[counter]] <- doc
    counter <- counter + 1
  }

  corpus <- structure(list(content = tdl,
          meta = CorpusMeta(source = x, readerControl = readerControl, postFUN = postFUN),
          dmeta = data.frame(row.names = seq_along(tdl))),
      class = c("WebCorpus", "VCorpus", "Corpus"))
  if(retryEmpty){
    corpus <- getEmpty(corpus)
  }
  corpus
}

# TODO: Tell Ingo to export CorpusMeta
CorpusMeta <-
    function(..., meta = NULL)
{
  if (is.null(meta))
    meta <- list(...)
  
  stopifnot(is.list(meta))
  
  structure(meta, class = "CorpusMeta")
}

# TODO: Tell Ingo to export prepareReader
prepareReader <- 
function(readerControl, reader = NULL, ...)
{
  if (is.null(readerControl$reader))
    readerControl$reader <- reader
  if (inherits(readerControl$reader, "FunctionGenerator"))
    readerControl$reader <- readerControl$reader(...)
  if (is.null(readerControl$language))
    readerControl$language <- "en"
  readerControl
}


#' @noRd
#' @export
`[.WebCorpus` <- function(x, i) {
	if (missing(i)) return(x)
	corpus <- NextMethod("[")
	class(corpus) <- c("WebCorpus", class(corpus))
	corpus
}

#' @title Update/Extend \code{\link{WebCorpus}} with new feed items.
#' @description The \code{corpus.update} method ensures, that the original 
#' \code{\link{WebCorpus}} feed sources are downloaded and checked against
#' already included \code{TextDocument}s. Based on the \code{ID} included
#' in the  \code{TextDocument}'s meta data, only new feed elements are
#' downloaded and added to the \code{\link{WebCorpus}}.
#' All relevant information regariding the original source feeds are stored
#' in the \code{\link{WebCorpus}}' meta data (\code{\link[tm]{meta}}).
#' @param x object of type \code{\link{WebCorpus}}
#' @param ... 
#' \describe{
#' \item{fieldname}{name of \code{\link{Corpus}} field name to be used as ID, defaults to "ID"}
#' \item{retryempty}{specifies if empty corpus elements should be downloaded again, defaults to TRUE}
#' \item{...}{additional parameters to \code{\link{Corpus}} function}
#' }
#' @export corpus.update
#' @aliases corpus.update.WebCorpus
corpus.update <- function(x, ...){
	UseMethod("corpus.update", x)	
}

#' Update/Extend \code{\link{WebCorpus}} with new feed items.
#' @param x \code{\link{WebCorpus}}
#' @param fieldname name of \code{\link{Corpus}} field name to be used as ID, defaults to "ID"
#' @param retryempty specifies if empty corpus elements should be downloaded again, defaults to TRUE
#' @param ... additional parameters to \code{\link{Corpus}} function
#' @importFrom tm Corpus
#' @importFrom NLP meta
#' @noRd
#' @export
corpus.update.WebCorpus <- 
function(x, fieldname = "id", retryempty = TRUE, verbose = FALSE, ...) {
	cm <- x$meta
	
	newsource <- source.update(cm$source)
	
  #WebCorpus
	newcorpus <- WebCorpus(newsource, readerControl = cm$MetaData$ReaderControl, 
      retryEmpty = FALSE, ...)
	#intersect on ID
	id_old <- sapply(x, meta, fieldname)
	if(any(sapply(id_old, length) == 0))
		stop(paste("Not all elements in corpus to update have field '", fieldname, "' defined", sep = ""))

	id_new <- sapply(newcorpus, meta, fieldname)
	if(any(sapply(id_new, length) == 0))
		stop(paste("Not all elements in corpus to update have field '", fieldname, "' defined", sep = ""))
	
	newcorpus <- newcorpus[!id_new %in% id_old]
	
	if(length(newcorpus) > 0){
		if(!is.null(cm$postFUN)){
			newcorpus <- cm$postFUN(newcorpus)
		}
		corpus <- c(x, newcorpus)
		#attr(corpus, "CMetaData") <- CMetaData(x)
		class(corpus) <- c("WebCorpus", class(corpus))
	}else{
		corpus <- x
	}
	
	if(retryempty){
		corpus <- getEmpty(corpus)
	}
	
	if(verbose){
		cat(length(newcorpus), " corpus items added.\n")
	}
		
	corpus
}


#' @title Retrieve Empty Corpus Elements through \code{$postFUN}. 
#' @description Retrieve content of all empty (textlength equals zero) corpus elements. If 
#' corpus element is empty, \code{$postFUN} is called (specified in \code{\link{meta}})
#' @param x object of type \code{\link{WebCorpus}}
#' @param ... additional parameters to PostFUN
#' @seealso \code{\link{WebCorpus}}
#' @export getEmpty
#' @aliases getEmpty.WebCorpus
getEmpty <- function(x, ...){
	UseMethod("getEmpty", x)	
}
	


#' @importFrom NLP content
#' @noRd
#' @export
getEmpty.WebCorpus <- 
function(x, nChar = 0, ...){
	cm <- x$meta
	noContent <- which(sapply(x, function(y){
            cy <- content(y)
            if(length(cy) == 0L) 0
            else nchar(content(y)) 
          }) <= nChar)
	if(length(noContent) > 0){
		corp_nocontent <- x[noContent]
		if(!is.null(cm$postFUN)){
			corp_nocontent <- cm$postFUN(corp_nocontent, ...)
		}
    # TODO: stupid construct because direct assignment of corpus does not work
    for(i in 1:length(noContent)){
      x[[noContent[i]]] <- corp_nocontent[[i]]
    }
	}
	x
}


