#' @title Get main content for corpus items, specified by links. 
#' @description \code{getLinkContent} downloads and extracts content from weblinks for \code{\link[tm]{Corpus}} objects.
#' Typically it is integrated and called as a post-processing function (field:\code{$postFUN}) for most \code{\link{WebSource}}
#' objects. \code{getLinkContent} implements content download in chunks which has been proven to be a stabler approach for
#' large content requests. 
#' @param corpus object of class \code{\link[tm]{Corpus}} for which link content should be downloaded
#' @param links character vector specifyinig links to be used for download, defaults to 
#' sapply(corpus, meta, "Origin")
#' @param timeout.request timeout (in seconds) to be used for connections/requests, defaults to 30
#' @param curlOpts curl options to be passed to \code{\link{getURL}}
#' @param chunksize Size of download chunks to be used for parallel retrieval, defaults to 20
#' @param verbose Specifies if retrieval info should be printed, defaults to getOption("verbose")
#' @param retry.empty Specifies number of times empty content sites should be retried, defaults to 3
#' @param sleep.time Sleep time to be used between chunked download, defaults to 3 (seconds)
#' @param extractor Extractor to be used for content extraction, defaults to extractContentDOM
#' @param ... additional parameters to \code{\link{getURL}}
#' @param .encoding encoding to be used for \code{\link{getURL}}, defaults to integer() (=autodetect)
#' @return corpus including downloaded link content
#' @seealso \code{\link{WebSource}} \code{\link[RCurl]{getURL}} \code{\link[boilerpipeR]{Extractor}} 
#' @importFrom NLP content
#' @importFrom RCurl getURL
#' @export
getLinkContent <- function(corpus, links = sapply(corpus, meta, "origin"),
		timeout.request = 30, chunksize = 20, verbose = getOption("verbose"),
		curlOpts = curlOptions(verbose = FALSE,
				followlocation = TRUE, 
				maxconnects = 5,
				maxredirs = 20,
				timeout = timeout.request,
				connecttimeout = timeout.request,
				ssl.verifyhost=FALSE,
				ssl.verifypeer = FALSE,
				useragent = "R", 
				cookiejar = tempfile()),  
		retry.empty = 3, 
		sleep.time = 3, 
		extractor = ArticleExtractor, 
		.encoding = integer(),
		...){
	
	if(length(corpus) != length(links))
		stop("Corpus length not equal to links length\n")
	
	#content_urls <- unlist(sapply(content_parsed, linkreader))
	if(verbose){
		cat("Starting URL Download ...\n")
	}
	retries <- 0
	while(any(empty <- sapply(corpus, function(x) identical(content(x), character(0)))) & (retries <= retry.empty)){
			retries <- retries + 1
			emptycontent.ids <- which(empty)
			
			if(verbose){
				cat("Run ", retries, ", retrieving ", length(emptycontent.ids), " content items\n")
			}
			
			#for(cstart in seq(from = 1, to =  length(links), by = chunksize)){
			for(cstart in seq(from = 1, to =  length(emptycontent.ids), by = chunksize)){
				if(sleep.time > 0){
					if(verbose){
						cat("Sleeping ", sleep.time, " seconds...\n")
					}
					Sys.sleep(sleep.time)
				}
				
				cend <- min(cstart[1] + chunksize-1, length(emptycontent.ids))
				chunk.ids <- emptycontent.ids[cstart:cend]
				chunk <- links[chunk.ids]
				
				# TODO Enable chunk download
				content <- tryCatch({
							getURL(chunk, .opts = curlOpts, .encoding = .encoding, ...)
						},
						error=function(e){
							print(e)
							# TODO: Check if single retrieval part is really necessary
							cat("\nError on retrieval, single retrieval fallback... \n")
							content <- list()
							for(i in 1:length(chunk)){
								content[[i]] <- tryCatch({
											getURL(chunk[i], .opts = curlOpts, .encoding = .encoding, ...)
										},error = function(f) {
											print(f)
											""})
							}
							#cat("Done\n")
							do.call(c, content)})
				
				
				# Extract Content
				extract <- sapply(content, extractor)

				# Put Content Into Corpus
				for(i in 1:length(chunk.ids)){
					cid <- chunk.ids[i]
					content(corpus[[cid]]) <- extract[i]
					
				}
				if(verbose){
					progress <- floor(cend/length(links)*100)
					cat(paste(progress, "% (",cend,"/",length(emptycontent.ids), ") ", Sys.time(), "\n",sep = ""))
				}
			}
	}
	corpus
}