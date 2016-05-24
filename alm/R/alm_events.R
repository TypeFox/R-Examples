#' Retrieve PLoS article-level metrics (ALM) events.
#'
#' Events are the details of the metrics that are counted related to PLoS papers.
#'
#' @importFrom reshape sort_df
#' @importFrom plyr rbind.fill
#' @export
#' @param doi Digital object identifier for an article in PLoS Journals (character)
#' @param pmid PubMed object identifier (numeric)
#' @param pmcid PubMed Central object identifier (numeric)
#' @param wos Web of Science identifier (character)
#' @param scp Scopus identifier (character)
#' @param url Canonical URL (character)
#' @param source_id (character) Name of source to get ALM information for. One source only.
#'    You can get multiple sources via a for loop or lapply-type call.
#' @param publisher_id (character) Metrics for articles by a given publisher, using the Crossref
#'    \code{member_id}.
#' @param compact (logical) Whether to make output compact or not. If TRUE (default), remove
#'    empty sources.
#' @param key (character) Your API key, either enter, or loads from .Rprofile. Only required for
#'    PKP source, not the others.
#' @param api_url API endpoint, defaults to http://alm.plos.org/api/v3/articles (character)
#' @param ... optional additional curl options (debugging tools mostly)
#' @details You can only supply one of the parmeters doi, pmid, pmcid, and mendeley.
#'
#' 		Query for as many articles at a time as you like. Though queries are broken
#' 		up in to smaller bits of 30 identifiers at a time.
#'
#' 		If you supply both the days and months parameters, days takes precedence,
#' 		and months is ignored.
#'
#' 		You can get events from many different sources. After calling alm_events,
#' 		then index the output by the data provider you want. The options are:
#' 		bloglines, citeulike, connotea, crossref, nature, postgenomic, pubmed,
#' 		scopus, plos, researchblogging, biod, webofscience, pmc, facebook,
#' 		mendeley, twitter, wikipedia, and scienceseeker.
#'
#' 		Beware that some data source are not parsed yet, so there may be event data
#' 		but it is not provided yet as it is so messy to parse.
#'
#'    See more info on PLOS's relative metrics event source here
#'    \url{http://www.plosone.org/static/almInfo#relativeMetrics}
#' @return PLoS altmetrics as data.frame's.
#' @references See a tutorial/vignette for alm at
#' \url{http://ropensci.org/tutorials/alm_tutorial.html}
#' @examples \dontrun{
#' # For one article
#' out <- alm_events(doi="10.1371/journal.pone.0029797")
#' names(out) # names of sources
#' # remove those with no data
#' out <- out[!out %in% c("sorry, no events content yet","parser not written yet")]
#' out[["pmc"]] # get the results for PubMed Central
#' out[["twitter"]] # get the results for twitter
#' out[["plos_comments"]] # get the results for PLOS comments, sorta messy
#' out[c("twitter","crossref")] # get the results for two sources
#'
#' # Another example
#' (out <- alm_events(doi="10.1371/journal.pone.0001543"))
#' # remove those with no data
#' out <- out[!out %in% c("sorry, no events content yet","parser not written yet")]
#' names(out)
#' out[['scopus']]
#' out[['mendeley']]
#' out[['figshare']]
#' out[['pubmed']]
#'
#' # Two doi's
#' dois <- c('10.1371/journal.pone.0001543','10.1371/journal.pone.0040117')
#' out <- alm_events(doi=dois)
#' out[[1]]
#' out[[2]]
#' out[[1]][["figshare"]]$events
#'
#' # Many pmcid's
#' out <- alm_events(pmcid=c(212692,2082661))
#' names(out)
#' out['212692']
#'
#' # Many pmid's
#' out <- alm_events(pmid = c(19300479, 19390606, 19343216))
#' names(out)
#' out['19390606']
#'
#' # Specify two specific sources
#' ## You have to do so through lapply, or similar approach
#' lapply(c("crossref","twitter"),
#'    function(x) alm_events(doi="10.1371/journal.pone.0035869", source_id=x))
#'
#' # Figshare data
#' alm_events(doi="10.1371/journal.pone.0069841", source_id='figshare')
#'
#' # Datacite data
#' alm_events("10.1371/journal.pone.0012090", source_id='datacite')
#'
#' # Reddit data
#' alm_events("10.1371/journal.pone.0015552", source_id='reddit')
#'
#' # Wordpress data
#' alm_events("10.1371/journal.pcbi.1000361", source_id='wordpress')
#'
#' # Articlecoverage data
#' alm_events(doi="10.1371/journal.pmed.0020124", source_id='articlecoverage')
#'
#' # Articlecoveragecurated data
#' headfoo <- function(x) head(x$articlecoveragecurated$events)
#' headfoo(alm_events(doi="10.1371/journal.pone.0088278", source_id='articlecoveragecurated'))
#' headfoo(alm_events(doi="10.1371/journal.pmed.1001587", source_id='articlecoveragecurated'))
#'
#' # F1000 Prime data
#' alm_events(doi="10.1371/journal.pbio.1001041", source_id='f1000')
#' dois <- c('10.1371/journal.pmed.0020124','10.1371/journal.pbio.1001041',
#'            '10.1371/journal.pbio.0040020')
#' res <- alm_events(doi = dois, source_id='f1000')
#' res[[3]]
#'
#' # by source_id only
#' alm_events(source_id = "crossref")
#' alm_events(source_id = "reddit")
#'
#' # by publisher_id only
#' alm_events(publisher_id = 340)
#' 
#' # search the software lagotto sever
#' urls <- c("https://github.com/najoshi/sickle","https://github.com/lh3/wgsim",
#'    "https://github.com/jstjohn/SeqPrep")
#' dat <- alm_events(url = urls, api_url = "http://software.lagotto.io/api/v5/articles")
#' }
#'
#' @examples \dontest{
#' # Crossref article data
#' # You need to get an API key first, and pass in a different URL
#' api_url <- "http://alm.labs.crossref.org/api/v3/articles"
#' key <- getOption("crossrefalmkey")
#' # With wikipedia data
#' alm_events(doi='10.1371/journal.pone.0086859', api_url = api_url, key = key)
#' # With facebook data
#' alm(doi='10.1080/15459624.2013.816432', api_url = api_url, key = key)
#' alm_events(doi='10.1080/15459624.2013.816432', url = url, key = key)
#' # With CrossRef citation data - no events data for citations though...
#' alme(doi='10.1021/cr400135x', api_url = api_url, key = key)
#' alm_events(doi='10.1021/cr400135x', api_url = api_url, key = key)
#'
#' # Public Knowledge Project article data
#' # You need to get an API key first, and pass in a different URL
#' api_url <- 'http://pkp-alm.lib.sfu.ca/api/v3/articles'
#' alm_events(doi='10.3402/gha.v7.23554', api_url = api_url, key = getOption("pkpalmkey"))
#'
#' # Copernicus publishers article data
#' # You need to get an API key first, and pass in a different URL
#' api_url <- 'http://metricus.copernicus.org/api/v3/articles'
#' alm_events(doi='10.5194/acpd-14-8287-2014', api_url = api_url, key = getOption("copernicusalmkey"))
#' }

alm_events <- function(doi = NULL, pmid = NULL, pmcid = NULL, wos = NULL, scp = NULL, url = NULL,
  source_id = NULL, publisher_id = NULL, compact = TRUE, key = NULL,
  api_url='http://alm.plos.org/api/v5/articles', ...)
{
	id <- almcompact(list(doi=doi, pmid=pmid, pmcid=pmcid, wos=wos, scp=scp, url=url, source_id=source_id, publisher_id=publisher_id))
	# id <- almcompact(list(doi=doi, pmid=pmid, pmcid=pmcid, wos=wos, scp=scp, url=url))
	if(length(delsp(id)) > 1) {
	  stop("Only supply one of: doi, pmid, pmcid, wos, scp, url")
	}
	if(length(source_id) > 1) stop("You can only supply one source_id")
	if(length(publisher_id) > 1) stop("You can only supply one publisher_id")

	parse_events <- function() {
	  args <- almcompact(list(api_key = key, info = 'detail', source_id = source_id,
                            publisher_id=publisher_id, type = idtype(names(id))))
		if(length(almcompact(list(doi=doi, pmid=pmid, pmcid=pmcid, wos=wos, url=url))) == 0){
      if(length(id) == 0) stop("Please provide one of: doi, pmid, pmcid, wos, scp, url, source_id, or publisher_id")
      ttt <- alm_GET(api_url, args, ...)
      events <- lapply(ttt$data, function(x) x$sources)
		} else {
		  if(length(delsp(id)[[1]]) == 1){
		    passid <- if(names(delsp(id)) == "doi") gsub("/", "%2F", delsp(id)) else delsp(id)
		    ttt <- alm_GET(api_url, c(args, ids = passid[[1]]), ...)
		    events <- lapply(ttt$data, function(x) x$sources)
		  } else
		    if(length(delsp(id)[[1]]) > 1){
		      if(length(delsp(id)[[1]]) > 50){
		        slice <- function(x, n) split(x, as.integer((seq_along(x) - 1) / n))
		        idsplit <- slice(delsp(id)[[1]], 50)
		        repeatit <- function(y) {
		          id2 <- if(names(delsp(id)) == "doi") paste(sapply(y, function(x) gsub("/", "%2F", x)), collapse=",") else paste(delsp(id)[[1]], collapse=",")
		          alm_POST(api_url, c(args, ids = id2), ...)
		        }
		        temp <- lapply(idsplit, repeatit)
		        ttt <- do.call(c, lapply(temp, "[[", "data"))
		        events <- unname(lapply(ttt, function(x) x$sources))
		      } else {
		        id2 <- concat_ids(delsp(id))
		        ttt <- alm_POST(x = api_url, y = c(args, ids = id2), ...)
		        events <- lapply(ttt$data, function(x) x$sources)
		      }
		    }
		}

		# get juse the events data
    # events <- lapply(ttt$data, function(x) x$sources)

		# Function to extract and parse events data for each source
		getevents <- function(x, label=NULL){

			# Parser code
			parsers <- function(y){

        sorry <- "sorry, no events content yet"

				if(y$name == "counter"){
					if(length(y$events)==0){paste(sorry)} else
					{
						year <- as.numeric(sapply(y$events, `[[`, "year"))
            month <- as.numeric(sapply(y$events, `[[`, "month"))
						pdf_views <- as.numeric(sapply(y$events, `[[`, "pdf_views"))
						html_views <- as.numeric(sapply(y$events, `[[`, "html_views"))
						xml_views <- as.numeric(sapply(y$events, `[[`, "xml_views"))
						df <- data.frame(year, month, pdf_views, html_views, xml_views, stringsAsFactors = FALSE)
						list(events_url=y$events_url, events=df, csl=y$events_csl)
					}
				} else if(y$name == "citeulike"){
					if(length(y$events)==0){paste(sorry)} else
					{
            eventspar <- lapply(y$events, data.frame, stringsAsFactors = FALSE)
            csl <- y$events_csl
            list(events_url=y$events_url, events=eventspar, csl=csl)
					}
				} else if(y$name == "crossref"){
					if(length(y$events)==0){paste(sorry)} else
					{
						parsecrossref <- function(x) {
              if(is.null(x[[1]][["publication_type"]])){
                x[[1]][["publication_type"]] <- NA
              }
              if(!("contributors" %in% names(x[[1]]))){
                x[[1]][["contributors"]] <- list(contributor=NA)
                x[[1]]$issn <- paste(x[[1]]$issn, collapse="; ")
								data.frame(x[[1]])
              } else if(length(x[[1]]$contributors$contributor[[1]])>1){
								x[[1]]$contributors$contributor <-
									paste(sapply(x[[1]]$contributors$contributor,
															 function(x) paste(x[1:2], collapse=" ")), collapse="; ")
								x[[1]]$issn <- paste(x[[1]]$issn, collapse="; ")
								data.frame(x[[1]])
							} else {
								x[[1]]$contributors$contributor <-
									paste(x[[1]]$contributors$contributor[1:2], collapse=" ")
								x[[1]]$issn <- paste(x[[1]]$issn, collapse="; ")
								data.frame(x[[1]])
							}
						}
						eventspar <- ldply(y$events, parsecrossref)
            csl <- ldply(y$events_csl, parse_csl)
            list(events_url=y$events_url, events=eventspar, events_csl=csl)
					}
				} else if(y$name == "nature"){
					if(length(y$events)==0){paste(sorry)} else
					{
						parsenature <- function(x){
							temp <- x$event
							blog_ <- data.frame(temp$blog[names(temp$blog) %in% c('title','url')], stringsAsFactors = FALSE)
							names(blog_) <- c('blog_title','blog_url')
							post_ <- data.frame(temp[names(temp) %in% c('title','num_words','url','percent_complex_words','created_at')], stringsAsFactors = FALSE)
							names(post_) <- c('post_percent_complex_words','post_created_at','post_title','post_url','post_num_words')
							cbind(blog_, post_)
						}
						eventspar <- ldply(y$events, parsenature)
						csl <- ldply(y$events_csl, parse_csl)
						list(events_url=y$events_url, events=eventspar, events_csl=csl)
					}
				} else if(y$name == "researchblogging"){
					if(length(y$events)==0){paste(sorry)} else
					{
						parserblogging <- function(w){
							temp <- w$event
              ss <- temp[names(temp) %in% c('post_title','blog_name','blogger_name','published_date','post_url')]
              ss[sapply(ss, is.null)] <- NA
							bloginfo <- data.frame(ss, stringsAsFactors = FALSE)
							if(length(temp$citations$citation[[1]])>1){
								citations <- paste(sapply(temp$citations$citation, function(z) z$doi), sep="", collapse=",")
							} else
							{
								citations <- temp$citations$citation$doi
							}
							cbind(bloginfo, citations)
						}
            df <- do.call(rbind, lapply(y$events, parserblogging))
						list(events_url=y$events_url, events=df, csl=y$events_csl)
					}
				} else if(y$name == "biod"){
					if(length(y$events)==0){paste(sorry)} else
					{
						if(length(y$events) > 1){
							do.call(rbind, lapply(y$events, data.frame))
						} else
						{
							y$events
						}
					}
				} else if(y$name == "pubmed"){
					if(length(y$events)==0){paste(sorry)} else {
					  eventspar <- ldply(y$events, function(x) data.frame(x[c("event","event_url")]))
            list(events_url=y$events_url, events=eventspar, events_csl=y$events_csl)
				  }
				} else if(y$name == "facebook"){
					if(length(y$events)==0){paste(sorry)} else
					{
						parsefb <- function(x){
              x[sapply(x, is.null)] <- "none"
              data.frame(x, stringsAsFactors = FALSE)
						}
            df <- ldply(y$events, parsefb)
						list(events_url=y$events_url, events=df, csl=y$events_csl)
					}
				} else if(y$name == "mendeley"){
					if(length(y$events)==0){paste(sorry)} else
					{
# 						parsemendeley <- function(mm){
# 							readers <- data.frame(name="readers", value=mm$readers, stringsAsFactors = FALSE)
# 							disc <- if(length(mm$discipline) > 1){
#                 ldply(mm$discipline, function(x) data.frame(x, stringsAsFactors = FALSE))[,-1]
# 							} else { data.frame(mm$discipline, stringsAsFactors = FALSE)[,-1] }
# 							country <- ldply(mm$country, function(x) data.frame(x, stringsAsFactors = FALSE))
# 							status <- ldply(mm$status, function(x) data.frame(x, stringsAsFactors = FALSE))
# 							dfs <- list(readers = readers, discipline = disc, country = country, status = status)
# 							ldply(dfs)
# 						}
						# df <- parsemendeley(y$events)
					  list(events_url=y$events_url, events=y$events, csl=y$events_csl)
					}
				} else if(y$name == "twitter"){
					if(length(y$events)==0){paste(sorry)} else
					{
					  temp <- ldply(y$events, function(f){
                data.frame(f$event, event_url=f$event_url, event_time=f$event_time, stringsAsFactors = FALSE)
            })
            csl <- ldply(y$events_csl, parse_csl)
					  list(events_url=y$events_url, events=temp, csl=csl)
					}
				} else if(y$name == "wikipedia"){
					if(length(y$events)==0){paste(sorry)} else
					{
					  df <- ldply(y$events)
					  names(df) <- c("language","values")
					  list(events_url=y$events_url, events=df, csl=y$events_csl)
					}
				} else if(y$name == "bloglines"){
					if(length(y$events)==0){paste(sorry)} else
					{
						parsebloglines <- function(x){
							temp <- data.frame(t(x$event))
							if(any(names(temp) %in% "author")==TRUE && any(names(temp) %in% "site_name")==TRUE)
							{
								temp2 <- temp[,c("site_name","author")]
							} else
							{
								temp2 <- data.frame(site_name=temp$site_name, author="none")
							}
							cbind(temp2, event_url=x$event_url)
						}
						ldply(y$events, parsebloglines)
					}
				} else if(y$name == "postgenomic"){
					if(length(y$events)==0){paste(sorry)} else
						{
							temp <- y$events[[1]]
							name <- temp$event$blog_name
							eventurl <- temp$event_url
							dois <- sapply(temp$event$citing, function(x) x$doi_id )
							list(blog_name=name, event_url=eventurl, dois=dois)
						}
				} else if(y$name == "scopus"){
					if(length(y$events)==0){paste(sorry)} else
						{
						  csl <- ldply(y$events_csl, parse_csl)
						  list(events_url=y$events_url, events=y$events, events_csl=csl)
						}
				} else if(y$name == "wos"){
					if(length(y$events)==0){paste(sorry)} else
					{
						if(length(y$events) > 1){
							ldply(y$events, function(x) data.frame(t(x)))
						} else
						{
							y$events
						}
					}
				} else if(y$name == "pmc"){
					if(length(y$events)==0){paste(sorry)} else
					{
						df <- ldply(y$events, data.frame, stringsAsFactors=FALSE)
						list(events_url=y$events_url, events=df, csl=y$events_csl)
					}
				} else if(y$name == "connotea"){
					if(length(y$events)==0){paste(sorry)} else
					{ paste("parser not written yet") }
				} else if(y$name == "scienceseeker"){
					if(length(y$events)==0){paste(sorry)} else
					{
						 parsesciseeker <- function(xx){
						 	temp <- xx$event
						 	info <- temp[c('title','author')]
						 	recommendations <- data.frame(t(sapply(temp$`ss:community`$`ss:recommendations`, function(x) x[[2]])))
						 	names(recommendations) <- c("user","editor")
						 	categories <- paste(sapply(temp$category, function(x) x[[1]]), collapse=",")

						 	cbind(info, recommendations, categories=categories, event_url=xx$event_url)
						 }
						 df <- ldply(y$events, parsesciseeker)
						 list(events_url=y$events_url, events=df, csl=y$events_csl)
					}
				} else if(y$name == "relativemetric"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    meta <- y$events[names(y$events) %in% c("start_date","end_date")]
				    data <- do.call(rbind.fill,
				                    lapply(y$events$subject_areas, function(x)
				                      data.frame(subject_area=x[[1]], average_usage=t(data.frame(x[[2]])), stringsAsFactors = FALSE)
				                    )
				    )
				    list(events_url=y$events_url, events=list(meta=meta, data=data), csl=y$events_csl)
				  }
				} else if(y$name == "f1000"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    eventsdat <- ldply(y$events, function(b){
				      data.frame(lapply(b$event, function(bb){
				        tmp <- if(length(bb) > 1) paste(do.call(c, bb), collapse = "; ") else bb
				        if(length(tmp) == 1) unlist(tmp) else tmp
				      }), stringsAsFactors=FALSE)
				    })
				    list(events_url=y$events_url, events=eventsdat, csl=y$events_csl)
				  }
				} else if(y$name == "figshare"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
            eventsdat <- lapply(y$events, function(b){
                sapply(b, function(bb){
                  tmp <- if(length(bb) > 1) do.call(c, bb) else bb
                  if(length(tmp) == 1) unlist(tmp) else tmp
                })
              })
				    list(events_url=y$events_url, events=eventsdat, csl=y$events_csl)
				  }
				} else if(y$name == "wordpress"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    eventsdat <- ldply(y$events, function(b){
                data.frame(b$event, event_time=b$event_time, event_url=b$event_url, stringsAsFactors = FALSE)
            })
				    csl <- ldply(y$events_csl, parse_csl)
				    list(events_url=y$events_url, events=eventsdat, csl=csl)
				  }
				} else if(y$name == "pmceurope"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    y$events
				  }
				} else if(y$name == "pmceuropedata"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    list(events_url=y$events_url, events=y$events, csl=y$events_csl)
				  }
				} else if(y$name == "openedition"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    y$events
				  }
				} else if(y$name == "reddit"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
            eventsdat <-
              lapply(y$events, function(b){
                tmp <- b$event[ !names(b$event) %in% c('selftext_html','selftext') ]
                tmp[sapply(tmp, length)==0] <- NA
                list(data=data.frame(tmp, stringsAsFactors = FALSE),
                     details=list(selftext_html=b$event$selftext_html, selftext=b$event$selftext))
              })
            csl <- ldply(y$events_csl, parse_csl)
				    list(events_url=y$events_url, events=eventsdat, csl=csl)
				  }
				}  else if(y$name == "datacite"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
            eventsdat <- ldply(y$events, function(b){
              b$event <- sapply(b$event, function(v) if(length(v) > 1) paste(v, collapse = "; ") else v)
              data.frame(b$event, event_url=b$event_url, stringsAsFactors = FALSE)
            })
				    list(events_url=y$events_url, events=eventsdat, csl=y$events_csl)
				  }
				}  else if(y$name == "copernicus"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    y$events
				  }
				}  else if(y$name == "articlecoverage"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    ev <- do.call(rbind.fill, lapply(y$events, function(y){
              y[sapply(y, is.null)] <- NA
              data.frame(y, stringsAsFactors = FALSE)
            }))
				    csl <- ldply(y$events_csl, parse_csl)
				    list(events_url=y$events_url, events=ev, csl=y$events_csl)
				  }
				}  else if(y$name == "articlecoveragecurated"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    ev <- do.call(rbind.fill, lapply(y$events, function(b){
				      tmp <- b$event
				      tmp[sapply(tmp, length)==0|sapply(tmp, is.null)] <- NA
				      data.frame(tmp,
				                 event_time=b$event_time,
				                 event_url=if(!is.null(b$event_url)) b$event_url else NA,
                         stringsAsFactors = FALSE)
				    }))
				    csl <- ldply(y$events_csl, parse_csl)
				    list(events_url=y$events_url, events=ev, csl=csl)
				  }
				}  else if(y$name == "plos_comments"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
            eventsdat <- lapply(y$events, function(b){
                tmp <- b$event
                tmp[sapply(tmp, length)==0] <- NA
                tmp$replies <- NULL
                list(comment=data.frame(tmp,
                           event_time=b$event_time,
                           event_url=if(!is.null(b$event_url)) b$event_url else NA, stringsAsFactors = FALSE),
                     replies=b$event$replies)
              })
            csl <- ldply(y$events_csl, parse_csl)
				    list(events_url=y$events_url, events=y$events, csl=csl)
				  }
				}  else if(y$name == "twitter_search"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    y$events
				  }
				}  else if(y$name == "doi_resolution"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    y$events
				  }
				}  else if(y$name == "orcid"){
				  if(length(y$events)==0){paste(sorry)} else
				  {
				    y$events
				  }
				}
			}

			# Run the parsers on each element
			datout <- lapply(x, parsers)
			# Assign names to each list element
			names(datout) <- if(is.null(label)) sapply(events[[1]], "[[", "name") else label
			return( datout )
		}

		# Actually get the events data
		tmpout <- lapply(events, getevents, label=source_id)
		# byid <- names(almcompact(list(doi=doi, pmid=pmid, pmcid=pmcid, wos=wos, scp=scp, url=url, source_id=source_id, publisher_id=publisher_id)))
		byid <- names(almcompact(list(doi=doi, pmid=pmid, pmcid=pmcid, wos=wos, scp=scp, url=url)))
		if(length(byid) == 0) byid <- ""
		nmz <- if(byid == "") {
		  NULL
		} else if(!byid == 'doi') {
		  id[[1]]
		} else {
		  if(length(id[[1]]) > 15) {
		    vapply(ttt, "[[", character(1), "id")
		  } else {
		    vapply(ttt$data, "[[", character(1), "id")
		  }
		}
    setNames(tmpout, nmz)
	}
	safe_parse_events <- plyr::failwith(NULL, parse_events)
	finaldata <- safe_parse_events()
	if(length(finaldata) > 1){
	  lapply(finaldata, compact_events)
	} else {
	   compact_events(finaldata[[1]])
	 }
}

compact_events <- function(x){
  phrases <- c("sorry, no events content yet","parser not written yet")
  keep <- sapply(x, function(y) if(any(phrases %in% y)) FALSE else TRUE)
  x[keep]
}

try_date_parts <- function(w){
  tmp <- if(is.null(w[['date-parts']])) w[['date_parts']] else w[['date-parts']]
  paste(unlist(tmp), collapse="-")
}

parse_csl <- function(z){
  z[sapply(z, is.null)] <- NA
  aut <- paste(sapply(z$author, function(zz) paste(zz, collapse = " ")), collapse = "; ")
  data.frame(authors=aut,
             title=z$title,
             container_title=z$`container-title`,
             issued=try_date_parts(z$issued),
             url=z$url,
             type=z$type,
             stringsAsFactors = FALSE)
}

idtype <- function(x){
  x <- x[ !x %in% c("source_id","publisher_id") ]
  if(length(x) == 0){
    NULL
  } else {
    if( x %in% c("doi", "pmid", "pmcid", "wos", "scp", "url") ) x else NULL
  }
}

delsp <- function(x){
  x[ !names(x) %in% c("source_id","publisher_id") ]
}
