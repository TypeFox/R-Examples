### Guardian API Wrapper v.06
### by M.T. Bastos & C. Puschmann
### library(RCurl)
### library(RJSONIO)

get_json <- function(keywords, format="json", from.date, to.date, api.key)
{
  # pagination
  page.size <- 50
  this.page <- 1
  pages <- 1
  format="json"
  
  if(as.Date(as.character(to.date))-as.Date(as.character(from.date))>28) {
      warning("The requested period is potentially too long. To avoid errors make multiple request (e.g. weekly chunks)")
  }
  
  # prepare list for storing api responses
  api.responses <- NULL
  
  # call guardian API
  while (this.page <= pages)
  {
    request <- paste("http://content.guardianapis.com/search?q=", keywords, "&from-date=", from.date, "&to-date=", to.date, "&format=", format, "&show-fields=all&page=", this.page, "&pageSize=", page.size, "&api-key=", api.key, sep="")
    if(.Platform$OS.type == "windows") { if(!file.exists("cacert.perm")) download.file(url="https://curl.haxx.se/ca/cacert.pem", destfile="cacert.perm") }
    if(.Platform$OS.type == "windows") { json <- getURL(request, cainfo = "cacert.perm", timeout = 240, ssl.verifypeer = FALSE) }
    else { json <- getURL(request, timeout = 240) }
    json <- fromJSON(json, simplify=FALSE)
    this.api.response <- json$response
    stopifnot(!is.null(this.api.response))
    if(this.api.response$total==0){
      print(paste("No matches were found in the Guardian database for keyword '", keywords, "'", sep=""))
      this.page <- this.page + 1
    } else {
    stopifnot(!is.null(this.api.response))
    pages <- this.api.response$pages
    if (pages >= 1)
    {
      print(paste("Fetched page #", this.page, " of ", pages, sep=""))
      api.responses <- c(api.responses, this.api.response)
      this.page <- this.page + 1
    }
    else
    {
      print("Fetched page #1 of 1.")
    }
    api.responses <- c(api.responses, this.api.response)
    this.page <- this.page + 1
  }}
  return(api.responses)
  if(.Platform$OS.type == "windows") { file.remove("cacert.perm") }
}

# parse json to data frame
parse_json_to_df <- function(api.responses)
{
  api.df = data.frame(
    id=NULL,
    sectionId=NULL,
    sectionName=NULL,
    webPublicationDate=NULL,
    webTitle=NULL,
    webUrl=NULL,
    apiUrl=NULL,
    newspaperPageNumber=NULL,
    trailText=NULL,
    headline=NULL,
    showInRelatedContent=NULL,
    lastModified=NULL,
    hasStoryPackage=NULL,
    score=NULL,
    standfirst=NULL,
    shortUrl=NULL,
    wordcount=NULL,
    commentable=NULL,
    allowUgc=NULL,
    isPremoderated=NULL,
    byline=NULL,
    publication=NULL,
    newspaperEditionDate=NULL,
    shouldHideAdverts=NULL,
    liveBloggingNow=NULL,
    commentCloseDate=NULL,
    body=NULL)
  
  # determine number of pages
  pages <- try(seq(from=9, to=length(api.responses), by=9), silent=T)
  
  for (i in pages)
  {
    for (j in 1:length(api.responses[i]$results))
    {
      id <- api.responses[i]$results[[j]]$id
      sectionId <- api.responses[i]$results[[j]]$sectionId
      sectionName <- api.responses[i]$results[[j]]$sectionName
      webPublicationDate <- api.responses[i]$results[[j]]$webPublicationDate
      webTitle <- api.responses[i]$results[[j]]$webTitle
      webUrl <- api.responses[i]$results[[j]]$webUrl
      apiUrl <- api.responses[i]$results[[j]]$apiUrl
      newspaperPageNumber <- api.responses[i]$results[[j]]$fields$newspaperPageNumber
      trailText <- api.responses[i]$results[[j]]$fields$trailText
      headline <- api.responses[i]$results[[j]]$fields$headline
      showInRelatedContent <- api.responses[i]$results[[j]]$fields$showInRelatedContent
      lastModified <- api.responses[i]$results[[j]]$fields$lastModified
      hasStoryPackage <- api.responses[i]$results[[j]]$fields$hasStoryPackage
      score <- api.responses[i]$results[[j]]$fields$score
      standfirst <- api.responses[i]$results[[j]]$fields$standfirst
      shortUrl <- api.responses[i]$results[[j]]$fields$shortUrl
      wordcount <- api.responses[i]$results[[j]]$fields$wordcount
      commentable <- api.responses[i]$results[[j]]$fields$commentable
      allowUgc <- api.responses[i]$results[[j]]$fields$allowUgc
      isPremoderated <- api.responses[i]$results[[j]]$fields$isPremoderated
      byline <- api.responses[i]$results[[j]]$fields$byline
      publication <- api.responses[i]$results[[j]]$fields$publication
      newspaperEditionDate <- api.responses[i]$results[[j]]$fields$newspaperEditionDate
      shouldHideAdverts <- api.responses[i]$results[[j]]$fields$shouldHideAdverts
      liveBloggingNow <- api.responses[i]$results[[j]]$fields$liveBloggingNow
      commentCloseDate <- api.responses[i]$results[[j]]$fields$commentCloseDate
      body <- api.responses[i]$results[[j]]$fields$body
      
      if (is.null(id)) id <- NA
      if (is.null(sectionId)) sectionId <- NA
      if (is.null(sectionName)) sectionName <- NA
      if (is.null(webPublicationDate)) webPublicationDate <- NA
      if (is.null(webTitle)) webTitle <- NA
      if (is.null(webUrl)) webUrl <- NA
      if (is.null(apiUrl)) apiUrl <- NA
      if (is.null(newspaperPageNumber)) newspaperPageNumber <- NA 
      if (is.null(trailText)) trailText <- NA
      if (is.null(headline)) headline <- NA
      if (is.null(showInRelatedContent)) showInRelatedContent <- NA 
      if (is.null(lastModified)) lastModified <- NA
      if (is.null(hasStoryPackage)) hasStoryPackage <- NA
      if (is.null(score)) score <- NA
      if (is.null(standfirst)) standfirst <- NA
      if (is.null(shortUrl)) shortUrl <- NA
      if (is.null(wordcount)) wordcount <- NA
      if (is.null(commentable)) commentable <- NA
      if (is.null(allowUgc)) allowUgc <- NA
      if (is.null(isPremoderated)) isPremoderated <- NA
      if (is.null(byline)) byline <- NA
      if (is.null(publication)) publication <- NA
      if (is.null(newspaperEditionDate)) newspaperEditionDate <- NA
      if (is.null(shouldHideAdverts)) shouldHideAdverts <- NA
      if (is.null(liveBloggingNow)) liveBloggingNow <- NA
      if (is.null(commentCloseDate)) commentCloseDate <- NA
      if (is.null(body)) body <- NA
      
      this.api.df <- data.frame(
        id=id,
        sectionId=sectionId,
        sectionName=sectionName,
        webPublicationDate=webPublicationDate,
        webTitle=webTitle,
        webUrl=webUrl,
        apiUrl=apiUrl,
        newspaperPageNumber=newspaperPageNumber,
        trailText=trailText,
        headline=headline,
        showInRelatedContent=showInRelatedContent,
        lastModified=lastModified,
        hasStoryPackage=hasStoryPackage,
        score=score,
        standfirst=standfirst,
        shortUrl=shortUrl,
        wordcount=wordcount,
        commentable=commentable,
        allowUgc=allowUgc,
        isPremoderated=isPremoderated,
        byline=byline,
        publication=publication,
        newspaperEditionDate=newspaperEditionDate,
        shouldHideAdverts=shouldHideAdverts,
        liveBloggingNow=liveBloggingNow,
        commentCloseDate=commentCloseDate,
        body=body)		
      
      api.df <- rbind(api.df, this.api.df)
    }
  }
  return(unique(api.df))
}

get_guardian <- function(keywords, format="json", from.date, to.date, api.key)
{
  guardian.api.responses <- get_json(keywords, format="json", from.date, to.date, api.key)
  guardian.api.df <- parse_json_to_df(guardian.api.responses)
  return (guardian.api.df)
}
