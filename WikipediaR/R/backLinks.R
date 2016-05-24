##############################################################################

#' lists pages that link to the Wikipedia page
#' 
#' @description  lists all pages that link to a specific Wikipedia page. 
#' 
#' 
#' @param page numeric identifier or character title of the specific wikipedia page
#' @param domain a character value specifying the language of the wikipedia page.The default value is "en" for "english language".
#' 
#' @return an object of class \code{backLinksClass} containing:
#' \itemize{
#' \item{\code{call}}{ the command line}
#' \item{\code{page} }{ title, identification number and domain of the Wikipedia page, created only if the page exists.}
#'\item{\code{backLinks}}{ a data frame containing 
#'  \code{title}, the titles of the pages that link to the specific page, 
#'  \code{pageid}, the pages IDs,
#'  \code{ns}, numbers of the namespaces (identification of the type of pages, defining nscat et nssubj, the two next variables),  
#'  \code{nscat}, the categories of linked pages (Subject, Talk, or Virtual), 
#'  \code{nssubj}, the subjects of the linked pages (Main Article, User, Wikipedia, File, MediaWiki, Template, Help, Category, Protal, Book, Draft, Education Program, TimedText, Module, Topic, Special, Media, Other). 
#'  For more details about namespace, see \url{http://en.wikipedia.org/wiki/Wikipedia:Namespace#Subject_namespaces}
#' If the page has no back link or does not exist, this item is not created.
#'  }
#' \item{\code{testWikiPage}}{ a list of four elements,
#'    \itemize{
#'    \item{\code{takeOnlyFirst}}{ a boolean indicating if the class of \code{page} parameter is invalid, 
#'              for example vector, list, matrix..., and in that case, only the first element is considered.}
#'    \item{\code{redirPage}}{ title of the redirected page. This item is NULL if the page is not redirected.}
#'    \item{\code{test}}{ an integer with value: 
#'      \itemize{
#'      \item{4}{ for invalid domain,}
#'      \item{3}{ for an empty parameter page,}
#'      \item{2}{ when Wikipedia does not have an article with this exact name,}
#'      \item{1}{ for ambiguous page, direct or redirect,} 
#'      \item{0}{ for valid an unambiguous page, direct or redirect. }
#'      }
#'    }
#'    \item{\code{warnMessage}}{ is a vector of warning messages.}
#'    }
#'  }
#'}
#' @details 
#'  This function uses API query syntax: "list=backlinks". For more details, see \url{https://www.mediawiki.org/wiki/API:Backlinks} 
#' 
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#' 
#' @examples 
#' \dontrun{
#' #' # a simple example
#' backLinks("Louis Pasteur") 
#' 
#' backLinks.Baschet <- backLinks(page ="Cristal Baschet", domain ="fr")
#' table(backLinks.Baschet$backLinks$nscat)
#' 
#' ## example with no back link
#' backLinks(page = 976, domain = "en" )
#' 
#' # with a page that not exist (at the moment of the redaction of this help page)
#' backLinks("zzzzz")
#' }
#' @seealso print.backLinksClass links
#' 
#' @importFrom XML xmlToList xmlTreeParse htmlParse 
#' @importFrom httr GET
#' @importFrom utils URLencode
#' @export 
#' 
##############################################################################

backLinks <- function (page=NULL,domain="en")
  {
  # initialize output 
  out <- NULL
  out$call <- match.call()
  out$testWikiPage <- ""
  
  # first, test validity of the page 
  test <- testWikiPage(domain=domain, page =page)
  # if the domain is valid and the page parameter, and page is not empty :
  if(!(test$test %in% c(2,3,4))) 
  {
    pagebis <-page    
    # if page is a vector, a matrix or a list, take only the first element
    if (test$takeOnlyFirst) { pagebis <- unlist(page)[1]}
    # if the page is redirected, go to the redirected page 
    if(! is.null(test$redirPage)) { pagebis <-test$redirPage }
    
    if (is.character(pagebis)==TRUE)
    {
      # manage encoding and spaces
      pagebis <- gsub(" ",replacement ="_",x = pagebis)
      pagebis <- URLencode(iconv(pagebis,to="UTF-8"))
      url.link  <- GET(paste("http://",domain,
                         ".wikipedia.org/w/api.php?action=query&list=backlinks&blfilterredir=all&bllimit=250&format=xml&bltitle=",
                         pagebis,sep=""))
      url.info  <- GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&titles=",pagebis,"&prop=info&format=xml", sep=""))
      
     }
    else {
      url.link  <- GET(paste("http://",domain,
                         ".wikipedia.org/w/api.php?action=query&list=backlinks&blfilterredir=all&bllimit=250&format=xml&blpageid=",
                         pagebis,sep=""))    
      url.info  <- GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&pageids=",pagebis,"&prop=info&format=xml", sep=""))
      
    }
      
    
    # XML informations download for the specific URL
    # Parses an XML or HTML file or string containing XML/HTML content, and generates an R structure representing the XML/HTML tree
    # XMLInternalNode = TRUE : to separate pageid, ns and title
    xml.link <- xmlToList(xmlTreeParse(url.link,useInternalNodes=TRUE) ) # Convert an XML node/document to a more R-like list
    list.link <- xml.link$query$backlinks
    xml.linkbis  <- xml.link 
    
    # management of limited number of results 
    indice <- 1
    while (!is.null(xml.linkbis$"query-continue") & (indice < 40)) {
      indice <- indice + 1
      continue <-  xml.linkbis$"query-continue"$backlinks  
      url.link <- GET(paste(url.link, "&blcontinue=", continue, sep = "")) 
      xml.linkbis  <- xmlToList(xmlTreeParse(url.link,useInternalNodes=TRUE))
      list.link <-c(list.link, xml.linkbis$query$backlinks)
    }
    if(indice == 40) {warning("More than 10000 pages that links to this page.")}
    #initialize page informations
    gbl.title <- gbl.pageid <- gbl.ns <- nscat <- nssubj <-  NULL

    if (!is.null(list.link)){
    for (i in 1:length(list.link)) {
      if ( length(list.link[i]$bl)== 2){
        for(j in 1: length(list.link[i]$bl$redirlinks)){
          gbl.title <- c(gbl.title,iconv(as.matrix(list.link[i]$bl$redirlinks[j]$bl)[name="title",],"latin1","UTF-8"))
          gbl.pageid <- c(gbl.pageid,as.matrix(list.link[i]$bl$redirlinks[j]$bl)[name="pageid",]) 
          gbl.ns <- c(gbl.ns, as.matrix(list.link[i]$bl$redirlinks[j]$bl)[name="ns",])
        }
      } else {
        gbl.title <- c(gbl.title,iconv(as.matrix(list.link[i]$bl)[name="title",],"UTF-8","UTF-8"))
        gbl.pageid <- c(gbl.pageid,as.matrix(list.link[i]$bl)[name="pageid",]) 
        gbl.ns <- c(gbl.ns,as.matrix(list.link[i]$bl)[name="ns",]) 
      }
      
    } # end loop on lengh od list.link 
    

    for (i in 1 : length(gbl.ns))
    {
      if (gbl.ns[i] == 0) { nscat[i] <- "Subject" ;  nssubj[i] <- "Main Article"}      
      if (gbl.ns[i] == 2) { nscat[i]  <- "Subject" ; nssubj[i] <- "User"}      
      if (gbl.ns[i] == 4) { nscat[i]  <- "Subject" ; nssubj[i] <- "Wikipedia"}      
      if (gbl.ns[i] == 6) { nscat[i]  <- "Subject" ; nssubj[i] <- "File"}      
      if (gbl.ns[i] == 8) { nscat[i]  <- "Subject" ; nssubj[i] <- "MediaWiki"}      
      if (gbl.ns[i] == 10) { nscat[i]  <-   "Subject" ; nssubj[i] <- "Template"}      
      if (gbl.ns[i] == 12) { nscat[i] <-   "Subject" ; nssubj[i] <- "Help"}      
      if (gbl.ns[i] == 14) { nscat[i] <-   "Subject" ; nssubj[i] <- "Category"}      
      if (gbl.ns[i] == 100) { nscat[i]  <-   "Subject" ; nssubj[i] <- "Portal"}      
      if (gbl.ns[i] == 108) { nscat[i]  <-   "Subject" ; nssubj[i] <- "Book"}
      if (gbl.ns[i] == 118) { nscat[i]  <-   "Subject" ; nssubj[i] <- "Draft"}      
      if (gbl.ns[i] == 446) { nscat[i]  <-   "Subject" ; nssubj[i] <- "Education Program"}      
      if (gbl.ns[i] == 710) { nscat[i]  <-   "Subject" ; nssubj[i] <- "TimedText"}      
      if (gbl.ns[i] == 828) { nscat[i]  <-   "Subject" ; nssubj[i] <- "Module"}          
      if (gbl.ns[i] == 2600) { nscat[i]  <-   "Subject" ; nssubj[i] <- "Topic"}      
      if (gbl.ns[i] == 1) { nscat[i] <- "Talk" ; nssubj[i] <- "Talk"}      
      if (gbl.ns[i] == 3) { nscat[i] <- "Talk" ; nssubj[i] <- "User talk"}      
      if (gbl.ns[i] == 5) { nscat[i] <- "Talk" ; nssubj[i] <- "Wikipedia Talk"}      
      if (gbl.ns[i] == 7) { nscat[i] <- "Talk" ; nssubj[i] <- "File talk"}      
      if (gbl.ns[i] == 9) { nscat[i] <- "Talk" ; nssubj[i] <- "MediaWiki talk"}      
      if (gbl.ns[i] == 11) { nscat[i] <- "Talk" ; nssubj[i] <- "Template talk"}      
      if (gbl.ns[i] == 13) { nscat[i] <- "Talk" ; nssubj[i] <- "Help talk"}      
      if (gbl.ns[i] == 15) { nscat[i] <- "Talk" ; nssubj[i] <- "Category talk"}      
      if (gbl.ns[i] == 101) { nscat[i] <- "Talk" ;  nssubj[i] <- "Portal talk"}      
      if (gbl.ns[i] == 109) { nscat[i] <- "Talk" ; nssubj[i] <- "Book talk"}      
      if (gbl.ns[i] == 119) { nscat[i] <- "Talk" ; nssubj[i] <- "Draft talk"}      
      if (gbl.ns[i] == 447) { nscat[i] <- "Talk" ; nssubj[i] <- "Education talk"}      
      if (gbl.ns[i] == 711) { nscat[i] <- "Talk" ; nssubj[i] <- "TimedText talk"}      
      if (gbl.ns[i] == 829) { nscat[i] <- "Talk" ; nssubj[i] <- "Module talk"}      
      if (gbl.ns[i] == -1) { nscat[i] <- "Virtual" ; nssubj[i] <- "Special"}      
      if (gbl.ns[i] == -2) { nscat[i] <- "Virtual" ; nssubj[i] <- "Media"}
      if (!(gbl.ns[i] %in% c(0,2,4,6,8,10,12,14,100,108,118,446,710,828,2600,1,3,5,7,9,11,13,15,101,109,119,447,711,829))) 
      { nscat[i] <- "Other" ; nssubj[i] <- "Other"}
      
    } # end loop on back links
    
    out$backLinks <- cbind(gbl.title, gbl.pageid, gbl.ns, nscat, nssubj)
    colnames(out$backLinks ) <- c("title", "pageid", "ns", "nscat", "nssubj")
    rownames(out$backLinks) <- NULL
    out$backLinks <- as.data.frame(out$backLinks)
    out$backLinks <- unique(out$backLinks)
    
    
    } # end of at least one back link
  xml.info <- xmlToList(xmlTreeParse(url.info,useInternalNodes=TRUE) )
  out$page <- c(iconv(xml.info$query$pages$page["title"],"UTF-8","UTF-8"),xml.info$query$pages$page["pageid"],domain)
    
  } # end if page and domain valid
  
  
  
  
  out$testWikiPage <- test
  if(!is.null(test$warnMessage)){for( i in 1:length(test$warnMessage) ) {warning(test$warnMessage[i]) } }
  class(out) <- c("backLinksClass")
  
  return (out)
  
  

}
