##############################################################################

#' lists links on a Wikipedia page
#' 
#' @description lists all links (to wikipedia and to external url) that are present in a specific wikipedia page. 
#'  
#' @param page numeric identifier or character title of the specific wikipedia page
#' @param domain a character value specifying the language of the wikipedia url.The default value is "en" for "english language".
#' 
#' @return an object of class \code{linksClass} containing:
#' \itemize{
#' \item{\code{call}}{ the command line}
#' \item{\code{page} }{title, identification number, and domain of the Wikipedia page}
#' \item{\code{links}}{ a data frame containing \code{title}, the titles of the linked pages, 
#'  \code{ns}, numbers of the namespaces (identification of the type of pages, defining nscat et nssubj, the two next variables),  
#'  \code{nscat}, the categories of linked pages (Subject, Talk, or Virtual), 
#'  \code{nssubj}, the subjects of the linked pages (Main Article, User, Wikipedia, File, MediaWiki, Template, Help, Category, Protal, Book, Draft, Education Program, TimedText, Module, Topic, Special, Media, Other). 
#' For more details about namespace, see \url{http://www.mediawiki.org/wiki/Help:Namespaces} and \url{http://en.wikipedia.org/wiki/Wikipedia:Namespace#Subject_namespaces}
#'  }
#'\item{\code{extlinks}}{ a vector containing the list of url of external links. If there is no external link in the page, this item is not created.}

#' \item{\code{testWikiPage}}{a list of four elements,
#'    \itemize{
#'    \item{\code{takeOnlyFirst}}{ a boolean indicating if the class of \code{page} parameter is invalid, 
#'              for example vector, list, matrix..., and in that case, only the first element is considered.}
#'    \item{\code{redirPage}}{ title of the redirected url. This item is NULL if the page is not redirected.}
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
#' }
#'}
#'
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#' 
#' @importFrom XML xmlToList xmlTreeParse htmlParse 
#' @importFrom httr GET
#' @importFrom utils URLencode
#' 
#' @seealso print.linksClass backLinks
#' 
#' @details
#' This function uses the MediaWiki API query syntax: "prop=links". 
#' For more details, see \url{http://www.mediawiki.org/wiki/API:Properties#links_.2F_pl}. 
#' 

#' @examples 
#' \dontrun{
#' # a simple example
#' links("Louis Pasteur") # default domain : en
#' 
#' # with a redirected page
#' links.Obama <- links(page ="Obama")
#' links.Obama
#' # warning message
#' 
#' # a simple example with page specified by its page ID number
#' links(page = 976, domain = "fr" )
#' 
#' # with a page that not exist (at the moment of the redaction of this help page)
#' links("zzzzz") 
#' }
#' @export 
#' 
##############################################################################

links <- function (page=NULL,domain="en") 
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
      url.link  <- GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&titles=",pagebis,"&prop=links&pllimit=max&format=xml", sep=""))
      url.extlink  <- GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&titles=",pagebis,"&prop=extlinks&ellimit=max&format=xml", sep=""))
    }
    else {
      url.link  <- GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&pageids=",pagebis,"&prop=links&pllimit=max&format=xml", sep=""))
      url.extlink  <- GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&pageids=",pagebis,"&prop=extlinks&ellimit=max&format=xml", sep=""))
    }
    
    # XML informations download for the specific URL
    xml.link <- xmlToList(xmlTreeParse(url.link,useInternalNodes=TRUE) )
    list.link <- xml.link$query$pages$page$links
    xml.linkbis <- xml.link
    # Management of the pllimit argument (maximum item per page)
    while (!is.null(xml.linkbis$"query-continue")) {
      continue <-  xml.linkbis$"query-continue"$links      
      url.link <- GET(paste(url.link, "&plcontinue=", continue, sep = "") )
      xml.linkbis  <- xmlToList(xmlTreeParse(url.link,useInternalNodes=TRUE))
      list.link <-c(list.link, xml.linkbis$query$pages$page$links)
    }
    # Specification of Namespace information and the others output 
    ns <- title <- nscat <- nssubj <- NULL
    
    for (i in 1:length(list.link)) {
      
      title <- c(title,iconv(list.link[i]$pl[2],"UTF-8","UTF-8"))
      names(title) <- NULL
      ns <- c(ns,list.link[i]$pl[1])
      
      if (ns[i] == 0) { nscat <- c(nscat,"Subject") ;  nssubj <- c(nssubj,"Main Article")}      
      if (ns[i] == 2) { nscat <- c(nscat,"Subject") ; nssubj <- c(nssubj,"User")}      
      if (ns[i] == 4) { nscat <- c(nscat,"Subject") ; nssubj <- c(nssubj,"Wikipedia")}      
      if (ns[i] == 6) { nscat <- c(nscat,"Subject") ; nssubj <- c(nssubj,  "File")}      
      if (ns[i] == 8) { nscat <- c(nscat,"Subject") ; nssubj <- c(nssubj,  "MediaWiki")}      
      if (ns[i] == 10) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Template")}      
      if (ns[i] == 12) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Help")}      
      if (ns[i] == 14) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Category")}      
      if (ns[i] == 100) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Portal")}      
      if (ns[i] == 108) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Book")}
      if (ns[i] == 118) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Draft")}      
      if (ns[i] == 446) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Education Program")}      
      if (ns[i] == 710) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "TimedText")}      
      if (ns[i] == 828) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Module")}          
      if (ns[i] == 2600) { nscat <- c(nscat,  "Subject") ; nssubj <- c(nssubj,  "Topic")}      
      if (ns[i] == 1) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Talk")}      
      if (ns[i] == 3) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "User talk")}      
      if (ns[i] == 5) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Wikipedia Talk")}      
      if (ns[i] == 7) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "File talk")}      
      if (ns[i] == 9) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "MediaWiki talk")}      
      if (ns[i] == 11) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Template talk")}      
      if (ns[i] == 13) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Help talk")}      
      if (ns[i] == 15) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Category talk")}      
      if (ns[i] == 101) { nscat <- c(nscat,  "Talk") ;  nssubj <- c(nssubj,  "Portal talk")}      
      if (ns[i] == 109) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Book talk")}      
      if (ns[i] == 119) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Draft talk")}      
      if (ns[i] == 447) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Education talk")}      
      if (ns[i] == 711) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "TimedText talk")}      
      if (ns[i] == 829) { nscat <- c(nscat,  "Talk") ; nssubj <- c(nssubj,  "Module talk")}      
      if (ns[i] == -1) { nscat <- c(nscat,  "Virtual") ; nssubj <- c(nssubj,  "Special")}      
      if (ns[i] == -2) { nscat <- c(nscat,  "Virtual") ; nssubj <- c(nssubj,  "Media")}
      if (!(ns[i] %in% c(0,2,4,6,8,10,12,14,100,108,118,446,710,828,2600,1,3,5,7,9,11,13,15,101,109,119,447,711,829))) { nscat <- c(nscat,  "Other") ; nssubj <- c(nssubj,  "Other")}
    }
    ## same for external links
    # XML informations download for the specific URL
    xml.extlink <- xmlToList(xmlTreeParse(url.extlink,useInternalNodes=TRUE) )
    if( any(names(xml.extlink$query$pages$page) =="extlinks") ) 
        {list.extlink <- unlist(xml.extlink$query$pages$page$extlinks)[grep("http",unlist(xml.extlink$query$pages$page$extlinks))] }
    else {list.extlink <- ""}
    xml.extlinkbis <- xml.extlink
    
    # Management of the pllimit argument (maximum item per page)
    while (!is.null(xml.extlinkbis$"query-continue")) {
      continue <-  xml.extlinkbis$"query-continue"$extlinks      
      url.extlink <- paste(url.extlink, "&eloffset=", continue, sep = "") 
      xml.extlinkbis  <- xmlToList(xmlTreeParse(url.extlink,useInternalNodes=TRUE))
      list.extlink <-c(list.extlink, unlist(xml.extlinkbis$query$pages$page$extlinks)[grep("http",unlist(xml.extlinkbis$query$pages$page$extlinks))])
    }
    
    # Management of the output
    
    out$page <- c(iconv(xml.link$query$pages$page$.attrs[3],"UTF-8","UTF-8"),xml.link$query$pages$page$.attrs[1], domain)
    if(length(title) >0 & title[1] !="")
    {
      out$links <- cbind(title,as.numeric(ns),nscat,nssubj)
      rownames(out$links) <- NULL
      colnames(out$links) <- c("title","ns","nscat","nssubj")
      out$links <- as.data.frame(out$links)
    } # end at least one link
    if(length(list.extlink) >0 & !(list.extlink[1] == ""))
    {
      out$links <- as.data.frame(out$links)
      names(list.extlink) <- NULL
      out$extlinks <- list.extlink      
    } # end at least one external link
    
  }
  

  out$testWikiPage <- test
  if(!is.null(test$warnMessage)){for( i in 1:length(test$warnMessage) ) {warning(test$warnMessage[i]) } }
  class(out) <- c("linksClass")
  return (out)
  

}