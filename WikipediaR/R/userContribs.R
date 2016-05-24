##############################################################################
#' lists contributions for a specific user


#' @description lists contributions for a specific user: discussions for Talk pages, and revisions for Subject pages.
#' 
#' @param user.name a character value providing the name of the user
#' @param domain a character value providing the language of the wikipedia page. The default value is "en" for "english language".
#' @param ucprop Include informations (separate with '|') :
#' \itemize{
#' \item{ids}{ the page ID and revision ID}
#' \item{title}{ the title and namespace ID of the page. This parameter is not optional.}
#' \item{timestamp}{ the timestamp of the edit, i.e. day and time}
#' \item{comment}{ the comment of the edit}
#' \item{parsedcomment}{ the parsed comment of the edit}
#' \item{size}{ the new size of the edit}
#' \item{sizediff}{ the size delta of the edit against its parent}
#' \item{flags}{ flags of the edit}
#' \item{patrolled}{ patrolled edits}
#' \item{tags}{ tags for the edit}
#' }
#' Default: ids|title|timestamp|comment|sizediff|flags

#' @return an object of class \code{userContribsClass}:
#' \itemize{
#' \item{\code{call}}{ the command line}
#' \item{\code{user}}{ a vector containing the user name and the user identifier}
#' \item{\code{contribs}}{ a data frame containing the asked properties of the revisions. 
#' If ids is asked, the corresponding results are pageid, revid, and parentid. 
#' Three others informations are automatically added : 
#'  \code{ns}, numbers of the namespaces (identification of the type of pages, defining nscat et nssubj, the two next variables), 
#'  \code{nscat}, the categories of linked pages (Subject, Talk, or Virtual), 
#'  \code{nssubj}, the subjects of the linked pages (Main Article, User, Wikipedia, File, MediaWiki, Template, Help, Category, Protal, Book, Draft, Education Program, TimedText, Module, Topic, Special, Media, Other) . 
#' This item is provided if the user has at least one contribution.
#' For more details about namespace, see \url{http://www.mediawiki.org/wiki/Help:Namespaces} and \url{http://en.wikipedia.org/wiki/Wikipedia:Namespace#Subject_namespaces}

#'   }
#' \item{\code{testWikiUser}}{ A list of three elements. 
#'  The first is \code{takeOnlyFirst}, a boolean indicating if the class of \code{user.name} parameter is invalid, 
#'  for example vector, list, matrix..., and in that case, the only the first element is considered.
#'    The second element is \code{test}, an integer with value: 
#'    \itemize{
#'    \item{4}{ for invalid domain,}
#'    \item{3}{ for an empty parameter user,}
#'    \item{2}{ when Wikipedia does not have an user with this exact name,}
#'    \item{0}{ for valid existing user. }
#'    }
#'  The last element, \code{warnMessage}, is a vector of warning messages.
#'  }
#' }
#' @export
#' 
#' @details 
#' This function uses the API query syntax: "list=usercontribs". 
#' For more details, see \url{https://www.mediawiki.org/wiki/API:Usercontribs}
#' Additionnally to the titles of the modified pages, this function always returns in the contribs item the ns for namespace (identifcation of the type of pages), 
#' the nscat for the category of the pages (Subject, Talk, or Virtual), and the subject of the pages.
#' 
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#'
#' @seealso userInfo
#' 
#' @importFrom XML xmlToList xmlTreeParse htmlParse 
#' @importFrom httr GET
#' @examples 
#' \dontrun{
#' LouiseContribs <- userContribs(user.name = "Louise", domain = "en")
#' ## try a user that does not exist (at the moment of the redaction of this help page)
#' userContribs(user.name="Louise Baschet", domain ="fr")
#' }
#############################################################################

userContribs <- function (user.name=NULL,domain="en", ucprop ="ids|title|timestamp|comment|sizediff|flags" )
  {
  
  
  # initialize output 
  out <- NULL
  out$call <- match.call()
  ## verify ucprop parameters
  if (is.null(ucprop)) { stop("ucprop argument is required.")}
  if(!is.character(ucprop)){stop("ucprop argument must be a character string")}
  
  props <- unlist(strsplit(ucprop, split ="|", fixed = TRUE))
  
  if (! all(props %in% c("ids", "title", "timestamp", "comment", "parsedcomment","size","sizediff", "flags", "patrolled", "tags")))
  {stop("ucprop argument is not valid")}
  if (! any(props =="title"))   {stop("title is not optional in ucprop argument.")}
  if (any(props == "ids" )) { props <- c(props [-which(props == "ids")],"pageid","revid","parentid")}
  
  ## verify the domain and user.name validity
  test <- testWikiUser(user.name = user.name, domain = domain)

  if (test$test ==0) {

    if (test$takeOnlyFirst) { user.name <- unlist(user.name )[1]}
 
    # URL building
    
    url.contrib = GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&list=usercontribs&ucprop=",ucprop,
                        "&ucuser=", user.name, "&uclimit=500&ucdir=newer&format=xml", sep = "")) 
    
    # XML informations download for the specific URL
    xml.contrib <- xmlToList(xmlTreeParse(url.contrib,useInternalNodes=TRUE))
    
    list.contrib <- xml.contrib$query$usercontribs
    xml.contribbis <- xml.contrib

    # Management of the uclimit argument (maximum item per page)
    while (!is.null(xml.contribbis$"query-continue")) {
      continue <- xml.contribbis$"query-continue"$usercontribs
      
      url.contrib <-  GET(paste(url.contrib, "&uccontinue=", continue, sep = ""))
      xml.contribbis <- NULL
      xml.contribbis <- xmlToList(xmlTreeParse(url.contrib,useInternalNodes=TRUE))
      list.contrib <- c( list.contrib,xml.contribbis$query$usercontribs)
    }
    out$user <- c(test$user, domain)

    # Specification of Namespace information and others outputs
    
   if(length(list.contrib) >0) # at least one contribution
   {     
     out$contribs <- matrix(nrow = length(list.contrib), ncol = length(props)+3)
  
     for (i in 1:length(list.contrib))
     {  
       for (j in 1 : length(props))
       {
         out$contribs[i,j] <- list.contrib[i]$item[name = props[j]] 
       }
       out$contribs[i,length(props)+1] <-  list.contrib[i]$item["ns"] 
       if (out$contribs[i,length(props)+1] == 0) { out$contribs[i,length(props)+2] <- "Subject" ;  out$contribs[i,length(props)+3] <- "Main Article" }      
       if (out$contribs[i,length(props)+1] == 2) { out$contribs[i,length(props)+2] <- "Subject" ; out$contribs[i,length(props)+3] <- "User"}      
       if (out$contribs[i,length(props)+1] == 4) { out$contribs[i,length(props)+2] <- "Subject" ; out$contribs[i,length(props)+3] <- "Wikipedia"}      
       if (out$contribs[i,length(props)+1] == 6) { out$contribs[i,length(props)+2] <- "Subject" ; out$contribs[i,length(props)+3] <-   "File"}      
       if (out$contribs[i,length(props)+1] == 8) { out$contribs[i,length(props)+2] <- "Subject" ; out$contribs[i,length(props)+3] <-   "MediaWiki"}      
       if (out$contribs[i,length(props)+1] == 10) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Template"}      
       if (out$contribs[i,length(props)+1] == 12) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Help"}      
       if (out$contribs[i,length(props)+1] == 14) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Category"}      
       if (out$contribs[i,length(props)+1] == 100) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Portal"}      
       if (out$contribs[i,length(props)+1] == 108) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Book"}
       if (out$contribs[i,length(props)+1] == 118) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Draft"}      
       if (out$contribs[i,length(props)+1] == 446) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Education Program"}      
       if (out$contribs[i,length(props)+1] == 710) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "TimedText"}      
       if (out$contribs[i,length(props)+1] == 828) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Module"}          
       if (out$contribs[i,length(props)+1] == 2600) { out$contribs[i,length(props)+2] <-   "Subject" ; out$contribs[i,length(props)+3] <-   "Topic"}      
       if (out$contribs[i,length(props)+1] == 1) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Talk"}      
       if (out$contribs[i,length(props)+1] == 3) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "User talk"}      
       if (out$contribs[i,length(props)+1] == 5) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Wikipedia Talk"}      
       if (out$contribs[i,length(props)+1] == 7) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "File talk"}      
       if (out$contribs[i,length(props)+1] == 9) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "MediaWiki talk"}      
       if (out$contribs[i,length(props)+1] == 11) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Template talk"}      
       if (out$contribs[i,length(props)+1] == 13) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Help talk"}      
       if (out$contribs[i,length(props)+1] == 15) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Category talk"}      
       if (out$contribs[i,length(props)+1] == 101) { out$contribs[i,length(props)+2] <-   "Talk" ;  out$contribs[i,length(props)+3] <-   "Portal talk"}      
       if (out$contribs[i,length(props)+1] == 109) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Book talk"}      
       if (out$contribs[i,length(props)+1] == 119) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Draft talk"}      
       if (out$contribs[i,length(props)+1] == 447) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Education talk"}      
       if (out$contribs[i,length(props)+1] == 711) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "TimedText talk"}      
       if (out$contribs[i,length(props)+1] == 829) { out$contribs[i,length(props)+2] <-   "Talk" ; out$contribs[i,length(props)+3] <-   "Module talk"}      
       if (out$contribs[i,length(props)+1] == -1) { out$contribs[i,length(props)+2] <-   "Virtual" ; out$contribs[i,length(props)+3] <-   "Special"}      
       if (out$contribs[i,length(props)+1] == -2) { out$contribs[i,length(props)+2] <-   "Virtual" ; out$contribs[i,length(props)+3] <-   "Media"}
       if (!(out$contribs[i,length(props)+1] %in% c(0,2,4,6,8,10,12,14,100,108,118,446,710,828,2600,1,3,5,7,9,11,13,15,101,109,119,447,711,829))) 
       { out$contribs[i,length(props)+2] <-   "Other" ; out$contribs[i,length(props)+3] <-   "Other"}
       
     } # end loop on contributions 
     
     colnames(out$contribs) <- c(props,"ns","nscat","nssubj")
     rownames(out$contribs) <- NULL
     out$contribs <- as.data.frame(out$contribs)
     out$contribs$title <- iconv(out$contribs$title,from="UTF-8",to="UTF-8")
     if (any(props =="comment")) {out$contribs$comment <- iconv(out$contribs$comment,from="UTF-8",to="UTF-8")}
   #Encoding(out$contribs$title)
   } # end if at least one contribution
   
  } # end if user.name and domain are valid
  
  out$testWikiUser <- list(test$takeOnlyFirst,test$test,test$warnMessage) # do not keep twice user info, if the user exists
  if(!is.null(test$warnMessage)){for( i in 1:length(test$warnMessage) ) {warning(test$warnMessage[i]) } }
  class(out) <- c("userContribsClass")
  
  return(out)
  
}
