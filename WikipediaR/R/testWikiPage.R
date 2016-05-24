 #' internal function testWikiPage
#' @description internal function to test validity of a wikidepia page, called in \code{links}, \code{backLinks}, and \code{contribs} functions.
#' @param page numeric identifier or character title of the specific wikipedia page
#' @param domain a character value providing the language of the wikipedia page.The default value is "en" for "english language".
#' 
#' @return a list of four elements, the first is \code{takeOnlyFirst}, a boolean indicating if the class of \code{page} parameter is invalid, 
#' for example vector, list, matrix..., and in that case, the only the first element is considered.
#'  The second element is \code{redirPage}, the title page of the redirected page.
#'  The third element is \code{test}, an integer with value: 
#' \itemize{
#' \item{4}{ for invalid domain,}
#' \item{3}{ for an empty parameter page,}
#' \item{2}{ when Wikipedia does not have an article with this exact name,}
#' \item{1}{ for ambiguous page, direct or redirect,} 
#' \item{0}{ for valid an unambiguous page, direct or redirect. }
#' }
#' The last element, \code{warnMessage}, is a vector of warning messages.
#' 
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#' 
#' @importFrom XML xmlToList xmlTreeParse htmlParse 
#' @importFrom httr GET
#' @importFrom utils URLencode
#' 
testWikiPage <- function(page =NULL, domain = "en") 
{
  # initialize results
  test <- 0  
  takeOnlyFirst <- FALSE
  warnMessage <- NULL
  redirpage <- NULL
  #-------------------------------#  
  # verify valid type of page : vector, list, matrix... are invalid
  # do not stop but take the first element and warn
  
  if(length(page)> 1 | is.list(page))  
  {
    takeOnlyFirst <- TRUE
    warnMessage <- c(warnMessage,"Invalid dimension for the page, only the first element is considered.") 
    page <- unlist(page) [1]
  }
  #-------------------------------#
  ## verify that domain is valid
  
  if(is.character(domain)==FALSE | length(domain)>1)         
  {
    test <- 4 
    warnMessage <- ifelse(length(domain)==1, 
                      c(warnMessage,paste("Invalid domain:",domain,"is not a character string.")),
                      c(warnMessage,"Invalid domain: dimensions are not valid.")
                     )
  }
  else{
    # try to connect to the welcome page of the specified domain
    tryconnect <- tryCatch(htmlParse(GET(paste("http://",domain,".wikipedia.org/", sep=""))), error=function(e) e)
    if (any(class(tryconnect) == "XML_IO_LOAD_ERROR")  )              
    {
      test<-4 
      warnMessage <- c(warnMessage,paste("Invalid domain:", paste("Impossible to connect to http://",domain,".wikipedia.org/",sep=""))) 
    }
    else{

      #-------------------------------#  
      ## empty parameter page
      miss.page <- FALSE
      if(is.null(page)){miss.page <- TRUE}
      else{ if (as.character(page) =="") { miss.page <- TRUE }}
      if(miss.page)    
      {test <- 3 ; warnMessage <- c(warnMessage,"The parameter 'page' is empty")}
      else{
        
        if (is.numeric(page))  {url.part1 <- paste("http://",domain,".wikipedia.org/w/api.php?action=query&pageids=",page, sep="")}
        if (is.character(page))
        {
          # manage encoding and spaces     
          page1 <- gsub(" ",replacement ="_",x = page)
          page2 <- URLencode(iconv(page1,to="UTF-8"))
          url.part1 <- paste("http://",domain,".wikipedia.org/w/api.php?action=query&titles=",page2, sep="")
        }
        
        xml.page <- xmlToList(xmlTreeParse(GET(paste(url.part1, "&prop=info&format=xml", sep= "")),useInternalNodes = TRUE))
        
        #-------------------------------#   
        ## verify that the page exists
        if(!any(names(xml.page$query$pages$page) =="pageid") | any(names(xml.page$query$pages$page) =="missing"))
        {test <- 2 ; warnMessage <- c(warnMessage,"Wikipedia does not have an article with this exact name") }
        else
        {
          #-------------------------------#  
          #for redirected page
          
          if(any(names(xml.page$query$pages$page) =="redirect"))
          {
            redirpage <- xmlToList(xmlTreeParse(GET(paste(url.part1,"&prop=links&format=xml",sep="")), useInternalNodes = TRUE))$query$pages$page$links$pl [2]
            warnMessage <- c(warnMessage,paste(page,"is redirected to", redirpage))
            page2 <- gsub(" ",replacement ="_",x = redirpage)
            url.part1 <- paste("http://",domain,".wikipedia.org/w/api.php?action=query&titles=",page2, sep="")
          }
          
          #-------------------------------#     
          ## page is ambiguous - or homonymy
          url.cat <- GET(paste(url.part1,"&prop=categories&continue&cllimit=500&format=xml",sep=""))
          tree.page <- xmlTreeParse(url.cat,useInternalNodes = TRUE )
          xml.page.2 <- xmlToList(tree.page) 
   
          if(any(names(unlist(xml.page.2$query$pages$page))=="categories"))
          {
            if (any(length(grep("DISAMBIG",toupper(unlist(xml.page.2$query$pages$page$categories) ))) )            
               |  any(length(grep("HOMONYMIE",toupper(unlist(xml.page.2$query$pages$page$categories) ))) )    
              # special french case : homonymie for ambibuous
            )
          
          {  test <- 1 ; warnMessage <- c(warnMessage,"This disambiguation page lists articles associated with the same title.") } 
          } # categories in the result of url.cat
        } # end existing page
        
      } # end not missing page
      
    } # end valid domain
    
    
  } # end domain not missing
 
  return (list(takeOnlyFirst=takeOnlyFirst , redirPage = redirpage, test = test, warnMessage = warnMessage))
}
