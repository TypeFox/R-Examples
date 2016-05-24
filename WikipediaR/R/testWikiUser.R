#' internal function testWikiUser
#' 
#' @description 
#' internal function to test validity of a wikidepia domain and user, 
#' called in \code{userInfo} and \code{userContribs} functions.
#' 
#' @param user.name a character value providing the name of the specific wikipedia user
#' @param domain a character value providing the language of the wikipedia page.The default value is "en" for "english language".
#' 
#' @return a list of four elements, the first is \code{takeOnlyFirst}, a boolean indicating if the class of \code{user.name} parameter is invalid, 
#' for example vector, list, matrix..., and in that case, the only the first element is considered.
#'  The second element is \code{test}, an integer with value: 
#' \itemize{
#' \item{4}{ for invalid domain,}
#' \item{3}{ for an empty parameter user.name,}
#' \item{2}{ when Wikipedia does not have an user with this exact name,}
#' \item{0}{ for valid existing user.name. }
#' }
#' The third element is \code{user}, a vector containing the user name, the user identifier, and domain. If the user does not exist, this item is not created.
#' The last element, \code{warnMessage}, is a vector of warning messages.

#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#' 
#' @importFrom XML xmlToList xmlTreeParse htmlParse 
#' @importFrom httr GET 
#---------------------------------------------------------------

testWikiUser <- function(user.name =NULL, domain = "en") 
{
  # initialize results
  out <- NULL
  test <- 0  
  takeOnlyFirst <- FALSE
  warnMessage <- NULL
  #-------------------------------#  
  # verify valid type of user.name : vector, list, matrix... are invalid
  # do not stop but take the first element and warn
  
  if(length(user.name)> 1 | is.list(user.name))  
  {
    takeOnlyFirst <- TRUE
    warnMessage <- c(warnMessage,"Invalid dimension for the user.name, only the first element is considered.") 
    user.name <- unlist(user.name) [1]
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
      warnMessage <- c(warnMessage,paste("Invalid domain: Impossible to connect to", paste("http://",domain,".wikipedia.org/",sep=""))) 
    }
    else{
      
      #-------------------------------#  
      ## empty parameter user.name
      miss.user.name <- FALSE
      if(is.null(user.name)) {miss.user.name <- TRUE}
      else{        
        user.name <- as.character(user.name)
        if (user.name =="") { miss.user.name <- TRUE }
      }
      if(miss.user.name)    
      {test <- 3 ; warnMessage <- c(warnMessage,"The parameter 'user.name' is empty")}
      else{
        # manage encoding and spaces
        Encoding(user.name) <- "latin1"
        # user.name2 <- gsub(" ",replacement ="_",x = user.name)
        url.user <- GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&&format=xml&list=users&ususers=",user.name, sep=""))
        xml.user <- xmlToList(xmlTreeParse(url.user, useInternalNodes = TRUE))
        if (any(names(xml.user$query$users$user)=="missing")) 
        { 
          test <- 2
          user <- NULL
          warnMessage <- c(warnMessage ,"Wikipedia does not have an user with this exact name.")
        }
        else
        {
          user <- vector (length=2)
          user[1] <- xml.user$query$users$user[2]
          user[2] <- xml.user$query$users$user[1]
          names(user) <- c("user", "userid")
          out$user <- user
        }
      } # end parameter user not empty
    } # end domain exists
   } # end valid domain
  
  out$takeOnlyFirst <-takeOnlyFirst
  out$test <- test  
  out$warnMessage <- warnMessage
  return (out)
}