#############################################################################
#' provides general information about a Wikipedia user
#' 
#' @description provides general information about a specific Wikipedia user
#' 
#' @param user.name a character value for the name of a user
#' @param domain a character value providing the language of the wikipedia page.The default value is "en" for "english language".
#' @param  usprop What pieces of information to include (separate with '|'):
#' \itemize{
#' \item{\code{groups}}{: lists all the groups the user(s) belongs to}
#' \item{\code{implicitgroups}}{: lists all the groups a user is automatically a member of}
#' \item{\code{rights}}{: lists all the rights the user(s) has}
#' \item{\code{editcount}}{: adds the user's edit count}
#' \item{\code{registration}}{: adds the user's registration timestamp}
#' \item{\code{emailable}}{: tags if the user can and wants to receive email through [[Special:Emailuser]]}
#' \item{\code{gender}}{: tags the gender of the user. Returns "male", "female", or "unknown"}
#' }
#' Default: groups|implicitgroups|rights|editcount|registration|emailable|gender

#' @return an object of class \code{userInfoClass} containing:
#' \itemize{
#' \item{\code{call}}{ the command line}
#' \item{\code{rights}}{ only if asked: a vector of rights. If not asked, this item does not exist.}
#' \item{\code{groups}}{ only if asked: a vector of groups names. If not asked, this item does not exist.}
#' \item{\code{implicitgroups}}{ only if asked: a vector of implicit groups names. If not asked, this item does not exist.}
#' \item{\code{info}}{ a data frame containing the other asked properties of user (at least userid and name). }
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
#' 
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#' 
#' @details
#' This function uses the MediaWiki API query syntax: "list=users". 
#' For more details, see \url{https://www.mediawiki.org/wiki/API:Users}
#'
#' @examples 
#' \dontrun{
#' LouiseInfo <- userContribs(user.name = "Louise", domain = "en")
#' LouiseInfo
#' bobInfo <- userInfo(user.name = "bob", domain = "en")
#' bobInfo
#' ## try a user that does not exist (at the moment of the redaction of this help page)
#' userInfo(user.name="Louise Baschet", domain ="fr")
#' }
#' @importFrom XML xmlToList xmlTreeParse htmlParse 
#' @importFrom httr GET
#' 
#' @export 
#' 
#############################################################################

userInfo <- function (user.name=NULL, domain="en", usprop = "groups|implicitgroups|rights|editcount|registration|emailable|gender") 
  {
    
  # initialize output 
  out <- NULL
  out$call <- match.call()
  out$rights <- NULL
  out$groups <- NULL
  out$implicitgroups <- NULL
  
  ## verify usprop parameters
  if (is.null(usprop)) { stop("usprop argument is required.")}
  if(!is.character(usprop)){stop("usprop argument must be a character string")}
  
  props <- unlist(strsplit(usprop, split ="|", fixed = TRUE))
  
  if (! all(props %in% c("groups","implicitgroups","rights","editcount","registration","emailable","gender")))
  {stop("usprop argument is not valid")}
  
  ## verify the domain and user.name validity
  test <- testWikiUser(user.name = user.name, domain = domain)
  
  if (!test$test %in% c(4,3,2)) {
    
    if (test$takeOnlyFirst) { user.name <- unlist(user.name )[1]}
    
    # if anonymous 
    print(test)
    if (test$user[1] =="") {
      test$warnMessage <- c(test$warnMessage,"This user is anonymous.")
    }
    else{
      # URL building
      url.user <- GET(paste("http://",domain,".wikipedia.org/w/api.php?action=query&list=users&ususers=", as.character(user.name), 
                                "&continue&format=xml&usprop=",usprop, sep = ""))
      print(url.user)
      # XML informations download for the specific URL
      xml.user <- xmlToList(xmlTreeParse(url.user ,useInternalNodes=TRUE) )
      print(unlist(xml.user$query$users$user))
      if (any(props == "rights")) {out$rights <- as.vector(unlist(xml.user$query$users$user$rights))}
      if (any(props == "groups")) {out$groups <- as.vector(unlist(xml.user$query$users$user$groups))}
      if (any(props == "implicitgroups")) {out$implicitgroups <- as.vector(unlist(xml.user$query$users$user$implicitgroups))}
      
      if (any(props %in% c("rights","groups","implicitgroups")))
        # others informations are in $.attrs
        {
        out$info <- as.data.frame(t(xml.user$query$users$user$.attr) )
        rownames(out$info) <- NULL
        if(any(props =="emailable")) {
          if (any(names(xml.user$query$users$user$.attr)=="emailable"))
          { out$info["emailable"] <- TRUE  }
          else{ out$info["emailable"] <- FALSE 
          colnames(out$info)[ncol(out$info)] <-"emailable"
          }
        } # end if emailable is asked
      } 
      else {
        # informations are in $user directly
        out$info <- as.data.frame(t(xml.user$query$users$user) )
        rownames(out$info) <- NULL
        if(any(props =="emailable")) {
          if (any(names(xml.user$query$users$user)=="emailable"))
          { out$info["emailable"] <- TRUE  }
          else{ out$info["emailable"] <- FALSE 
                colnames(out$info)[ncol(out$info)] <-"emailable"
          }
        } # end if emailable is asked
      } 
      
    # reorder info
    out$info <- as.data.frame(c(out$info[c(2,1)], domain, out$info[3:length(out$info)]))
    names(out$info)[3] <- "domain"
    rownames(out$info) <- NULL
    } # end if user.name and domain are valid
  
  }
  out$testWikiUser <- list(test$takeOnlyFirst,test$test,test$warnMessage) # do not keep twice user info
  
  if(!is.null(test$warnMessage)){for( i in 1:length(test$warnMessage) ) {warning(test$warnMessage[i]) } }
 
  
  class(out) <- c("userInfoClass")
  
  out
  
  
}
