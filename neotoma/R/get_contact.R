
#' @title Get contact information.
#'
#' @description A function to obtain contact information for data contributors from the Neotoma Paleoecological Database.
#' @importFrom jsonlite fromJSON
#' @importFrom httr content GET
#' @param contactid Contact ID is a numerical value associated with the Neotoma
#'    Contact table's numerical Contact ID.
#' @param contactname A character string indicating the data contributors' project,
#'    organization or personal name.  May be a partial string and can include wildcards.
#' @param contactstatus The current status of the contact.  Possible values include:
#'    active, deceased, defunct, extant, inactive, retired, unknown.
#' @param familyname A character string.  Full or partial string indicating the
#'    contact's last name.
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return The function takes parameters defined by the user and returns a list
#'    of contact information supplied by the Neotoma Paleoecological Database.
#'    The user may define all or none of the possible fields.  The function contains
#'    data chacks for each defined parameter.
#'
#'    The function returns either a single item of class \code{"try-error"} describing
#'    the reason for failure (either mis-defined parameters or an error from the Neotoma API),
#'    or a table of contacts, with rows corresponding to the number of individual
#'    contacts returned by the Neotoma API.  Each row entry includes the following parameters:
#'
#'  \item{ \code{contact.name} }{Full name of the person, last name first (e.g. \code{"Simpson, George Gaylord"}) or name of organization or project (e.g. \code{"Great Plains Flora Association"}).}
#'  \item{ \code{contact.status} }{Current status of the person, organization, or project. Field links to the ContactStatuses lookup table.}
#'  \item{ \code{family.name} }{Family or surname name of a person.}
#'  \item{ \code{leading.initials} }{Leading initials for given or forenames without spaces (e.g. \code{"G.G."}).}
#'  \item{ \code{given.names} }{Given or forenames of a person (e.g. \code{"George Gaylord"}). Initials with spaces are used if full given names are not known (e.g. \code{"G. G")}.}
#'  \item{ \code{suffix} }{Suffix of a person's name (e.g. \code{"Jr."}, \code{"III"}).}
#'  \item{ \code{title} }{A person's title (e.g. \code{"Dr."}, \code{"Prof."}, \code{"Prof. Dr"}).}
#'  \item{ \code{phone} }{Telephone number.}
#'  \item{ \code{fax} }{Fax number.}
#'  \item{ \code{email} }{Email address.}
#'  \item{ \code{url} }{Universal Resource Locator, an Internet World Wide Web address.}
#'  \item{ \code{address} }{Full mailing address.}
#'  \item{ \code{notes} }{Free form notes or comments about the person, organization, or project.}
#'  \item{ \code{contact.id} }{Unique database record identifier for the contact.}
#'  \item{ \code{alias.id} }{The ContactID of a person's current name. If the AliasID is different from the ContactID, the ContactID refers to the person's former name.}
#'
#' @examples \dontrun{
#' #  To find all data contributors who are active:
#' active.cont <- get_contact(contactstatus = 'active')
#'
#' # To find all data contributors who have the last name "Smith"
#' smith.cont <- get_contact(familyname = 'Smith')
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords IO connection
#' @export
get_contact <- function(contactid, contactname, contactstatus, familyname){

  base.uri <- 'http://api.neotomadb.org/v1/data/contacts'

  cl <- as.list(match.call())
  cl[[1]] <- NULL
  cl <- lapply(cl, eval, envir = parent.frame())

  #  Pass the parameters to param_check to make sure everything is kosher.
  error_test <- param_check(cl)
  if(error_test[[2]]$flag == 1){
    stop(paste0(unlist(error_test[[2]]$message), collapse='\n  '))
  } else{
    cl <- error_test[[1]]
  }

  neotoma_content <- content(GET(base.uri, query = cl), as = "text")
  
  if (identical(neotoma_content, "")) stop("")
  
  aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)
  
  if (aa[[1]] == 0){
    stop(paste('Server returned an error message:\n', aa[[2]]), call. = FALSE)
  }
  if (aa[[1]] == 1){
    aa <- aa[[2]]
    message(paste0('The API call was successful, you have returned ',
                   length(aa), ' records.\n'))
  }

  if (inherits(aa, "try-error")) {
      output <- neotoma_content
  } else {
    
    # replace NULL values:
    aa <- lapply(aa, function(x) ifelse(x == "NULL", NA, x))
    
    output <- data.frame(contact.name = sapply(aa, '[[', "ContactName"),
                         contact.status = sapply(aa, '[[', "ContactStatus"), 
                         family.name = sapply(aa, '[[', "FamilyName"),
                         leading.initials = sapply(aa, '[[', "LeadingInitials"),
                         given.names = sapply(aa, '[[', "GivenNames"),
                         suffix = sapply(aa, '[[', "Suffix"),
                         title = sapply(aa, '[[', "Title"),
                         phone = sapply(aa, '[[', "Phone"),
                         fax = sapply(aa, '[[', "Fax"),
                         email = sapply(aa, '[[', "Email"),
                         url = sapply(aa, '[[', "URL"),
                         address = sapply(aa, '[[', "Address"),
                         notes = sapply(aa, '[[', "Notes"),
                         contact.id = sapply(aa, '[[', "ContactID"),
                         stringsAsFactors = FALSE)
  }
  
  output
}
