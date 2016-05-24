#' Guess Italian Fiscal Code
#' 
#' The function tries to guess regular fiscal code, extracting relevant
#' alphanumeric digits from surname, name, birth date, gender and 'codice
#' catastale' (computing the last character, the control digit).
#' 
#' @param surname character, surname
#' @param name character, names
#' @param birth_date Date, date of birth
#' @param female logical, female indicator variable (\code{FALSE} = man,
#' \code{TRUE} = woman)
#' @param codice_catastale  Italian 'codice catastale' (an identifier) of
#' the 'comune' of birth
#' @return The function return a character vector of fiscal code. 
#' @examples
#' 
#' ## using fictious data
#' Surnames <- c("Rossi", "Bianchi")
#' Names <- c("Mario", "Giovanna")
#' Birthdates <- as.Date(c("1960-01-01", "1970-01-01"))
#' Female <- c(FALSE, TRUE)
#' Comune_of_birth <- c("F205", # milan
#'                      "H501") # rome
#' guess_fc(Surnames, Names, Birthdates, Female, Comune_of_birth)
#' 
#' @export 
guess_fc <- function(surname = NULL,
                     name = NULL,
                     birth_date = NULL,
                     female = NULL,
                     codice_catastale = NULL) 
{

  ## validate input
  if(! is.character(surname))
    stop("surname must be a character vector.")
  if(! is.character(name))
    stop("name must be a character vector.")
  if(! inherits(birth_date, "Date"))
    stop("birth_date must be a Date vector.")
  if(! is.logical(female))
    stop("female must be a logical vector.")
  if(! is.character(codice_catastale))
    stop("codice_catastale must be a character vector.")

  ## validate lengths
  Len <- unlist(lapply(list(surname, name, birth_date,
                            female, codice_catastale),
                length))
  if (stats::var(Len) > 0){
    stop("all vector must be of the same length")
  }
  
  ## normalize input
  surname <- keepAlpha(surname)
  surname[surname %in% ""] <- NA
  name <- keepAlpha(name)
  name[name %in% ""] <- NA
  codice_catastale[codice_catastale %in% ""]  <- NA  
  female <- as.integer(female)
  year <- as.integer(format(birth_date, "%Y"))
  month <- as.integer(format(birth_date, "%m"))
  day <- as.integer(format(birth_date, "%d"))

  ## locate rows with NAs
  NAS <- apply(do.call("rbind",
                       lapply(list(surname, name, birth_date,
                                   female, codice_catastale),
                              is.na)),
               2, any)

  ## return value
  rval <- rep(NA_character_, Len[1])
  rval[!NAS] <- .Call("reg_guess_fc",
                      surname[!NAS],
                      name[!NAS],
                      year[!NAS],
                      month[!NAS],
                      day[!NAS],
                      female[!NAS],
                      codice_catastale[!NAS],
                      PACKAGE = "ifctools")
  rval

}

keepAlpha <- function(x)
  gsub("[^[:alpha:]]", "", toupper(x), perl =TRUE)
