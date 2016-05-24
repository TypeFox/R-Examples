#' @title
#' Parse organizational identity numbers
#'
#' @description
#' Check and convert a vector of organizational identity numbers. 
#' 
#' @details 
#' The following format is accepted: 
#' \itemize{
#'   \item character: \code{GNNNNN-NNNC}
#' }
#' 
#' @param oin Vector with swedish organizational identity numbers in character format. See details.
#' 
#' @references 
#' \href{http://www.skatteverket.se/download/18.70ac421612e2a997f85800040284/1359707510840/70909.pdf}{SKV 709}
#' 
#' @return
#' Character vector (of class \code{oin} and \code{AsIs}) with swedish organizational identity numbers.
#'
#' @examples
#' ex_oin <- c("556000-4615", "232100-0156", "802002-4280", "8020024280", "AA2002-4280")
#' as.oin(ex_oin)
#' 
#' @export
as.oin <- function(oin){
  UseMethod("as.oin")
}

#' @export
as.oin.character <- function(oin){
  suppressWarnings(
    correct <- 
      is.character(oin) &
        grepl(pattern = "^[0-9]{6}-[0-9]{4}$", oin) &
        as.numeric(substr(oin,3,3)) >= 2 
  )
  newoin <- oin
  newoin[!correct] <- NA
  
  # Warning for incorrect pin
  isna <- is.na(newoin)
  if(any(isna)) {
    warning("Erroneous oin(s) (set to NA).")
  }
  
  # Add class
  class(newoin) <- c("AsIs", "oin", "character")
  
  return(newoin)
}

#' @export
as.oin.oin <- function(oin){
  oin
}

#' @export
as.oin.factor <- function(oin){
  as.oin(as.character(oin))
}

#' @export
as.oin.default <- function(oin){
  stop("Object of class ", paste(class(oin), collapse = ", "), 
       " can not be coerced to oin!"
  )
}



#' @title
#' Test if a character vector contains correct \code{oin}
#' 
#' @description
#' Test which elements in a text vector that contains organization identity 
#' number.
#' 
#' @param oin Character vector to be tested if it is an \code{oin} of the right format.
#' 
#' @return
#' Logical vector indicating if the elements can be an organization identity number.
#'
#' @examples
#' ex_oin <- c("556000-4615", "232100-0156", "802002-4280", "8020024280", "AA2002-4280")
#' is.oin(ex_oin)
#' 
#' @export
is.oin <- function(oin){
  inherits(oin, "oin")
}

#' @title
#' Check the control numbers for \code{oin}
#' 
#' @description
#' Calculates the control number using the Luhn algorithm and compare it with the control number in the organization identity number (oin).
#' 
#' @param oin A vector of class \code{oin}. See \link{as.oin}.
#' @param force_logical If TRUE, force all NA in oin to be FALSE. Default is FALSE.
#' 
#' @references 
#' \href{http://www.skatteverket.se/download/18.70ac421612e2a997f85800040284/1302507382017/70909.pdf}{SKV 709}
#' 
#' @return
#' Logical vector indicating if a oin is correct (\code{TRUE}) or not (\code{FALSE})
#'
#' @examples
#' ex_oin <- c("556000-4615", "232100-0156", "802002-4280", "232100-0157", "802002-4281")
#' oin_ctrl(ex_oin)
#' 
#' @export
oin_ctrl <- function(oin, force_logical = FALSE){
  if(force_logical){
    if(!is.oin(oin)) oin <- suppressWarnings(as.oin(oin))
  } else {
    if(!is.oin(oin)) oin <- as.oin(oin)
  }
  
  oin_char <- paste0(substr(oin,1,6), substr(oin,8,11))
  oin_char[is.na(oin)] <- NA
    
  res <- vapply(oin_char, luhn_algo, integer(1), USE.NAMES = FALSE, 
                multiplier = c(2, 1, 2, 1, 2, 1, 2, 1, 2, 0))
  res <- as.integer(substr(oin, 11, 11)) == res
  if(force_logical) res[is.na(res)] <- FALSE
  
  res
}

#' @title
#' Calculate organization group from \code{oin}
#' 
#' @description
#' Calculates the organization group from the organization number.
#' 
#' @inheritParams oin_ctrl
#' 
#' @references 
#' \href{http://www.skatteverket.se/download/18.70ac421612e2a997f85800040284/1302507382017/70909.pdf}{SKV 709}
#' 
#' @return
#' Factor with organization categories.
#' 
#' @examples
#' ex_oin <- c("556000-4615", "232100-0156", "802002-4280")
#' oin_group(ex_oin)
#'
#' @export
oin_group <- function(oin){
  if(!is.oin(oin)) oin <- as.oin(oin)
  
  as.factor(vapply(X = oin, 
                   FUN = oin_group_element, 
                   FUN.VALUE = character(1), 
                   USE.NAMES = FALSE))
}
