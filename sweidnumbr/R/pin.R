#' @title
#' Parse personal identity numbers to ABS format
#' 
#' @description
#' \code{as.pin} Converts personal identity numbers of different formats to standard (ABS) 
#' pin format \code{YYYYMMDDNNNC} where \code{YYYYMMDD} is the date of birth, \code{NNN} 
#' is the birth number and \code{C} is the
#' control number.
#' \code{is.pin} checks wether an R object is of class "pin".
#' 
#' @details
#' \code{as.pin} converts different formats of swedish personal identity numbers to
#' the standard ABS format. The formats that can be converted are:
#' \itemize{
#'   \item numeric: \code{YYYYMMDDNNNC}
#'   \item numeric: \code{YYMMDDNNNC} (assuming < 100 years of age)
#'   \item character: \code{"YYYYMMDDNNNC"}
#'   \item character: \code{"YYMMDD-NNNC"},  \code{"YYMMDD+NNNC"}
#'   \item character: \code{"YYYYMMDD-NNNC"}
#'   \item character: \code{"YYMMDDNNNC"} (assuming < 100 years of age)
#' }
#' (where "C" can be substituted by characters "A", "T" or "X" if "YYYY" < 1967).
#' 
#' @param pin Vector with swedish personal identity numbers in character or numeric format. 
#' See details.
#' 
#' @references 
#' \itemize{
#'  \item \href{https://www.skatteverket.se/download/18.8dcbbe4142d38302d74be9/1387372677724/717B06.pdf}{Population registration in Sweden}
#'  \item \href{https://www.skatteverket.se/download/18.1e6d5f87115319ffba380001857/1285595720207/70408.pdf}{SKV 704}
#'  \item \href{http://www.riksdagen.se/sv/Dokument-Lagar/Utredningar/Statens-offentliga-utredningar/Personnummer-och-samordningsnu_GWB360/}{SOU 2008:60 : Personnummer och samordningsnummer}
#'  \item \emph{Personnummer: information fran Centrala folkbokförings- och uppbördsnämnden.} (1967). Stockholm
#'  \item \emph{Den svenska folkbokföringens historia under tre sekel.} (1982). Solna: Riksskatteverket \href{http://www.skatteverket.se/privat/folkbokforing/omfolkbokforing/folkbokforingigaridag/densvenskafolkbokforingenshistoriaundertresekler.4.18e1b10334ebe8bc80004141.html}{URL}
#' }
#' @return
#' \code{as.pin} returns a vector of class "pin" (with additional classes "AsIs" and character) 
#' with swedish personal identity numbers with standard ABS format \code{"YYYYMMDDNNNC"}.
#' \code{is.pin} returns \code{TRUE} if \code{pin} is of class "pin", otherwise \code{FALSE}.
#'
#' @examples
#' # Examples taken from SKV 704 (see references)
#' ex_pin1 <- c("196408233234", "640823-3234", "19640823-3234")
#' as.pin(pin = ex_pin1)
#' ex_pin2 <- c("6408233234")
#' as.pin(ex_pin2)
#' ex_pin3 <- c(6408233234, 196408233234)
#' as.pin(ex_pin3)
#' ex_pin4 <-rep(c("20121209-0122", "201212090122", "121209-0122", "1212090122"),250)
#' as.pin(ex_pin4)
#' ex_pin5 <-c("205012090122", "186512090122", "121209-0122", "121209-012A")
#' as.pin(pin = ex_pin5)
#' pin <-c("201212090122", "201212090122", "121209-0122", "1212090122")
#' 
#' @export
#' @name as.pin
as.pin <- function(pin){
  UseMethod("as.pin")
}

#' @export
as.pin.numeric <- function(pin){
  pin <- as.character(pin)
  pin[!is.na(pin)] <- stringr::str_pad(pin[!is.na(pin)], 10, pad = "0")
  as.pin(pin)
}

#' @export
as.pin.pin <- function(pin){
  pin
}

#' @export
as.pin.factor <- function(pin){
  as.pin(as.character(pin))
}

#' @export
as.pin.default <- function(pin){
  stop("Object of class ", paste(class(pin), collapse = ", "), 
       " can not be coerced to pin!"
  )
}

# Vector of only NA:s can also get the class attribute pin
#' @export
as.pin.logical <- function(pin){
  if (all(is.na(pin))){
    structure(pin, class = c("AsIs", "pin", "character"))
  } else{
    NextMethod()
  }
}
  

#' @export
as.pin.character <- function(pin){
  all_pins <- pin
  pin <- all_pins[!is.na(all_pins)]
  
  formats <- character(8)
  # format 1: "YYYYMMDDNNNC"
  formats[1] <- "^(18|19|20)[0-9]{2}(0[1-9]|1[0-2])([06][1-9]|[1278][0-9]|[39][0-1])[0-9]{4}$"
  # format 2: "YYYYMMDD-NNNC"
  formats[2] <- "^(18|19|20)[0-9]{2}(0[1-9]|1[0-2])([06][1-9]|[1278][0-9]|[39][0-1])[-][0-9]{4}$"
  # format 3: "YYMMDD-NNNC"
  formats[3] <- "^[0-9]{2}(0[1-9]|1[0-2])([06][1-9]|[1278][0-9]|[39][0-1])[-+][0-9]{4}$"
  # format 4: "YYMMDDNNNC"
  formats[4] <- "^[0-9]{2}(0[1-9]|1[0-2])([06][1-9]|[1278][0-9]|[39][0-1])[0-9]{4}$"
  
  #  Additional formats for old "pins" for people deceased 1947 - 1967 (i.e. ctrl numbr is missing/replaced with A,T or X)
  # format 5: "YYYYMMDDNNNC"
  formats[5] <- "^(18[0-9]{2}|19([0-5][0-9]|6[0-6]))(0[1-9]|1[0-2])([06][1-9]|[1278][0-9]|[39][0-1])[0-9]{3}[ATX ]$"
  # format 6: "YYYYMMDD-NNNC"
  formats[6] <- "^(18[0-9]{2}|19([0-5][0-9]|6[0-6]))(0[1-9]|1[0-2])([06][1-9]|[1278][0-9]|[39][0-1])[-][0-9]{3}[ATX ]$"
  # format 7: "YYMMDD-NNNC"
  formats[7] <- "^([0-5][0-9]|6[0-6])(0[1-9]|1[0-2])([06][1-9]|[1278][0-9]|[39][0-1])[-+][0-9]{3}[ATX ]$"
  # format 8: "YYMMDDNNNC"
  formats[8] <- "^([0-5][0-9]|6[0-6])(0[1-9]|1[0-2])([06][1-9]|[1278][0-9]|[39][0-1])[0-9]{3}[ATX ]$"
  
  # Convert
  newpin <- rep(as.character(NA), length(pin))
  
  logi_format <- logical(length(pin))
  msg <- NA
  for (i in seq_along(formats)){
    logi_format <- grepl(formats[i], x = pin)
    newpin[logi_format] <- pin_convert(pin[logi_format], format = i - (i %/% 5) * 4)
    if (any(logi_format)) {
      if (i %in% c(4, 8)) {
        msg[1] <- "pin of format YYMMDDNNNC is assumed to be less than 100 years old"
      } 
      if (i %in% 5:8) {
        msg[2] <- paste("people with birth year before 1967 and",
                        "character 'A', 'T' or 'X' instead of control number",
                        "assumed deceast before 1967.")
      }
    }
  }
  # Maximum one of each message is enough, messages are therefore stored and possibly 
  # overwritten but not printed inside the loop
  if (!isTRUE(is.na(msg))) {
    msg <- paste(stats::na.omit(msg), collapse = " and ")
    message(paste("Assumption:", paste(toupper(substring(msg, 1, 1)), substring(msg, 2), sep = "", collapse = " ")))
    }

  # Check dates
  date <- as.Date(pin_coordn_correct(structure(newpin, class = "pin")),"%Y%m%d")
  suppressWarnings( 
    correct_date <-
      !is.na(date) &
      date <= Sys.Date() &
      date >= as.Date("1830-01-01")
  )
  newpin[!correct_date] <- NA
  
  # Warning for incorrect pin
  isna <- is.na(newpin)
  if(any(isna)) {
    warning("Erroneous pin(s) (set to NA).")
  }

  all_pins[!is.na(all_pins)] <- newpin    
  class(all_pins) <- c("AsIs", "pin", "character") 
  all_pins
}

#' @rdname as.pin
#' @export
is.pin <- function(pin) inherits(pin, "pin")

#' @title
#' Check control number from \code{pin}
#' 
#' @description
#' Calculates the control number using the Luhn algorithm and compare it with the 
#' control number in the personal identity number.
#' 
#' @param pin A vector of class \code{pin}. See \link{as.pin}.
#' @param force_logical If TRUE, force all NA in pin to be FALSE. Default is FALSE.
#' 
#' @references 
#' \href{https://www.skatteverket.se/download/18.8dcbbe4142d38302d74be9/1387372677724/717B06.pdf}{Population registration in Sweden}
#' \href{https://www.skatteverket.se/download/18.1e6d5f87115319ffba380001857/1285595720207/70408.pdf}{SKV 704}
#' \href{http://www.riksdagen.se/sv/Dokument-Lagar/Utredningar/Statens-offentliga-utredningar/Personnummer-och-samordningsnu_GWB360/}{SOU 2008:60 : Personnummer och samordningsnummer}
#' 
#' @return
#' Logical vector indicating if a pin is correct (\code{TRUE}) or not (\code{FALSE})
#' 
#' @examples
#' # Examples taken from SKV 704 (see references)
#' ex_pin <- c("196408233234", "196408233235")
#' pin_ctrl(ex_pin)
#' 
#' @export
pin_ctrl <- function(pin, force_logical = FALSE){
  if(force_logical){
    if(!is.pin(pin)) pin <- suppressWarnings(as.pin(pin))
  } else {
    if(!is.pin(pin)) pin <- as.pin(pin)
  }

  res <- vapply(pin, luhn_algo, integer(1), USE.NAMES = FALSE, 
                multiplier = c(0, 0, 2, 1, 2, 1, 2, 1, 2, 1, 2, 0))
  old_pin_format <- format(pin_to_date(pin), format = "%Y") <= "1967" & grepl("*[ATX]$", pin)
  res <- as.integer(substr(pin, 12, 12)) == res | old_pin_format
  if(force_logical) res[is.na(res)] <- FALSE
  res
}

#' @title
#' Calculate sex from \code{pin}
#' 
#' @description
#' Calculates the sex from the personal identification number.
#' 
#' @inheritParams pin_ctrl
#' 
#' @references 
#' \href{https://www.skatteverket.se/download/18.8dcbbe4142d38302d74be9/1387372677724/717B06.pdf}{Population registration in Sweden}
#' \href{https://www.skatteverket.se/download/18.1e6d5f87115319ffba380001857/1285595720207/70408.pdf}{SKV 704}
#' \href{http://www.riksdagen.se/sv/Dokument-Lagar/Utredningar/Statens-offentliga-utredningar/Personnummer-och-samordningsnu_GWB360/}{SOU 2008:60 : Personnummer och samordningsnummer}
#' 
#' @return
#' Factor with label 'Male' and 'Female'.
#' 
#' @examples
#' # Examples taken from SKV 704 (see references)
#' ex_pin <- c("196408233234", "186408233224")
#' pin_sex(ex_pin)
#'
#' @export
pin_sex <- function(pin){
  if(!is.pin(pin)) pin <- as.pin(pin)
  female <- as.numeric(substr(pin,11,11)) %% 2 == 0
  output <- factor(ifelse(female, "Female", "Male"))
  return(output)
}


#' @title
#' Check if \code{pin} is a coordination number
#' 
#' @description
#' Calculate if the personal identity number is a coordination number.
#' 
#' @inheritParams pin_ctrl
#' 
#' @references 
#' \href{https://www.skatteverket.se/download/18.8dcbbe4142d38302d74be9/1387372677724/717B06.pdf}{Population registration in Sweden}
#' \href{https://www.skatteverket.se/download/18.1e6d5f87115319ffba380001857/1285595720207/70408.pdf}{SKV 704}
#' \href{http://www.riksdagen.se/sv/Dokument-Lagar/Utredningar/Statens-offentliga-utredningar/Personnummer-och-samordningsnu_GWB360/}{SOU 2008:60 : Personnummer och samordningsnummer}
#' 
#' @return
#' Logical vector indicating if the pin is a coordination number (\code{TRUE}) or pin (\code{FALSE}).
#'
#' @examples
#' # Examples taken from SKV 704 (see references)
#' ex_pin <- c("196408233234", "196408833224")
#' pin_coordn(ex_pin)
#'
#' @export
pin_coordn <- function(pin) {
  if(!is.pin(pin)) pin <- as.pin(pin)
  as.numeric(substr(pin,7,8)) > 60
}


#' @title
#' Calculate age of \code{pin} for a given date
#' 
#' @description
#' Calculate the age in full years for a given date.
#' 
#' @inheritParams pin_ctrl
#' @param date Date at which age is calculated. If a vector is provided it must be
#'  of the same length as the \code{pin} argument.
#' @param timespan Timespan to use to calculate age. The actual timespans are:
#' \itemize{
#'   \item \code{years} (Default)
#'   \item \code{months}
#'   \item \code{weeks}
#'   \item \code{days}
#' }
#'
#' @references 
#' \href{https://www.skatteverket.se/download/18.1e6d5f87115319ffba380001857/1285595720207/70408.pdf}{SKV 704}
#'   
#' @return
#' Age as an integer vector.
#'
#' @examples
#' # Example with someone born today
#' today_pin <- 
#'   paste(paste(unlist(strsplit(as.character(Sys.Date()),split = "-")), collapse = ""),
#'         "0000",sep="")
#' pin_age(today_pin)
#' 
#' # Examples taken from SKV 704 (see references)
#' ex_pin <- c("196408233234", "186408833224")
#' pin_age(ex_pin, date = "2012-01-01")
#'
#' @export
pin_age <- function(pin, date=Sys.Date(), timespan = "years") {
  if (length(date) == 1) {
    message("The age has been calculated at ", as.character(date), 
            ".")
  } 
  else if (length(date) == length(pin)){
    message("The age is calculated relative to the '", deparse(substitute(date)), "' date")
  }
  else {
    stop("Multiple dates used.")
  }
  
  date <- as.Date(date)
  if(!is.pin(pin)) pin <- as.pin(pin)
  
  all_pins <- pin
  if (length(date) > 1){
    valid_diff <- !is.na(all_pins) & !is.na(date)
  }else{
    valid_diff <- !is.na(all_pins)
  }
  pin <- all_pins[valid_diff]
  
  diff <- lubridate::interval(pin_to_date(pin),
                   lubridate::ymd(date))

  timespan_lubridate <-
    switch(timespan,
           "years" = lubridate::years(1),
           "months" = lubridate::new_period(month=1),
           "weeks" = lubridate::weeks(1),
           "days" = lubridate::days(1))
  
  age <- as.integer(diff %/% timespan_lubridate)
  if(any(age < 0)) warning("Negative age(es).")
  
  all_age <- rep(as.integer(NA), length(all_pins))
  all_age[valid_diff] <- age
  all_age
}


#' @title
#' Calculate the date of birth from a \code{pin}
#' 
#' @description
#' Calculates the date of birth in date format.
#' 
#' @inheritParams pin_ctrl
#' 
#' @return
#' Date of birth as a vector in date format.
#' 
#' @examples
#' # Examples taken from SKV 704 (see references)
#' ex_pin <- c("196408233234", "186408833224")
#' pin_to_date(ex_pin)
#' 
#' @export
pin_to_date <- function(pin) {
  if(!is.pin(pin)) pin <- as.pin(pin)
  pin <- pin_coordn_correct(pin)
  lubridate::ymd(substr(pin,1,8))
}


#' @title
#' Calculate the birthplace of \code{pin}
#' 
#' @description
#' Calculate the birthplace for a given personal identity number born before 1990. See details.
#' 
#' @details
#' It is possible to calculate where people where born (and/or if a person has immigrated) 
#' through their personal identity number. This is possible for people that was born 
#' before 1990 and after 1945. 
#' 
#' For people born before 1946 the birthplace identifier contains information on where
#' one where registered the 1st of november 1946.
#' 
#' Personal identity numbers for people born after 1989 do not contain any information
#' on birthplace.
#' 
#' During the period 1946 - 1989 the pin also contains information on whether one has 
#' immigrated to Sweden during the period.
#' 
#' @inheritParams pin_ctrl
#' 
#' @references
#' \href{http://www.riksdagen.se/sv/Dokument-Lagar/Utredningar/Statens-offentliga-utredningar/Personnummer-och-samordningsnu_GWB360/}{SOU 2008:60 : Personnummer och samordningsnummer}
#' 
#' @return
#' Birthplace as factor.
#'
#' @examples
#' # Example with someone born today and from SKV 704 (see references)
#' today_pin <- paste0(format(Sys.Date(),"%Y%m%d"), "0000")
#' ex_pin <- c("196408233234", today_pin)
#' pin_birthplace(ex_pin)
#'
#' @export
pin_birthplace <- function(pin){
  if(!is.pin(pin)) pin <- as.pin(pin)
  birth_vector <- 
    c(rep("Stockholm stad",10),
      rep("Stockholms l\u00E4n", 4),
      rep("Uppsala l\u00E4n", 2),
      rep("S\u00F6dermanlands l\u00E4n", 3),
      rep("\u00D6sterg\u00F6tlands l\u00E4n", 5),
      rep("J\u00F6nk\u00F6pings l\u00E4n", 3),
      rep("Kronobergs l\u00E4n", 2),
      rep("Kalmar l\u00E4n", 3),
      rep("Gotlands l\u00E4n", 1),
      rep("Blekinge l\u00E4n", 2),
      rep("Kristianstads l\u00E4n", 4),
      rep("Malm\u00F6hus l\u00E4n", 7),
      rep("Hallands l\u00E4n", 2),
      rep("G\u00F6teborgs och Bohus l\u00E4n", 7),
      rep("\u00C4lvsborgs l\u00E4n", 4),
      rep("Skaraborgs l\u00E4n", 3),
      rep("V\u00E4rmlands l\u00E4n", 3),
      rep("Extra number", 1),
      rep("\u00D6rebro l\u00E4n", 3),
      rep("V\u00E4stmanlands l\u00E4n", 2),
      rep("Kopparbergs l\u00E4n", 3),
      rep("Extra number", 1),
      rep("G\u00E4vleborgs l\u00E4n", 3),
      rep("V\u00E4sternorrlands l\u00E4n", 4),
      rep("J\u00E4mtlands l\u00E4n", 3),
      rep("V\u00E4sterbottens l\u00E4n", 4),
      rep("Norrbottens l\u00E4n", 4),
      rep("Extra number and immigrants (immigrated after 1946)", 7))
  birth_other_text <- "Born after 31 december 1989"  
  
  to_na <- pin_coordn(pin)
  to_na[is.na(to_na)] <- TRUE
  
  res <- factor(vapply(X = pin, 
                       FUN = pin_birthplace_internal, 
                       FUN.VALUE = character(1), 
                       birth_vector = birth_vector, 
                       birth_other_text = birth_other_text,
                       USE.NAMES = FALSE), levels = c(unique(birth_vector), birth_other_text))
  res[to_na] <- NA  
  return(res)
}

