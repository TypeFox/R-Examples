#' Does the character vector contain CAS registry numbers? 
#' 
#' Checks that the input contains Chemical Abstract Service registry numbers.
#' 
#' @param x Input to check.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note CAS numbers take the form of 1 to 7 digits followed by a hyphen,  
#' followed by 2 digits, another hyphen and a final check digit.
#' @return A logical vector that is \code{TRUE} when the input contains valid  
#' CAS registry numbers.
#' @examples
#' x <- c(
#'   water            = "7732-18-5", 
#'   d_glucose        = "50-99-7",
#'   l_glucose        = "921-60-8",
#'   no_hyphens       = "7732185", 
#'   two_check_digits = "7732-18-55",
#'   bad_check_digit  = "7732-18-4",
#'   missing          = NA
#' )
#' is_cas_number(x)
#' assert_any_are_cas_numbers(x)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_all_are_cas_numbers(x))
#' @references Chemspider (\url{https://www.chemspider.com}) is a good service 
#' for looking up CAS numbers.
#' @importFrom assertive.base bapply
#' @importFrom assertive.strings character_to_list_of_integer_vectors
#' @export
is_cas_number <- function(x)
{
  x <- coerce_to(x, "character", get_name_in_parent(x))
  #Check format
  rx <- c(
    assertive.strings:::d(1, 7), 
    assertive.strings:::d(2), 
    assertive.strings:::d(1)
  )
  rx <- assertive.strings:::create_regex(rx, sep = "\\-")
    
  ok <- matches <- assertive.strings:::matches_regex(x, rx)

  
  #Check checkdigit
  not_missing_and_ok <- !is.na(ok) & ok
  x[not_missing_and_ok] <- suppressWarnings(assertive.strings:::strip_non_numeric(x[not_missing_and_ok]))
  ok[not_missing_and_ok] <- bapply(
    character_to_list_of_integer_vectors(x[not_missing_and_ok]), 
    function(x)
    {
      lenx <- length(x)
      actual_check_digit <- x[lenx]
      x <- x[-lenx]
      expected_check_digit <- sum(rev(x) * seq_along(x)) %% 10L
      expected_check_digit == actual_check_digit
    }
  )
  set_cause(ok, ifelse(matches, "bad checkdigit", "bad format"))
}

#' Does the character vector contain credit card numbers? 
#' 
#' Checks that the input contains credit card numbers.
#' 
#' @param x Input to check.
#' @param type Type of credit card.  Multiple types can be selected.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note Legacy card numbers, for example 13 digit Visa numbers and 15 digits 
#' JCB numbers are not supported.
#' @return A logical vector that is \code{TRUE} when the input contains valid 
#' credit card numbers.
#' @examples
#' x <- c(
#'   #visa
#'   "4111 1111 1111 1111",    #spaces are allowed where they 
#'                             #would occur on the card
#'   "4012888888881881",       #though they can be omitted
#'   "4111 1111 1111 11111",   #too many digits
#'   "4012888888881882",       #bad check digit
#'   #mastercard
#'   "5555 5555 5555 4444",
#'   "5105 1051 0510 5100",
#'   "5655 5555 5555 4443",    #starts 56
#'   "51051 051 0510 5100",    #bad spacing
#'   #amex
#'   "3782 822463 10005",
#'   "3714 496353 98431",
#'   "3787 344936 71000", 
#'   "3782 822463 1005",       #not enough digits
#'   #diners
#'   "3056 930902 5904",
#'   "3852 000002 3237",
#'   #discover
#'   "6011 1111 1111 1117",
#'   "6011 0009 9013 9424",
#'   #jcb
#'   "3530 1113 3330 0000",
#'   "3566 0020 2036 0505"
#' )
#' is_credit_card_number(x)
#' assert_any_are_credit_card_numbers(x)
#' assertive.base::dont_stop(assert_all_are_credit_card_numbers(x))
#' @references \url{http://www.regular-expressions.info/creditcard.html} 
#' contains the regexes used by this function.
#' The example valid card numbers are from
#' \url{http://www.paypalobjects.com/en_US/vhelp/paypalmanager_help/credit_card_numbers.htm}
#' @importFrom assertive.base bapply
#' @importFrom assertive.strings character_to_list_of_integer_vectors
#' @export
is_credit_card_number <- function(x, 
  type = c("visa", "mastercard", "amex", "diners", "discover", "jcb"))
{
  x <- coerce_to(x, "character", get_name_in_parent(x))
  #Check format
  type <- match.arg(type, several.ok = TRUE)
  d1 <- assertive.strings:::d(1)
  d2 <- assertive.strings:::d(2)
  d3 <- assertive.strings:::d(3)
  d4 <- assertive.strings:::d(4)
  d5 <- assertive.strings:::d(5)
  d6 <- assertive.strings:::d(6)
  
  rx <- list(
    visa       = c(paste0("4", d3), rep.int(d4, 3)),
    mastercard = c(paste0("5[1-5]", d2), rep.int(d4, 3)),
    amex       = c(paste0("3[47]", d2), d6, d5),
    diners     = c("3(0[0-5]|[68][[:digit:]])[[:digit:]]", d6, d4),
    discover   = c(paste("6011", paste0("65", d2), sep = "|"), rep.int(d4, 3)),
    jcb        = c(paste0("35", d2), rep.int(d4, 3))
  )
  rx <- assertive.strings:::create_regex(l = rx[type], sep = " ?")
  ok <- matches <- assertive.strings:::matches_regex(x, rx)
  not_missing_and_ok <- !is.na(ok) & ok
  
  x[not_missing_and_ok] <- suppressWarnings(assertive.strings:::strip_non_numeric(x[not_missing_and_ok]))
  
  #Check check digit with Luhn algorithm
  ok[not_missing_and_ok] <- bapply(
    character_to_list_of_integer_vectors(x[not_missing_and_ok]),
    function(x)
    {
      lenx <- length(x)
      actual_check_digit <- x[lenx]
      x <- rev(x[-lenx])
      doubled <- suppressWarnings(x * 2:1L)
      total <- sum(doubled) - 9 * sum(doubled > 9)
      expected_check_digit <- (9 * total) %% 10L
      expected_check_digit == actual_check_digit
    }   
  )  
  set_cause(
    ok, 
    ifelse(matches, "bad checkdigit", "bad format")
  )
}

#' Does the character vector contain email addresses?
#' 
#' Checks that the input contains email addresses.  (It does not check the the 
#' address exists, merely that the string is in a suitable format.)
#' 
#' @param x Input to check.
#' @param method Name of method to check for validity.  See notes below.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note Each method specifies a regular expression (see 
#' \code{\link[base]{regex}}) to match against. The \code{simple} method matches 
#' most email addresses in use, and is quite good at filtering out typos and 
#' nonsense.  It won't match \emph{every} email address however.  For example, 
#' emails from a top level domain longer than 4 characters won't pass.  The 
#' \code{rfc5322} method implements the official standard for emails.  Thus all 
#' genuine emails will pass, but since the spec is very broad, it isn't as good 
#' at filtering out nonsense.
#' @return A logical vector that is \code{TRUE} when the input contains valid 
#' email addresses.
#' @examples
#' addresses <- c(
#'   ok       = "a@@b.com", 
#'   no_at    = "a_at_b.com", 
#'   no_dot   = "a@@bcom", 
#'   long_tld = "a@@b.comma", 
#'   punct    = "a!$&@@b.com", 
#'   missing  = NA
#' )
#' is_email_address(addresses)
#' is_email_address(addresses, method = "rfc5322")
#' @references \url{http://www.regular-expressions.info/email.html} contains the 
#' regexes used by this function and a good discussion of the pros and cons of 
#' each.
#' @importFrom assertive.base coerce_to
#' @export
is_email_address <- function(x, method = c("simple", "rfc5322"))
{
  method <- match.arg(method)
  x <- coerce_to(x, "character", get_name_in_parent(x))
  rx <- switch(
    method,
    simple = "^[a-z0-9._%+-]+@[a-z0-9.-]+\\.[a-z]{2,4}$",
    rfc5322 = "(?:[a-z0-9!#$%&'*+/=?^_`{|}~-]+(?:\\.[a-z0-9!#$%&'*+/=?^_`{|}~-]+)*|\"(?:[\\x01-\\x08\\x0b\\x0c\\x0e-\\x1f\\x21\\x23-\\x5b\\x5d-\\x7f]|\\\\[\\x01-\\x09\\x0b\\x0c\\x0e-\\x7f])*\")@(?:(?:[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\\.)+[a-z0-9](?:[a-z0-9-]*[a-z0-9])?|\\[(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?|[a-z0-9-]*[a-z0-9]:(?:[\\x01-\\x08\\x0b\\x0c\\x0e-\\x1f\\x21-\\x5a\\x53-\\x7f]|\\\\[\\x01-\\x09\\x0b\\x0c\\x0e-\\x7f])+)\\])"
  )
  ok <- assertive.strings:::matches_regex(x, rx, perl = TRUE)
  set_cause(ok, "bad format")
}

#' Does the character vector contain hex colors?
#'
#' Checks that the input contains hexadecimal colors.
#' 
#' @param x Input to check.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note A string is considered to represent a hexadecimal colour when contains 
#' a hash followed by six hex values.  That is, digits or the letters from a to 
#' f (case insensitive).
#' @return A logical vector that is \code{TRUE} when the input contains hex 
#' colours.
#' @examples
#' x <- c(
#'   "#012345", "#789abc", "#defDEF",  #ok
#'   "012345",                         #no hash
#'   "#g12345",                        #bad letter
#'   "#01 23 45",                      #contains spaces
#'   "#12345", "#1234567"              #wrong length
#' )
#' is_hex_color(x)
#' assert_any_are_hex_colors(x)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_all_are_hex_colors(x))
#' @export
is_hex_color <- function(x)
{
  x <- coerce_to(x, "character", get_name_in_parent(x))
  rx <- assertive.strings:::create_regex("#[0-9a-f]{6}")
  ok <- assertive.strings:::matches_regex(x, rx)
  set_cause(ok, "bad format")
}

#' @rdname is_hex_color
#' @export
is_hex_colour <- is_hex_color

#' Is the string an honorific?
#' 
#' Checks that the input contains honorifics (a.k.a. titles 
#' or salutations).
#' @param x Input to check.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_honorific} returns \code{TRUE} if the input string contains
#' a valid UK postcode. The {assert_*} function returns nothing but throws an 
#' error when the \code{is_*} function returns \code{FALSE}. 
#' @note Single full stops (periods) following a word boundary 
#' and preceding a space or the end of the string are stripped.  
#' Case is ignored.  There is no formal list of official salutations,
#' so this should only be used as a guide, rather than giving a 
#' definitive result.  Especially note that cultural conventions
#' differ across the world and this function has a UK bias.
#' @references Many possibilities borrowed from the Salutation dropdown on
#' the MathWorks account creation page.
#' \url{https://www.mathworks.com/accesslogin/createProfile.do}
#' @examples
#' x <- c("Mr", "MR", "mr.", "Mister", "masTer", "Mr!", "M.r", ".Mr")
#' is_honorific(x)
#' @export
is_honorific <- function(x)
{
  x <- coerce_to(x, "character", get_name_in_parent(x))
  #Strip single dots after words
  x <- gsub("(?<=\\b)\\.(?=\\s|$)", "", x, perl = TRUE)  
  rx <- assertive.strings:::create_regex(
    # standard
    "m([ia]ste)?r", "mrs", "m(is)?s", "d(octo)?r", 
    # academic
    "((assoc)?iate)prof(essor)?", "dean", 
    # religious
    "rev(erend)?", "ft", "father", "bro(ther)?", "s(iste)?r", "(arch)?bishop", 
    # politics and nobility
    "r(igh)?t hon(ourable)", "sir", "lord", "lady", "dame", "prince(ss)?", "king", "queen",
    # military 
    "adm(iral)", "brig(adier )?gen(eral)?", "capt(ain)?", "cdr", "commander", 
    "cpt", "eur ing", "gen(eral)?", "ltc?", "lieutenant( commander| general)?", 
    "lt cdr", "ltjg", "maj(or)?( general)?", "mgen", "radm", "rag", "(staf )?sgt", 
    "sergeant",  "[12]lt",
    # french
    "mlle", "mme", "m",
    # german
    "herr", "frau",
    # dutch
    "heer", "mevrouw",
    # italian
    "sig( ra)?",
    # spanish
    "sr(ta|a)?"
  )
  ok <- assertive.strings:::matches_regex(x, rx)
  set_cause(ok, "bad format")
}

#' Does the character vector contain IP addresses?
#' 
#' Checks that the input contains IP addresses.  (It does not check that the 
#' address exists, merely that the string is in a suitable format.)
#' 
#' @param x Input to check.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note Valid IP addresses are considered to be four integers in the range 0 to 
#' 255, separated by dots, or the string "localhost".
#' @return A logical vector that is \code{TRUE} when the input contains valid IP 
#' addresses.
#' @examples
#' x <- c(
#'   localhost     = "localhost", 
#'   valid_address = "255.0.255.0", 
#'   out_of_range  = "1.2.3.256",
#'   five_blocks   = "1.2.3.4.5",
#'   non_numeric   = "1.2.3.Z",
#'   missing_block = "1.2.3.NA",
#'   missing       = NA
#' )
#' is_ip_address(x)
#' assert_any_are_ip_addresses(x)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_all_are_ip_addresses(x))
#' @importFrom assertive.base bapply
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.strings character_to_list_of_integer_vectors
#' @export
is_ip_address <- function(x)
{
  x <- coerce_to(x, "character", get_name_in_parent(x))  
  rx <- assertive.strings:::create_regex(
    rep.int(assertive.strings:::d(1, 3), 4), 
    sep = "\\."
  )
  ok <- matches <- assertive.strings:::matches_regex(x, rx)
  not_missing_and_ok <- !is.na(ok) & ok
  
  blocks <- strsplit(x[not_missing_and_ok], ".", fixed = TRUE)
  ok[not_missing_and_ok] <- bapply(
    blocks,
    function(b) 
    {  
      all(suppressWarnings(as.integer(b) %in% 0:255))
    }
  )
  
  ok <- ok | (x == "localhost")
  set_cause(ok, ifelse(matches, "big numbers", "bad format"))
}

#' @rdname is_isbn_code
#' @importFrom assertive.base bapply
#' @importFrom assertive.strings character_to_list_of_integer_vectors
#' @export
is_isbn10_code <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  #Check basic format
  rx <- assertive.strings:::create_regex(
    c(rep.int(assertive.strings:::d(1, Inf), 3), "[[:digit:]X]")
  )
  ok <- matches <- assertive.strings:::matches_regex(x, rx)
  not_missing_and_ok <- !is.na(ok) & ok
  
  #Check correct amount of numbers
  x[not_missing_and_ok] <- suppressWarnings(
    assertive.strings:::strip_non_numeric(x[not_missing_and_ok], allow_x = TRUE)
  )
  ok[not_missing_and_ok] <- len <- nchar(x[not_missing_and_ok]) == 10L
  not_missing_and_still_ok <- !is.na(ok) & ok
  
  #Check checkdigit
  ok[not_missing_and_still_ok] <- bapply(
    character_to_list_of_integer_vectors(x[not_missing_and_still_ok]),
    function(x)
    {
      actual_check_digit <- x[10L]
      x <- x[-10L]
      expected_check_digit <- (sum(x * 1:9L) %% 11L)
      if(expected_check_digit == 10L) 
      {
        return(is.na(actual_check_digit))
      }
      expected_check_digit == actual_check_digit
    }
  )  
  set_cause(
    ok, 
    ifelse(
      not_missing_and_still_ok,
      "bad check digit",
      ifelse(not_missing_and_ok, "bad length", "bad format")
    )
  )
}

#' @rdname is_isbn_code
#' @importFrom assertive.base bapply
#' @importFrom assertive.strings character_to_list_of_integer_vectors
#' @export
is_isbn13_code <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  #Check basic format
  rx <- assertive.strings:::create_regex(
    c(rep.int(assertive.strings:::d(1, Inf), 4), "[[:digit:]X]")
  )
  ok <- matches <- assertive.strings:::matches_regex(x, rx)
  not_missing_and_ok <- !is.na(ok) & ok
  
  #Check correct amount of numbers
  x[not_missing_and_ok] <- suppressWarnings(
    assertive.strings:::strip_non_numeric(x[not_missing_and_ok], allow_x = TRUE)
  )
  ok[not_missing_and_ok] <- len <- nchar(x[not_missing_and_ok]) == 13L
  not_missing_and_still_ok <- !is.na(ok) & ok
  
  #Check checkdigit
  ok[not_missing_and_still_ok] <- bapply(
    character_to_list_of_integer_vectors(x[not_missing_and_still_ok]),
    function(x)
    {
      (sum(suppressWarnings(x * c(1, 3))) %% 10L) == 0L
    }
  )  
  set_cause(
    ok, 
    ifelse(
      not_missing_and_still_ok,
      "bad check digit",
      ifelse(not_missing_and_ok, "bad length", "bad format")
    )
  )
}

#' Does the character vector contain ISBN book codes?
#' 
#' Checks that the input contains ISBN-10 or ISBN-13 book codes.
#' 
#' @param x Input to check.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param type Either "isbn10", "isbn13" or both (for matching either type).
#' @param .xname Not intended to be called directly.
#' @return  A logical vector that is \code{TRUE} when the input contains valid 
#' ISBN book codes.
#' @examples
#' x10 <- c(
#'   hyphens             = "0-387-98503-4",
#'   spaces              = "0 387 98503 4",
#'   just_numbers        = "0387985034",
#'   too_long            = "00-387-98503-4",
#'   too_short           = "0-387-9850-4",
#'   non_numeric         = "Z-387-98503-4",
#'   invalid_check_digit = "0-387-98503-5",
#'   missing             = NA
#' )
#' x13 <- c(
#'   hyphens             = "978-0-387-98503-9",
#'   spaces              = "978 0 387 98503 9",
#'   just_numbers        = "9780387985039",
#'   too_long            = "9978-0-387-98503-9",
#'   too_short           = "978-0-387-9850-9",
#'   non_numeric         = "Z78-0-387-9850-9",
#'   invalid_check_digit = "978-0-387-98503-8",
#'   missing             = NA
#' )
#' is_isbn_code(x10, type = "10")
#' assert_any_are_isbn_codes(x10, type = "10")
#' is_isbn_code(x13, type = "13")
#' assert_any_are_isbn_codes(x13, type = "13")
#' #These tests should fail.
#' assertive.base::dont_stop(assert_all_are_isbn_codes(x10, type = "10"))
#' assertive.base::dont_stop(assert_all_are_isbn_codes(x13, type = "13"))
#' @export
is_isbn_code <- function(x, type = c("10", "13"))
{
  type <- match.arg(type, several.ok = TRUE)
  .xname <- get_name_in_parent(x)
  x <- coerce_to(x, "character", .xname)
  ok <- lapply(
    type, 
    function(isbn) 
    {
      fn <- switch(
        isbn,
        "10" = is_isbn10_code,
        "13" = is_isbn13_code
      )
      fn(x, .xname)
    }
  )
  Reduce(`|`, ok)
}

