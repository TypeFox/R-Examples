#' Is the string a valid UK car licence plate number?
#'
#' Checks that the input contains UK car licence plate numbers.
#'
#' @param x Input to check.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note A single space, in the appropriate place, is allowed but not 
#' compulsory.
#' @return \code{is_uk_national_insurance_number} returns \code{TRUE} if the 
#' input string contains a valid UK car licence plate number The {assert_*} 
#' function returns nothing but throw an error when the \code{is_*} function 
#' returns \code{FALSE}.
#' @examples
#' licences <- c(
#'   #1903 to 1931
#'   "A 1", "AA 9999",                   #ok
#'   "A 01",                             #zero prefix on number
#'   "S0", "G0", "RG0", "LM0",           #ok, special plates
#'   #1931 to 1963
#'   "AAA 1", "AAA 999",                 #ok
#'   "III 1", "QQQ 1", "ZZZ 1",          #disallowed letters
#'   "AAA 01",                           #zero prefix on number
#'   #1931 to 1963 alt
#'   "1 AAA", "9999 AAA",                #ok
#'   "1 III", "1 QQQ", "1 ZZZ",          #disallowed letters
#'   "01 AAA",                           #zero prefix on number
#'   #1963 to 1982
#'   "AAA 1A", "AAA 999A",               #ok
#'   "AAA 1I", "AAA 1O", "AAA 1Q",       #disallowed letters
#'   "AAA 1U", "AAA 1Z", 
#'   "AAA 01A",                          #zero prefix on number
#'   #1982 to 2001
#'   "A1 AAA", "A999 AAA",               #ok    
#'   "I1 AAA", "O1 AAA",                 #disallowed letters
#'   "U1 AAA", "Z1 AAA",
#'   "A01 AAA",                          #zero prefix on number
#'   #2001 to 2051
#'   "AA00 AAA", "AA99 AAA",             #ok
#'   "II00 AAA", "QQ00 AAA", "ZZ00 AAA", #disallowed letters
#'   "AA00 III", "AA00 QQQ"
#' )
#' is_uk_car_licence(licences)
#' assert_any_are_uk_car_licences(licences)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_all_are_uk_car_licences(licences))
#' @references Regex taken from 
#' \url{http://www.regexlib.com/REDetails.aspx?regexp_id=527}.
#' @export
is_uk_car_licence <- function(x)
{
  x <- coerce_to(x, "character", get_name_in_parent(x))
  #http://regexlib.com/REDetails.aspx?regexp_id=617
  #http://www.dreamincode.net/code/snippet3031.htm
  one_to_999 <- paste0("[1-9]", assertive.strings:::d(0, 2))
  rx <- assertive.strings:::create_regex(
    `1903 to 1932`         = c("[A-Z]{1,2}", paste0("[1-9]", assertive.strings:::d(0, 3))),
    `1903 to 1932 special` = c("S|G|RG|LM", "0"),
    `1932 to 1963`         = c("[A-HJ-PR-Y]{3}", one_to_999),
    `1932 to 1963 alt`     = c(paste0("[1-9]", assertive.strings:::d(0, 3)), "[A-HJ-PR-Y]{3}"),
    `1963 to 1982`         = c("[A-Z]{3}", one_to_999, "[A-HJ-NPR-TV-Y]"),
    `1983 to 2001`         = c("[A-HJ-NP-TV-Y]", one_to_999, "[A-Z]{3}"),
    `2001 to 2051`         = c(paste0("[A-HJ-PR-Y]{2}", assertive.strings:::d(2)), "[A-HJ-PR-Z]{3}"),
    sep = " ?"
  )
  ok <- assertive.strings:::matches_regex(x, rx)
  set_cause(ok, "bad format")
}
#' @rdname is_uk_car_licence
#' @export
is_uk_car_license <- is_uk_car_licence

#' Is the string a valid UK national insurance number?
#'
#' Checks that the input contains UK national insurance numbers.
#'
#' @param x Input to check.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note A single space is allowed at the appropriate points (after the first 
#' two letters and after each pair of numbers) but not compulsory.
#' @return \code{is_uk_national_insurance_number} returns \code{TRUE} if the 
#' input string contains a valid UK national insurance number.  The {assert_*} 
#' function returns nothing but throw an error when the \code{is_*} function 
#' returns \code{FALSE}.
#' @examples
#' ni_numbers <- c(
#'   "AA 00 00 00 A", "AA 00 00 00", "AA000000A",                #ok
#'   "ZZ 99 99 99 M", "ZZ 99 99 99", "ZZ999999M",                
#'   "DA 00 00 00", "FA 00 00 00", "IA 00 00 00",                #bad first letter
#'   "QA 00 00 00", "UA 00 00 00", "VA 00 00 00",
#'   "AD 00 00 00", "AF 00 00 00", "AI 00 00 00", "AO 00 00 00", #bad second letter
#'   "AQ 00 00 00", "AU 00 00 00", "AV 00 00 00",
#'   "AA 00 00 00 E", "AA 00 00 00 G", "AA 00 00 00 H",          #bad final letter
#'   "AA 00 00 00 I", "AA 00 00 00 J", "AA 00 00 00 K",
#'   "AA 00 00 00 L", "AA 00 00 00 N", "AA 00 00 00 O",
#'   "AA 00 00 00 P", "AA 00 00 00 Q", "AA 00 00 00 R",
#'   "AA 00 00 00 S", "AA 00 00 00 T", "AA 00 00 00 U",
#'   "AA 00 00 00 V", "AA 00 00 00 W", "AA 00 00 00 X",
#'   "AA 00 00 00 Y", "AA 00 00 00 Z"    
#' )
#' is_uk_national_insurance_number(ni_numbers)
#' assert_any_are_uk_national_insurance_numbers(ni_numbers)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_all_are_uk_national_insurance_numbers(ni_numbers))
#' @references Regex taken from 
#' \url{http://www.regexlib.com/REDetails.aspx?regexp_id=527}.
#' @export
is_uk_national_insurance_number <- function(x)
{
  x <- coerce_to(x, "character", get_name_in_parent(x))
  rx <- assertive.strings:::create_regex(
    c(
      "[A-CEGHJ-PR-TW-Z]{1}[A-CEGHJ-NPR-TW-Z]{1}", 
      rep.int(assertive.strings:::d(2), 3), 
      "[A-DFM]?"
    ),
    sep = " ?"
  )
  ok <- assertive.strings:::matches_regex(x, rx)
  set_cause(ok, "bad format")
}

#' Is the string a valid UK postcode?
#' 
#' Checks that the input contains UK postcodes.
#' 
#' @param x Input to check.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_uk_postcode} returns \code{TRUE} if the input string 
#' contains a valid UK postcode. The {assert_*} function returns nothing but 
#' throws an error when the \code{is_*} function returns \code{FALSE}.
#' @note The function doesn't guarantee that the postcode actually exists.  It 
#' should correctly return \code{TRUE} for genuine postcodes, and will weed out 
#' most badly formatted strings and non-existent areas, but some non-existent 
#' districts may incorrectly return \code{TRUE}.  If you need 100% accuracy, 
#' check against an up-to-date postcode database.
#' @examples
#' postcodes <- c(
#'   "SW1A 1AA", "SK11 9DW", "M34FP", "Le45ns", "TS25 2BZ", "gir 0aa"
#' )
#' is_uk_postcode(postcodes)
#' assert_all_are_uk_postcodes(postcodes)
#' @references Regexes taken from 
#' \url{https://en.wikipedia.org/wiki/Postcodes_in_the_United_Kingdom#Validation}.
#' @export
is_uk_postcode <- function(x)
{
  x <- coerce_to(x, "character", get_name_in_parent(x))
  #Alternative regex, not used, at 
  #http://www.regexlib.com/REDetails.aspx?regexp_id=1064  
  standard_area <- "(A[BL]|B[ABDHLNRSTX]?|C[ABFHMORTVW]|D[ADEGHLNTY]|E[HNX]?|F[KY]|G[LUY]?|H[ADGPRSUX]|I[GMPV]|JE|K[ATWY]|L[ADELNSU]?|M[EKL]?|N[EGNPRW]?|O[LX]|P[AEHLOR]|R[GHM]|S[AEGKLMNOPRSTY]?|T[ADFNQRSW]|UB|W[ADFNRSV]|YO|ZE)[1-9]?[0-9]"
  london_area <- "((E|N|NW|SE|SW|W)1|EC[1-4]|WC[12])[A-HJKMNPR-Y]|(SW|W)([2-9]|[1-9][0-9])|EC[1-9][0-9]"
  district <- "[0-9][ABD-HJLNP-UW-Z]{2}"
  
  rx <- assertive.strings:::create_regex(    
    c(standard_area, district),
    c(london_area, district),
    c("GIR", "0AA"),
    sep = " ?"
  )
  ok <- assertive.strings:::matches_regex(x, rx)
  set_cause(ok, "bad format")
}

#' Is the string a valid UK telephone number?
#' 
#' Checks that the input contains UK telephone numbers.
#' 
#' @param x Input to check.
#' @return \code{is_uk_telephone_number} returns \code{TRUE} if the input string 
#' contains a valid UK telephone number. The {assert_*} function returns nothing 
#' but throws an error when the \code{is_*} function returns \code{FALSE}.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note The function doesn't guarantee that the phone number is in use, but 
#' checks that the format is correct, and that the area code exists.
#' Spaces, hyphens and round brackets are allowed to appear in arbitrary places.  
#' The international UK prefix of 0044 or +44 is allowed.
#' @examples
#' phone_nos <- c("+44 207 219 3475", "08457 90 90 90")
#' is_uk_telephone_number(phone_nos)
#' assert_all_are_uk_telephone_numbers(phone_nos)
#' @references The regex is adapted from one on the now defunct 
#' aa-asterisk.org.uk site with some additional consultation from
#' \url{https://en.wikipedia.org/wiki/List_of_United_Kingdom_dialling_codes}
#' @export
is_uk_telephone_number <- function(x)
{
  x <- coerce_to(x, "character", get_name_in_parent(x))
  #Spaces and round brackets appear in arbitrary places; ignore them.
  x <- suppressWarnings(
    assertive.strings:::strip_invalid_chars(x, invalid_chars="[ -()]")
  )
  
  #All numbers should begin with 0 or the country code, 0044. Check and remove.
  start <- "(0|0044|\\+44)"
  d1 <- assertive.strings:::d(1)
  d2 <- assertive.strings:::d(2)
  d3 <- assertive.strings:::d(3)
  d4 <- assertive.strings:::d(4)
  d6 <- assertive.strings:::d(6)
  d7 <- assertive.strings:::d(7)
  d8 <- assertive.strings:::d(8)
  
  first_rx <- assertive.strings:::create_regex(
    c(start, d7, assertive.strings:::d(2, 3, TRUE)), #country prefix + 7, 9 or 10 digits
    sep = ""
  )
  ok <- assertive.strings:::matches_regex(x, first_rx) 
  not_missing_and_ok <- !is.na(ok) & ok
  
  #remove country code prefix
  x[not_missing_and_ok] <- sub(paste0("^", start), "", x[not_missing_and_ok]) 
  
  regional <- paste0("[2-9]", assertive.strings:::d(4, 5))
  second_rx <- assertive.strings:::create_regex(
    #new style city
    c("20[01378]", d7),
    c("23[0189]", d7),
    c("24[017]", d7),
    c("28[0-46-9]", d7),
    c("29[012]", d7),
    #old style city
    c("113[0-48]", d6),
    c("11[46][0-4]", d6),
    c("115[012789]", d6),
    c("117[0-39]", d6),
    c("118[01349]", d6),
    c("121[0-7]", d6),
    c("131[0-8]", d6),
    c("1[459]1", d7),
    c("161[0-46-9]", d6),
    #regional (4+6)
    c("120[024-9]", d6),
    c("122[3-9]", d6),
    c("123[3-79]", d6),
    c("124[1-689]", d6),
    c("12[58][02-9]", d6),
    c("126[0-4789]", d6),
    c("127[013-9]", d6),
    c("129", d7),
    c("130", d7),
    c("13[25][02-9]", d6),
    c("133[02-579]", d6),
    c("13[468][0-46-9]", d6),
    c("137[1235679]", d6),
    c("139[24578]", d6),
    c("140[03-9]", d6),
    c("142[02-5789]", d6),
    c("14[37]", d7),
    c("144[02-69]", d6),
    c("145[0-8]", d6),
    c("14[69][0-79]", d6),
    c("150[1235-9]", d6),
    c("152[024-9]", d6),
    c("153[0145689]", d6),
    c("154[02-9]", d6),
    c("155[03-9]", d6),
    c("156", d7),
    c("157[0-35-9]", d6),
    c("158[0-468]", d6),
    c("159[0-5789]", d6),
    c("160[034689]", d6),
    c("162[0-689]", d6),
    c("16[38][013-9]", d6),
    c("164[1-467]", d6),
    c("165[0-69]", d6),
    c("166[13-9]", d6),
    c("167[0-8]", d6),
    c("169[0124578]", d6),
    c("170[0246-9]", d6),
    c("172", d7),
    c("173[023678]", d6),
    c("174[03-9]", d6),
    c("175[0-46-9]", d6),
    c("176[013-9]", d6),
    c("177[0-35-9]", d6),
    c("178[024-9]", d6),
    c("179[02-9]", d6),
    c("180[35-9]", d6),
    c("182[1-5789]", d6),
    c("183[02-578]", d6),
    c("184[0-578]", d6),
    c("185[124-9]", d6),
    c("186[2-69]", d6),
    c("187", d7),
    c("188[02-9]", d6),
    c("189[02569]", d6),
    c("190[02-589]", d6),
    c("192[02-689]", d6),
    c("193[1-5789]", d6),
    c("194[2-9]", d6),
    c("195[0-579]", d6),
    c("196[234789]", d6),
    c("197[0124578]", d6),
    c("198", d7),
    c("199[2-57]", d6),
    #regional (6+3)
    c("12046[1-4]", d3),
    c("12087[2-9]", d3),
    c("12545[1-79]", d3),
    c("12762", d4),
    c("12763[1-8]", d3),
    c("12766[1-6]", d3),
    c("12972[0-4]", d3),
    c("12973[2-5]", d3),
    c("12982[2-8]", d3),
    c("12987[0-4789]", d3),
    c("12988[345]", d3),
    c("13638[2-5]", d3),
    c("13647[23]", d3),
    c("13847[04-9]", d3),
    c("13864[015789]", d3),
    c("14044[1-7]", d3),
    c("14202[23]", d3),
    c("14208", d4),
    c("146030", d3),
    c("14605[2-57]", d3),
    c("14606[1-8]", d3),
    c("14607[2-8]", d3),
    c("146140", d3),
    c("148052", d3),
    c("14887[123]", d3),
    c("15243[2-79]", d3),
    c("15246", d4),
    c("15276", d4),
    c("15626[06-9]", d3),
    c("156686", d3),
    c("16064", d4),
    c("16067[4-79]", d3),
    c("16295[567]", d3),
    c("1635[34]", d4),
    c("164724", d3),
    c("164761", d3),
    c("16595[08]", d3),
    c("16596[67]", d3),
    c("165974", d3),
    c("16955[0-4]", d3),
    c("17266[13-9]", d3),
    c("17267[0-7]", d3),
    c("17442", d4),
    c("17502[0-3]", d3),
    c("1750[3-68]2", d3),
    c("175076", d3),
    c("1827[56]", d4),
    c("18375[2-5]", d3),
    c("18378[239]", d3),
    c("18843[2-58]", d3),
    c("19006[1-8]", d3),
    c("190085", d3),
    c("19052", d4),
    c("193583", d3),
    c("19466[1-8]", d3),
    c("19492[01]", d3),
    c("194981", d3),
    c("196323", d3),
    c("19633[1-4]", d3),
    c("199561", d3),    
    #special regional
    c("176888[234678]", d2),
    c("16977[23]", d3),  
    #mobiles
    c("7[1-4]", d8),
    c("750[0-8]", d6),
    c("75[13-9]", d7),
    c("752[0-35-9]", d6),
    c("7624", d6),
    c("770[1-9]", d6),
    c("77[1-7]", d7),
    c("778[02-9]", d6),
    c("779[0-689]", d6),
    c("78[014-9]", d7),
    c("78[23][0-8]", d6),
    c("79[04-9]", d7),
    c("791[02-9]", d6),
    c("792[0-35-9]", d6),
    c("793[0-689]", d6),
    #pagers
    c("760[012]", d6),
    c("762[356]", d6),
    c("764[0134]", d6),
    c("765[49]", d6),
    c("766[0-369]", d6),
    c("7677", d6),
    c("7681", d6),
    c("769[39]", d6),
    #free
    "8001111",
    c("800", assertive.strings:::d(6,7)),         
    c("808", d7),
    c("500", d6),
    #premium
    c("87[123]", d7),
    c("9[01]", d8),
    c("98[123]", d7),
    #shared
    c("845464", d1),
    c("84[2-5]", d7),
    c("870", d7),
    #personal
    c("70", d8),
    #VoIP
    c("56", d8),
    #UAN
    c("55", d8),
    c("3[0347]", d8),
    sep = ""
  )
  ok[not_missing_and_ok] <- assertive.strings:::matches_regex(x[not_missing_and_ok], second_rx)
  ok[!is.na(x) & x == "999"] <- TRUE  #Emergency number
  set_cause(ok, "bad format")
}
