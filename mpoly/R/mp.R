#' Define a multivariate polynomial.
#' 
#' mp is a smart function which attempts to create a formal mpoly
#' object from a character string containing the usual
#' representation  of a multivariate polynomial.
#' 
#' @param string a character string containing a polynomial, see
#'   examples
#' @param varorder (optional) order of variables in string
#' @return An object of class mpoly.
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @seealso \code{\link{mpoly}}
#' @export mp
#' @examples
#' ( m <- mp("x + y + x y") )
#' is.mpoly( m )
#' unclass(m)
#' 
#' 
#' mp("x + 2 y + x^2 y + x y z") 
#' mp("x + 2 y + x^2 y + x y z", varorder = c("y", "z", "x")) 
#' # mp("x + 2 y + x^2 y", varorder = c("q", "p")) # -> error
#' 
#' ( ms <- mp(c("x + y", "2 x")) )
#' is.mpolyList(ms)
#' 
#' 
#' gradient( mp("x + 2 y + x^2 y + x y z") ) 
#' gradient( mp("(x + y)^10") ) 
#' 
#' # mp and the print methods are kinds of inverses of each other
#' ( polys <- mp(c("x + y", "x - y")) )
#' strings <- print(polys)
#' strings
#' mp(strings)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
mp <- function(string, varorder){
  
  stopifnot(is.character(string))
  
  # if string is a vector of polys, return mpolyList
  if(length(string) > 1){
    if(missing(varorder)){
      mpolyList <- lapply(as.list(string), mp)
    } else {
      mpolyList <- lapply(as.list(string), mp, varorder = varorder)
    }
    class(mpolyList) <- "mpolyList"
    return(mpolyList)
  }
  
  # check for unmatched parentheses
  unmatched_parentheses_stop(string)
  
  # compute
  out <- parse_parenthetical_polynomial(string)
  
  # check varorder argument
  if(!missing(varorder)){
    
    vars <- vars(out)
    
    if(!all(vars %in% varorder)){
      error <- paste(
        "if specified, varorder must contain all computed vars - ",
        paste(vars, collapse = ", "),
        sep = ""
      )
      stop(error, call. = FALSE)
    }
    
    # order vars appropriately
    vars <- intersect(varorder, vars)
    out  <- reorder.mpoly(out, varorder = vars)
  } 
  
  # return
  out
}









# parse_parenthetical_polynomial("x ((x+y)+2)")
# parse_parenthetical_polynomial("x ((x+y) + 2)")
# parse_parenthetical_polynomial("-(x + y) + 2 x (x + y)^2 + 3 y")
parse_parenthetical_polynomial <- function(string){
  
  # trim up
  string <- str_trim(string)
  
  # fix term joins
  terms <- extract_polynomial_terms(string)
  
  # parse into mpolys
  mpolys <- lapply(as.list(terms), parse_parenthetical_term)
  
  # add and return
  Reduce(`+.mpoly`, mpolys)
}



























# parse_parenthetical_term("3 y")
# parse_parenthetical_term("-(x + y)")
# parse_parenthetical_term(" 3 (x + y) 4 (x - y) ")
# parse_parenthetical_term("(x + y) (x - y)")
# parse_parenthetical_term("-(x + y) (x - y)")
# parse_parenthetical_term("(x + y)")
# parse_parenthetical_term("(x + y)^2")
# parse_parenthetical_term("(x + y)(x-y)")
# parse_parenthetical_term("-2 (x + y)^2 (x - y)^0 4 (1+1)^3")
# 
# # more complex usage
# parse_parenthetical_term("((x^2))") 
# parse_parenthetical_term("((5^2))") 
# string <- "(1+1) (2^3 z (x+y)^2)^2"
# parse_parenthetical_term(string)
# parse_parenthetical_term("6 (x)")
# parse_parenthetical_term("6.18033988749895 (x)")
parse_parenthetical_term <- function(string){
  
  # short circuit if simpler
  if(!contains_parenthetical_expression(string)) return(parse_nonparenthetical_term(string))
  
  # break into parenthetical pieces ("bubbles")
  pieces <- term_parentheticals(string)
  
  # mpoly pieces
  mpolys <- lapply(pieces, function(piece){
    
    # identify expression and exponent components
    expr <- str_extract(piece, "\\([()[:alnum:][+-]\\s\\^\\.]+\\)")
    expr <- str_sub(expr, 2, -2) # take off parens
    
    # check for exponent on the outer parenthetical
    last_paren_ndx     <- nchar(piece) - str_locate(str_rev(piece), fixed(")"))[[1]] + 1
    string_after_paren <- str_sub(piece, last_paren_ndx+1) # "" or "^3"
  
    # if "^3", extract, otherwise 1
    if(str_detect(string_after_paren, fixed("^"))){
      exponent <- as.numeric(str_rev(str_extract(str_rev(string_after_paren), "[\\d]+"))) # gets first
    } else {
      exponent <- 1
    }
    
    # parse
    if(contains_nested_parenthetical_expression(piece)){
      parse_parenthetical_polynomial(expr)^exponent  
    } else {
      parse_nonparenthetical_polynomial(expr)^exponent
    }
    
  })
  
  # product and return
  Reduce(`*.mpoly`, mpolys)
}











































# parse_nonparenthetical_polynomial(" -1")
# parse_nonparenthetical_polynomial("x-1")
# parse_nonparenthetical_polynomial("5-2x")
# parse_nonparenthetical_polynomial("5 - 2     x")
# parse_nonparenthetical_polynomial("5 + -2x")
# parse_nonparenthetical_polynomial("1--1")
# parse_nonparenthetical_polynomial("1 - - 1")
# parse_nonparenthetical_polynomial("5^2x")
# parse_nonparenthetical_polynomial("5^2-x")
# parse_nonparenthetical_polynomial("-x")
# parse_nonparenthetical_polynomial("-1")
# parse_nonparenthetical_polynomial("1+-x-x")
# parse_nonparenthetical_polynomial("1 - -3")
#
# parse_nonparenthetical_polynomial("-x + 2y - 4x - -4")
#
# string <- "-4 + 2+2 x +   1 x y^4 -3 prq^3 -y - 3 x 2 - 3 y -2"
# parse_nonparenthetical_polynomial(string)
# parse_nonparenthetical_polynomial("x    +       y")
# parse_nonparenthetical_polynomial("x    -       y+-xy")
# parse_nonparenthetical_polynomial("1e-2 x")
parse_nonparenthetical_polynomial <- function(string){
  
  # trim
  string <- str_trim(string)
  
  # check to see if it's a single term
  if(
    !str_detect(string, "[+]") && 
    !str_detect(str_sub(string, 2), "[-]")
  ) {
    return(parse_nonparenthetical_term(string))
  }
  
  # regularize term joins (deal with minuses)
  string <- fix_term_joins(string)  
  
  # split polynomial
  terms <- str_trim(str_split(string, fixed(" + "))[[1]])
  
  # parse terms
  mpolyTerms <- lapply(as.list(terms), parse_nonparenthetical_term)
  
  # combine and return
  Reduce(`+.mpoly`, mpolyTerms)
}









# parse_nonparenthetical_term("12var 2 y 2x")
# parse_nonparenthetical_term("-2      7")
# parse_nonparenthetical_term("2 x y^2 3 2           3^2")
# parse_nonparenthetical_term("2 x -2") # -> warn
# parse_nonparenthetical_term("x")
# parse_nonparenthetical_term("-x")
# parse_nonparenthetical_term("-5x")
# parse_nonparenthetical_term("-0x")
# parse_nonparenthetical_term("1.5x")
# parse_nonparenthetical_term("1.5^2x")
# parse_nonparenthetical_term("1e-2 x")
parse_nonparenthetical_term <- function(string){
  
  # trim
  string <- str_trim(string)
  
  # fix spaces around exponents "x ^ 2" -> "x^2"
  if(str_detect(string, fixed("^"))){
    string <- str_replace_all(string, "[\\s]*\\^[\\s]*", "^")
  }
  
  # fix spaces around minuses "x ^ 2" -> "x^2"
  if(str_detect(string, fixed("-"))){
    string <- str_replace_all(string, "[\\s]*\\-[\\s]*", "-")
  }
  
  # split based on spaces
  parts <- str_split(string, " ")[[1]]
  parts <- parts[nchar(parts) > 0] # for "2        -2"
  
  # if more than one negative provided error
  if(str_detect(str_sub(string, 2), fixed("-"))){
    stop("negative signs are only allowed at the beginning of terms.", call. = FALSE)
  }
  
  # fix, e.g. "2x"
  smashed_var_bool <- str_detect(parts, "[0-9]+[:alpha:]")
  if(any(smashed_var_bool)){
    places_to_break <- str_locate(parts[smashed_var_bool], "[:alpha:]")[,1]
    for(k in seq_along(places_to_break)){
      parts[smashed_var_bool][k] <- str_c(
        str_sub(parts[smashed_var_bool][k], 1, places_to_break[k]-1),
        "|",
        str_sub(parts[smashed_var_bool][k], places_to_break[k])
      )
    }
    parts <- unlist(str_split(parts, fixed("|")))
  }
  
  # fix, e.g. "-y"
  minus_var_bool <- str_detect(parts, "\\-[:alpha:]")
  if(any(minus_var_bool)){
    parts[minus_var_bool] <- str_c("-1 ", str_sub(parts[minus_var_bool], 2))
    parts <- unlist(str_split(parts, " "))
  }
  
  # collect numeric elements
  parts_with_vars <- str_detect(parts, "[:alpha:]")
  if(all(parts_with_vars)){
    coef <- 1L
  } else {
    coef <- prod(
      vapply(
        as.list(parts[which(!parts_with_vars)]), 
        function(.) eval(parse(text = .)),
        double(1)
      )  
    ) # this multiplies even, e.g., 5^2
  }
  
  # if only coefs are given, return
  if(all(parts_with_vars == FALSE)) return(mpoly(list(c(coef = coef))))
  
  # parse variable exponents
  var_parts <- parts[parts_with_vars]
  var_parts_with_exps_bool <- str_detect(var_parts, fixed("^"))
  var_parts[!var_parts_with_exps_bool] <- str_c(var_parts[!var_parts_with_exps_bool], "^1")
  var_parts <- str_split(var_parts, fixed("^"))
  vars <- vapply(var_parts, `[`, character(1), 1L)
  exps <- as.integer(vapply(var_parts, `[`, character(1), 2L))
  names(exps) <- vars
  
  # mpoly and return
  mpoly(list(c(coef = coef, exps)))  
}













# fix_term_joins("-2 - -2x + y - -3 y - 2")
# fix_term_joins("1-1")
# fix_term_joins("x[1]")
# fix_term_joins("x[1,1]")
# fix_term_joins("1--1")
# fix_term_joins("1 - - 1")
# fix_term_joins("5 - 2     x")
# fix_term_joins("5^2x - 1")
# fix_term_joins("1+-xx-x")
# fix_term_joins("-1-1")
# fix_term_joins("1e-2 x")
# fix_term_joins("1e+2 x")
# fix_term_joins("-1-1-") # error
# fix_term_joins("-1-1+") # error
fix_term_joins <- function(string){
  
  # trim
  string <- str_trim(string)
  
  # make sure last char is not a sign
  if(str_detect(str_sub(string, -1), "[+-]")){
    stop(paste0("term ", string, " does not terminate."), call. = FALSE)
  }
  
  # zero trick for leading symbol, e.g. "-1 + x" -> "0 + -1 + x"
  if(str_detect(str_sub(string, 1, 1), fixed("-"))){
    if(str_detect(str_sub(string, 1, 1), fixed("--"))){
      stop('"--" cannot lead an expression.', call. = FALSE)
    }
    string <- str_c("0 + ", string)
  }
  
  # fix scientific notation
  while(str_detect(string, fixed("e+"))){
    stringToReplace <- str_extract(string, "[0-9\\.]+e\\+[0-9]+")
    string <- str_replace(string, "[0-9\\.]+e\\+[0-9]+", format(as.numeric(stringToReplace)))
  }
  
  while(str_detect(string, fixed("e-"))){
    stringToReplace <- str_extract(string, "[0-9\\.]+e-[0-9]+")
    string <- str_replace(string, "[0-9\\.]+e-[0-9]+", format(as.numeric(stringToReplace)))    
  }
    
  # break string into pieces of terms and joins
  terms <- str_extract_all(string, "[[:alnum:]\\^\\|\\.\\[\\]\\,]+")[[1]]
  joins <- str_split(string, "[[:alnum:]\\^\\.\\[\\]\\,|]+")[[1]]
  if(joins[1] == "") joins <- joins[-1] 
  if(joins[length(joins)] == "") joins <- joins[-length(joins)] 
  if(length(joins) == 0L) return(string)
  
  # fix joins
  pureJoins <- str_replace_all(joins, "\\s", "")
  pureJoins[pureJoins == ""] <- "|"
  if(any(nchar(pureJoins) > 3)) stop(
    "arithmetic sign sequence of more than two detected.", call. = FALSE
  )
  cleanJoinMap <- c("-" = " - ", "+" = " + ", "--" = " + ", 
    "++" = " + ", "+-" = " - ", "-+" = " - ", "|" = " "
  )
  cleanedJoins <- unname(cleanJoinMap[pureJoins]) # cbind(joins, cleanedJoins)

  # reconstruct
  n <- length(terms) + length(joins) # n always odd, first term always a [:alnum:]
  temp <- character(n)
  temp[seq(1, n, 2)] <- terms
  temp[seq(2, n-1, 2)] <- cleanedJoins    
  string <- paste(temp, collapse = "")
  
  # kill minuses
  string <- str_replace_all(string, fixed("-"), "+ -1")
  
  # strip leading "0 + " if needed
  if(str_sub(string, 1, 4) == "0 + ") string <- str_sub(string, 5)
   
  # return
  string
}
















# string <- "-(x + y) + 2 x (x + y) + 3 y"
# extract_polynomial_terms(string)
extract_polynomial_terms <- function(string){
  
  # run fix_term_joins on blanked strings to get protect parentheticals
  blanked_string <- blank_parentheticals(string, "|")
  piped_string <- fix_term_joins(blanked_string)
  
  # change +"s to *"s for breaking later
  # they distinguish polynomial terms
  piped_string <- str_replace_all(piped_string, fixed("+"), "*")
  
  # unprotect
  string_ndcs <- str_locate_all(blanked_string, "[|]+")[[1]]
  piped_ndcs  <- str_locate_all(piped_string,   "[|]+")[[1]]  
  if(nrow(string_ndcs) > 0){
    for(k in 1:nrow(string_ndcs)){
      str_sub(piped_string, piped_ndcs[k,1], piped_ndcs[k,2]) <- 
        str_sub(string, string_ndcs[k,1], string_ndcs[k,2])  
    }
  }
  
  # split
  str_trim(str_split(piped_string, fixed("*"))[[1]])
}






















# an inner parenthetical is one that does not contain parentheticals
# extract_leftmost_inner_parenthetical("(x + 5)")
# extract_leftmost_inner_parenthetical("(x + 5)", contents_only = TRUE)
#
# extract_leftmost_inner_parenthetical("(x + 5)^10")
# extract_leftmost_inner_parenthetical("(x + 5)^10", contents_only = TRUE)
# extract_leftmost_inner_parenthetical("((x + 5)^10+2)^2")
# extract_leftmost_inner_parenthetical("((x + 5)^10+2)", contents_only = TRUE)
extract_leftmost_inner_parenthetical <- function(string, contents_only = FALSE){
  string <- str_extract(string, "\\([^()]*\\)[\\^\\d]*")
  if(contents_only){
    str_sub(string, 2, str_locate(string, fixed(")"))[1,1]-1)
  } else {
    string
  }  
}





# blank_parentheticals(" -1 1 x (3 x + -1 (7 + -1 2 x))^2 7 (x + 1) -3 ")
# blank_parentheticals(" -1 1 x (3 x + -1 (7 + -1 2 x))^2 7 (x + 1) -3 ", "*")
# blank_parentheticals(" -1 1 x (3 x + -1 (7 + -1 2 x))^2 7 (x + 1) -3 ", "_")
blank_parentheticals <- function(string, char = "-"){
  # " -1 1 x (3 x + -1 (7 + -1 2 x))^2 7 (x + 1) -3 " ->
  # " -1 1 x ------------------------- 7 ------- -3 "
  # this blanks parentheticals from the inside out
  # inside parentheticals are done first
  while(contains_parenthetical_expression(string)){
    bad <- extract_leftmost_inner_parenthetical(string)
    string <- str_replace(
      string, 
      "\\([^()]*\\)[\\^\\d]*", 
      str_dup(char, nchar(bad))
    )
  }
  string
}







# string <- " -3 -(x + y)^2 4 (x - y)x 4 "
# collect_nonparenthetical_elements(string)
# string <- " (x + y)^2   (x - y)(x)  "
# collect_nonparenthetical_elements(string)
collect_nonparenthetical_elements <- function(string){
  
  blanked_string <- blank_parentheticals(string, "|")  
  blanked_string <-str_replace_all(blanked_string, fixed("-|"), "-1 |")
  elements <- str_replace_all(blanked_string, "[|]+", "")
  elements <- str_trim(str_replace_all(elements, "[\\s]+", " "))
  if(str_detect(elements, "[:alnum:]")){
    return(str_c("(", elements, ")"))
  } else {
    return("")
  }
  
}

# string <- " -3 (x + y)^2 4 (x - y)x 4 "
# delete_nonparenthetical_elements(string)
# string <- " (x + y)^2   (x - y)(x)  "
# delete_nonparenthetical_elements(string)
# string <- ".2 (x)"
# delete_nonparenthetical_elements(string)
delete_nonparenthetical_elements <- function(string){
  
  blanked_string <- blank_parentheticals(string, "*")  
  erase_ndcs <- str_locate_all(blanked_string, "[[:alnum:][-][//^][.]]")[[1]][,1]
  for(k in erase_ndcs) str_sub(string, k, k) <- " "
  string
  
}




# string <- " -3 (x + y)^2 4 (x - y)(x) 4 "
# term_parentheticals(string)
# string <- " -(x + y)^2   3(x - y)(x)  "
# term_parentheticals(string)
# string <- ".2 (x)"
# term_parentheticals(string)
term_parentheticals <- function(string){
  
  nonpars <- collect_nonparenthetical_elements(string)
  pars <- delete_nonparenthetical_elements(string)
  pars <- str_replace_all(pars, fixed(")("), ") (") # )( -> ) (
  pars <- str_trim(str_replace_all(pars, "[\\s]+", " "))
  
  if(nonpars != ""){
    bubbles <- paste(nonpars, pars)
  } else {
    bubbles <- pars
  }
  
  blank_bubbles <- blank_parentheticals(bubbles)
  real_spaces <- str_locate_all(blank_bubbles, fixed(" "))[[1]][,1]
  for(k in real_spaces) str_sub(bubbles, k, k) <- "|"
  str_split(bubbles, fixed("|"))[[1]]
}










contains_parenthetical_expression <- function(string){ 
  any(str_detect(string, fixed(c("(",")"))))
}




# contains_nested_parenthetical_expression("5+5")
# contains_nested_parenthetical_expression("(5+5)")
# contains_nested_parenthetical_expression("((5+5))")
# contains_nested_parenthetical_expression("x + (5 y) + 2")
# contains_nested_parenthetical_expression("x + ((5 y) + 2)")
contains_nested_parenthetical_expression <- function(string){
  only_parentheses <- str_replace_all(string, "[^()]", "")
  str_detect(only_parentheses, fixed("(("))
}



unmatched_parentheses_stop <- function(string){
  if(contains_parenthetical_expression(string)){
    if(str_count(string, fixed("(")) > str_count(string, fixed(")"))){
      stop("not all parenthetical expressions closed.", call. = FALSE)
    } else if(str_count(string, fixed("(")) < str_count(string, fixed(")"))){
      stop('not all parenthetical expressions closed. (excess )"s detected)', call. = FALSE)
    }
  }
  invisible()
}


str_rev <- function(string) paste(rev.default(str_split(string, "")[[1]]), collapse = "")






