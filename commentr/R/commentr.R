
################################################################################
#                                                                              #
#                              General help page                               #
#                                                                              #
################################################################################


#' Create nicely formated comments
#'
#' Functions to automate creation of nicely formated comments. 
#' Comments are printed to the console and automaticaly pasted to clipboard 
#' (and can thereafter be copied to script files.
#'
#' The \code{header_comment} function uses global options for name, 
#' contact and textwidth as default. Width is a standard option in R
#' but name and contact are not. It might be
#' a good idea to define those as global options in a .Rprofile-script 
#' (use Google if you do not know the concept!).
#' Dates for creation and last update are both set to the same date by default.
#' The update date has to be manually updated in the script
#'
#' @section Functions:
#' There are 3 functions to help with script comments
#'
#' \itemize{
#'  \item header_comment: a header introduction to be included at the top of the script
#'  \item block_comment: a block comment to make a title for a new section in the script
#'  \item line_comment: a one line comment 
#' }
#'
#' @param title a short and descriptive title.
#' @param description additional information about the script. 
#' Could be a single line or a full paragraph. Empty character string as default.
#' @param author name of the script author(s). From \code{getOption("name")} as default.
#' @param contact How to contact the author of the script. 
#' From \code{getOption("contact")} as default.
#' @param client Name (and possible contact information) to the client for who 
#' the script was written. Sets author as default (assuming you wrote the script for yourself).
#' @param date_created A character string defining the date when the script was created.
#' Todays date as default. The date is specified as a string and does not need to be exact.
#' For example "Many years ago" would be accepted if you want to add a heading to an old script.
#' @param date_updated A character string defining the date when the script was last updated.
#' Same as date_created by default.
#' @param source is the path to where the script was originaly stored . 
#' This information might be good if the document is printed, forgotten and found. 
#' Current directory as default.
#' @param tab "tabulation length" to add space betwen left column titles and right column text.
#' @param empty_lines_first number of empty lines to be added before the text 
#' in a block comment. 1 by default.
#' @param empty_lines_last number of empty lines to be added after the text in a
#'  block comment. Same as \code{empty_lines_first} by default.
#' @param allign specify text allignment. "center" as default with other options 
#' "left" and "right". If the text is too long to be fit on one line, it is 
#' split to a left aligned paragraph (\code{allign} being ignored).
#' @param token a character string specifying the comment character to be used.
#' Default is \code{"#"}, which should be used in R scripts.
#' @param html Should the comment be used in a HTML- or R Markdown-document 
#' (FALSE by default). If so, the comment starts
#' with "<!--" and ends with "-->".
#' \code{token = "\%"} could be used to create comments for LaTeX and 
#' \code{token = "*"} for Stata etcetera. The length of the character string is 
#' not restricted but should normally
#' be just one character. If a longer character string is used, the width of the 
#' comment would be multiplied by the number of characters.
#' This might result in a quite messy output.
#' @param clipboard Should the output be copied to clipboard? (\code{TRUE} by default).
#' Corrently only supported on Mac. 
#' @param verbose Should the comment be printed to the console? 
#' (\code{TRUE} by default). This could be used if \code{clipboard = TRUE} and run on a Mac.
#' @param ... arguments passed to \link{comment_width} to specify the width of the comment.
#' @return There is no object returned from the function call. 
#' There is just a printed message to the console that could be copied to the 
#' beginning of a script).
#' @examples
#' \dontrun{
#'
#' # If global options specifies "author" and "contact", these do not need to be specified every time:
#' header_comment("Test", "This is a little test")
#'
#' header_comment("Test", "This is just a test!", width = "script_width")
#'
#' header_comment("Smaller block", "This is a smaller test block!", width = 55, tab = 17)
#'
#' header_comment("Smaller block", "This is a small test block but with a longer extra description
#'    that has to be split from a single line into a full paragraph.", width = 55, tab = 17)
#' }
#'
#' header_comment("Nice script", "This is a very nice script that is well documented",
#'      author       = "Jane Doe",
#'      contact      = "jane@@doe.se",
#'      client       = "John Doe",
#'      date_created = "2014-07-03",
#'      width        = "a4landscape")
#'
#' block_comment("A title for a new section in the script")
#' block_comment("A shorter box", width = 50)
#' block_comment("A compact title", empty_lines_first = 0, allign = "left")
#' \dontrun{
#' block_comment("A longer descriptive text that has to be
#'    separated into several lines in order to fit.
#'    Then it is no longer alligned to 'center' even if so specified!", allign = "center")
#' }
#' line_comment("A comment in the middle of a line")
#' line_comment("A comment in the middle of a shorter line", 50)
#' @name comment
#' @import stringr
NULL

#' @name commentr
#' @rdname comment
NULL

