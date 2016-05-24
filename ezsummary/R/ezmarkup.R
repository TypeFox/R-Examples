#' Easy way to "markup" a table before it is sent to be displayed
#'
#' @description The final step of an analysis is to export the tables generated
#' through analytical scripts to the desired platform, such as a pdf, a rmarkdown
#' document or even a Shiny page. Sometimes, a lot of people wants to reorganize
#' the data in a more human readable format. For example, people like to put the
#' standard deviation inside a pair of parentheses after the mean. People also
#' like to put the low and high ends of confidence interval inside a pair of
#' parenthese, separated by " ~ ", after the estimated average. However, as far
#' as I know, so far there isn't a straight forward function to deal with this
#' need. This function is built to address this issue.
#'
#' @param tbl The input table
#' @param pattern The grouping pattern. Each dot "." represent one column. If
#' two or more columns need to be combined in some certain formats, they should
#' be put together inside a pair of brackets "[ ]". You can add any special
#' characters, such as "(" and "~", inside the pair of brackets but please
#' don't leave those special characters outside the brackets. If you want to
#' add a dot as a special character. Please use "^.^" for every single dot
#' you would like to add.
#'
#' @examples
#' library(dplyr)
#' dt <- mtcars %>% group_by(cyl) %>% select(gear, carb) %>% ezsummary_categorical(n=TRUE)
#'
#' ezmarkup(dt, "...[.(.)]")
#' ezmarkup(dt, "..[. (. ~ .)]")
#'
#' @export
ezmarkup <- function(tbl, pattern){
  if(!is.character(pattern))stop('Please provide a vilid character string for the markup pattern. For example, "...[.(.)]". ')
  if((nchar(pattern) - nchar(gsub("\\.","", pattern)) - (nchar(pattern) - nchar(gsub('`\\.`', '', pattern)))/3) != ncol(tbl))stop('The number of dots(.) you entered does not match up with the number of columns in the table. Please review your pattern expression. \n Note: If you do want to Add an actual dot, please use `.` to denote that specific dot')
  if(nchar(gsub("\\.", "", gsub("\\[.*?\\]", "", pattern))) != 0) stop ('Please do not put special symbols outside of your bracket')
  general_pattern <- strsplit(pattern, "")[[1]]
  table_export <- data.frame(row.names = 1:nrow(tbl))
  table_export <- ezmarkup_(tbl=tbl, general_pattern = general_pattern, table_export = table_export)
  return(table_export)
}

ezmarkup_ <- function(tbl, general_pattern, table_export){
  if(length(general_pattern) == 0){return(table_export)}
  if(general_pattern[1] == "."){
    table_export <- cbind(table_export, tbl[1])
    general_pattern <- general_pattern[-1]
    tbl <- tbl[-1]
    ezmarkup_(tbl = tbl, general_pattern = general_pattern, table_export = table_export)
  }else{
    if(general_pattern[1] == "["){
      column_export <- ezmarkup_bracket(tbl = tbl, bracket_pattern = general_pattern[-1])
      table_export <- cbind(table_export, column_export)
      names(table_export)[ncol(table_export)] <- attributes(column_export)$column_name_export
      number_of_dots <- attributes(column_export)$number_of_dots
      number_of_characters <- attributes(column_export)$number_of_characters
      general_pattern <- general_pattern[(-1):(-number_of_characters)]
      tbl<-tbl[,c((-1):(-number_of_dots))]
      ezmarkup_(tbl = tbl, general_pattern = general_pattern, table_export = table_export)
    }
  }
}

ezmarkup_bracket <- function(tbl, bracket_pattern, column_export="", column_name_export="", number_of_dots=0, number_of_characters = 1){
  if(length(bracket_pattern) == 0)stop("Whenever you typed a [, please don't forget to close it by a ]")
  if(bracket_pattern[1] == "]") {
    number_of_characters <- number_of_characters + 1
    attributes(column_export)$number_of_dots <- number_of_dots
    attributes(column_export)$number_of_characters <- number_of_characters
    attributes(column_export)$column_name_export <- column_name_export
    return(column_export)
  }
  if(bracket_pattern[1] == ".") {
    column_export <- paste0(column_export, tbl[[1]])
    column_name_export <- paste0(column_name_export, names(tbl)[1])
    bracket_pattern <- bracket_pattern[-1]
    tbl <- tbl[-1]
    number_of_characters <- number_of_characters + 1
    number_of_dots <- number_of_dots+1
    ezmarkup_bracket(tbl = tbl, bracket_pattern = bracket_pattern, column_export = column_export, column_name_export=column_name_export, number_of_dots = number_of_dots, number_of_characters = number_of_characters)
  }else{
  if(bracket_pattern[1] == "^") {
    if (bracket_pattern[2] == "." & bracket_pattern[3] == "^"){
      column_export <- paste0(column_export, ".")
      column_name_export <- paste0(column_name_export, ".")
      bracket_pattern <- bracket_pattern[c(-1, -2, -3)]
      number_of_characters <- number_of_characters + 1
      ezmarkup_bracket(tbl = tbl, bracket_pattern = bracket_pattern, column_export = column_export, column_name_export=column_name_export, number_of_dots = number_of_dots, number_of_characters = number_of_characters)
    } else {
      column_export <- paste0(column_export, "^")
      column_name_export <- paste0(column_name_export, "^")
      bracket_pattern <- bracket_pattern[-1]
      number_of_characters <- number_of_characters + 1
      ezmarkup_bracket(tbl = tbl, bracket_pattern = bracket_pattern, column_export = column_export, column_name_export=column_name_export, number_of_dots = number_of_dots, number_of_characters = number_of_characters)
    }
  }
  column_export <- paste0(column_export, bracket_pattern[1])
  column_name_export <- paste0(column_name_export, bracket_pattern[1])
  bracket_pattern <- bracket_pattern[-1]
  number_of_characters <- number_of_characters + 1
  ezmarkup_bracket(tbl = tbl, bracket_pattern = bracket_pattern, column_export = column_export, column_name_export=column_name_export, number_of_dots = number_of_dots, number_of_characters = number_of_characters)}
}
