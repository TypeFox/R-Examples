#' Tidy data frame of Jane Austen's 6 completed, published novels
#' 
#' Returns a tidy data frame of Jane Austen's 6 completed, published novels with 
#' two columns: \code{text}, which contains the text of the novels divided into 
#' elements of up to about 70 characters each, and \code{book}, which contains the titles of
#' the novels as a factor in order of publication.
#' 
#' @return A data frame with two columns: \code{text} and \code{book}
#' 
#' @name austen_books
#' 
#' @import dplyr
#' 
#' @examples 
#' 
#' library(dplyr)
#'
#' austen_books() %>% group_by(book) %>%
#'      summarise(total_lines = n())
#'
#' @export

austen_books <- function(){
        ret <- bind_rows(data_frame(text = sensesensibility, book = "Sense & Sensibility"),
                  data_frame(text = prideprejudice, book = "Pride & Prejudice"),
                  data_frame(text = mansfieldpark, book = "Mansfield Park"),
                  data_frame(text = emma, book = "Emma"),
                  data_frame(text = northangerabbey, book = "Northanger Abbey"),
                  data_frame(text = persuasion, book = "Persuasion"))
        ret <- mutate(ret, book = factor(book, levels = unique(book)))
        ret
}