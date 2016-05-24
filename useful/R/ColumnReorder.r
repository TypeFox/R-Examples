#' @title colsToFront
#' @description Moves column names to the fron or back of the names
#' @details Moves column names to the fron or back of the names
#' @author Jared P. Lander
#' @export colsToFront
#' @param data data.frame or tbl
#' @param cols Columns that should be moved
#' @return Character vector of column names
#' @examples 
#' theDF <- data.frame(A=1:10, B=11:20, C=1:10, D=11:20)
#' colsToFront(theDF, c('B', 'C'))
#' colsToFront(theDF, c('C', 'B'))
#' colsToFront(theDF, c('C', 'C'))
#' colsToBack(theDF, c('C', 'C'))
#' colsToBack(theDF, c('C', 'B'))
#' colsToBack(theDF, c('C', 'C'))
#' 
colsToFront <- function(data, cols=names(data))
{
    allCols <- names(data)
    # get the columns that are not in cols
    back <- allCols[!allCols %in% cols]
    
    # return the new order
    c(cols, back)
}

#' @title colsToBack
#' @rdname colsToFront
#' @export colsToBack
#' @inheritParams colsToFront
#' 
colsToBack <- function(data, cols=names(data))
{
    allCols <- names(data)
    # get the columns that are not in cols
    back <- allCols[!allCols %in% cols]
    
    # return the new order
    c(back, cols)
}

#' @title moveToFront
#' @description Rearranges column order by moving specified columns to the fron or back.
#' @details Rearranges column order by moving specified columns to the fron or back.
#' @export moveToFront
#' @author Jared P. Lander
#' @importFrom dplyr select_
#' @importFrom magrittr "%>%"
#' @param data data.frame
#' @param cols Character vector specifying the columns to be moved to the front or back
#' @return A data.frame with the columns in the right order
#' @examples 
#' theDF <- data.frame(A=1:10, B=11:20, C=1:10, D=11:20)
#' moveToFront(theDF, c('B', 'C'))
#' moveToFront(theDF, c('C', 'B'))
#' moveToFront(theDF, c('C', 'C'))
#' moveToBack(theDF, c('C', 'C'))
#' moveToBack(theDF, c('C', 'B'))
#' moveToBack(theDF, c('C', 'C'))
#' 
moveToFront <- function(data, cols)
{
    colOrder <- colsToFront(data, cols)
    
    data %>% select_(.dots=colOrder)
}

#' @title moveToBack
#' @rdname moveToFront
#' @export moveToBack
#' @inheritParams moveToFront
#' 
moveToBack <- function(data, cols)
{
    colOrder <- colsToBack(data, cols)
    
    data %>% select_(.dots=colOrder)
}
