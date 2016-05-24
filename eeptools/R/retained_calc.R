##' Function to calculate whether a student has repeated a grade.
##' @description This function calculates whether or not a student has repeated 
##' a grade. It returns a \code{\link{data.frame}} with the student ID and a 
##' character vector with \code{Y} representing they repeated the grade and 
##' \code{N} that they had not.
##' @param df a data.frame containing minimally a student identifier and their 
##' grade.
##' @param sid a character that indicates the name of the student id attribute in 
##' \code{df}. The default value is \code{sid}.
##' @param grade a character that indicates the name of the student grade attribute in 
##' \code{df}. The default value is \code{grade}.
##' @param grade_val a numeric vector that contains the value of the grade that is 
##' being checked for retention. The default value is \code{9}.
##' @return a data.frame
##' @author Jason P. Becker
##' @export
##' @import data.table
##' @examples
##' x <- data.frame(sid = c(101, 101, 102, 103, 103, 103, 104),
##'                grade = c(9, 10, 9, 9, 9, 10, 10))
##' retained_calc(x)
retained_calc <- function(df, sid = 'sid', grade = 'grade', grade_val = 9){
  df$grade <- df[, grade]
  df$sid <- df[, sid]
  dt <- data.table(df, key=c("sid", "grade"))
  result <- dt[I(grade == grade_val), list(count = .N), by = key(dt)]
  result <- result[, list(sid, retained = ifelse(count > 1, 'Y', 'N'))]
  return(as.data.frame(result))
}
