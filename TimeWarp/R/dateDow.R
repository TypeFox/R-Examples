# dateDow format as character with day-of-week appended
dateDow <- function(date)
    UseMethod('dateDow')

dateDow.default <- function(date)
    paste(as.character(date), weekdays(date, abbreviate=TRUE))

dateDow.character <- function(date)
    paste(date, weekdays(dateParse(date, dross.remove=TRUE), abbreviate=TRUE))
