## Makes julian days from month and day number
## m: month number from 1 to 12
## d: day number within the month
## y: year
dayofyear <- function(m, d, y=1960){
as.numeric(mdy.date(m, d, y))-as.numeric(mdy.date(1,1,y))+1
}
