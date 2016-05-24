`as.OtherDate` <-
function(x, calendar,...) UseMethod("as.OtherDate")

as.OtherDate.default<-function(x, calendar,...){
     as.OtherDate(as.Date(x,...), calendar,...)
}
