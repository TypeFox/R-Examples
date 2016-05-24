summaryStats.matrix <-
function (object, ...) 
{
    summaryStats.data.frame(as.data.frame.matrix(object), ...)
}
