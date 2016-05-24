#' @S3method print Roc
print.Roc <- function(x,digits=2,...){
    summary(x,digits=digits,print.it=TRUE,...)
}
