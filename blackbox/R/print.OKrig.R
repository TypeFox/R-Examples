print.OKrig <-function(x,...) {summary(x,...)}
print.OKriglistplus <-function(x,...) {summary(x,...)}

summary.OKrig <- function (object, digits = 4, ...) {
  cat("OKrig object with elements:",paste(names(object),collapse=", "),"\n")
  cat("basic info:\n")
  basicinfo <- c(index=object$CKrigidx, 
                 Npoints=nrow(object$x), 
                 NuniqueX=object$nuniquerows
  )
  print(basicinfo)
  invisible(object)
}

summary.OKriglistplus <- function(object, ...) {
  sapply(object$Kriglist,summary.OKrig) 
  blob <- object
  blob$Kriglist <- NULL
  sapply(blob,summary) 
  invisible(object)
}
