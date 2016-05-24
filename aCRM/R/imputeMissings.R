imputeMissings <-
function(data) {

Mode <- function (x, na.rm) {
  xtab <- table(x)
  xmode <- names(which(xtab == max(xtab)))
  return(xmode[1])
}

imputed <- data.frame(sapply(data,function(x) {if (class(x) %in% c("character","factor")) x[is.na(x)] <- Mode(x)
                                               else if (class(x) %in% c("numeric","integer")) x[is.na(x)] <- median(x,na.rm=TRUE)
                                               return(x)}, simplify=FALSE))
imputed
}
