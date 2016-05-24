recoverHandler <- function(condition) {
    string1 <- function(what) if(length(what) > 1)
      paste(what[[1]], "...") else what
    message("A condition of class \"",
            string1(class(condition)), "\" occurred, with message:\n",
            conditionMessage(condition))
    call <- conditionCall(condition)
    if(!is.null(call))
      message(
            "The condition occurred in: ", string1(deparse(call)))
    recover()
}
