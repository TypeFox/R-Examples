## calling handlers

recoverHandler <- function(condition) {
    string1 <- function(what) if(length(what) > 1) paste(what[[1]], "...") else what
    message("Recover being started from a condition of class \"",
            string1(class(condition)), "\":\n", conditionMessage(condition),
            "\nCalled from: ", string1(deparse(conditionCall(condition))))
    recover()
}
