"checkList" <-
function (classList, classObjects, Args = NULL, list = Args[[1]], 
    txt = "in initialize") 
{
    "myPossibleExtends" <- function(class1, class2) {
        if (class1 == class2) 
            return(TRUE)
        else return(!(class(possibleExtends(class1, class2)) == "logical"))
    }
    result <- list()
    if (myPossibleExtends(class(list), "list")) {
        if (all(unlist(lapply(list, 
                  function(i) myPossibleExtends(class(i), classObjects))))) 
            result <- list
        else message(paste("Invalid list", txt, "'", classList, "'"))
        if (length(Args) > 1) 
            message(paste("Ignored arguments", txt, "'", classList, "'"))
    }
    else message(paste("Invalid argument", txt, "'", classList, "'"))
    return(result)
}
