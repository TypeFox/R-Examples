summary.influenceSSN <- function(object, ...)
{
    print(object)
}


print.influenceSSN <- function(x, ...)
{
    cat("Object of class influenceSSN\n\n")

    data <- getSSNdata.frame(x$ssn.object)
    data <- data[,"_resid_"]

    summary(data)
}
