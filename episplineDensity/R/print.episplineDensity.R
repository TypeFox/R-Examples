print.episplineDensity <- function (x, ...) 
{
cat ("Epispline Density Estimate\n")
status <- as.character(x$status)
if (status >= 0)
    cat ("Success: status is ", status, "\n") 
else cat ("Failure: status is ", status, "\n")
### cat ("\t", x$message, "\n")
}
