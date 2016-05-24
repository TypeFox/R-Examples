print.csFSSgr <-
function (x, ...) 
{
    cat("Class 'csFSSgr' : Graphics derived from csFSS obect:\n")
    cat("       ~~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("Contains scaling and hierarchical clustering solutions\n")
    cat("The contained csFSS object's details follow...\n")
    print.csFSS(x$x)
}
