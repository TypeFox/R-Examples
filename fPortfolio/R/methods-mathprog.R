

################################################################################
# FUNCTION:                     DESCRIPTION:
#  print.solver                  S3 Print method for solver objects
################################################################################


print.solver <- 
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   S3 Print method for solver objects

    # Arguments:
    #   x - a list as returned by solver functions

    # FUNCTION:

    cat("\nMathematical Programming:\n ")
    
    cat(x$version, sep="\n")
    cat("\n")
    
    cat(" Objective:       ", x$objective, "\n")
    cat(" Solution:        ", head(x$solution),  "...\n")
    cat(" Status Code:     ", x$status,    "\n")
    cat(" Message:         ", x$message,   "\n")
    cat(" Solver:          ", x$solver,    "\n")
    
    #cat(" Version:         ", x$version,   "\n")
    
    # Return Value:
    invisible(x)
}


################################################################################

