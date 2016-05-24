llgp <- function( coefficients, targets, achievements, maxiter=1000, verbose=FALSE )
{
###
### This function minimizes \eqn{ a'=[g_1(n,p), g_2(n,p), ..., g_K(n,p)] } subjective
### to C x + n - p = b, x >= 0, n >= 0 and p >= 0
###
### Parameters
### coefficients = a matrix with the coefficients of the linear objective functions
### targets = a vector of target values for the objective function
### achievements = a data frame with the weights of the deviation variables for each
###                objective along with the corresponding priority level
### maxiter = maximum number of iterations
### zero = number smaller than this value (in absolute terms) are set to zero
### verbose = an optional logic variable to indicate whether interm results are to be printed
###
### validate the argument objects
###
###
### create the tableau
###
    tab <- llgptab( coefficients, targets, achievements )
###
### reset the print and iteration countersls()
###
    prnt <- 0
    tab$iter <- 0
###
### check tableau for negative RHS target values and repair if necessary
###
    check.tb( tab )
###
### loop over priority levels
###
    for ( k in 1:tab$levels ) {
###
###     update the level in the tableau
###
        tab$level <- k
###
###    calculate the index rows for levels 1 to k
###
        tab <- calc.ti.k( tab, k )
###
###     calculate the achievements for levels 1 to k
###
        tab <- calc.ta.k( tab, k )
###
###     infinite loop while there a possibility of converging to a solution 
###
        sp <- ev.llgp( tab, k )
        while ( sp != 0 ) {
            tab$iter <- tab$iter + 1
            if ( tab$iter >= maxiter ) {
                prnt <- prnt + 1
                cat( paste( "Algorithm did not finish", tab$iter, "iterations at level", k ) )
                cat( "\n" )
                print( tab )
                out <- llgpout( tab, coefficients, targets )
                result <- list( tab=tab, out=out, converged=FALSE )
                class( result ) <- "llgp"
                return( result )
            }
###
###         get the index of the departing variable
###
            ip <- dv.llgp( tab, sp )
            if ( ip == 0 ) {
                cat( paste( "Failed pivot computation at level", k ) )
                cat( "\n" )
                prnt <- prnt + 1
                print( tab )
                out <- llgpout( tab, coefficients, targets )
                result <- list( tab=tab, out=out, converged=FALSE )
                class( result ) <- "llgp"
                return( result )
            }
###
###         swap the entering and departing variables
###
            tab <- piv.llgp( tab, sp, ip, verbose )
###
###         update the index rows and achievement functions
###
            tab <- calc.ti.k( tab, k )
            tab <- calc.ta.k( tab, k )
            sp <- ev.llgp( tab, k )
            if ( verbose ) print( tab )
        }
    }
    out <- llgpout( tab, coefficients, targets )
    result <- list( tab=tab, out=out, converged=TRUE )
    class( result ) <- "llgp"
    return( result )
}
