library(gtools)

##############################
### Examples from man page ###
##############################

### Toy examples

# search for x=10
s <- binsearch( function(x) x-10, range=c(0,20) )
stopifnot(s$where==10)

# search for x=10.1
s <- binsearch( function(x) x-10.1, range=c(0,20) )
stopifnot( s$where==c(10,11) )

### Classical toy example

# binary search for the index of 'M' among the sorted letters
fun <- function(X) ifelse(LETTERS[X] > 'M', 1,
                          ifelse(LETTERS[X] < 'M', -1, 0 ) )

s = binsearch( fun, range=1:26 )
stopifnot( LETTERS[s$where]=="M")

##################################
### Test boundary contiditions ###
##################################

s = binsearch(fun = function(x) x-10, range=c(1,10) )
with(s, stopifnot(where==10, value==0, flag=="Found") )

s = binsearch(fun = function(x) x-1, range=c(1,10) )
with(s, stopifnot(where==1, value==0, flag=="Found") )


checkWarning <- function( expr )
    {
        myEnv <- environment()

        catchWarning <- function(w) {
            assign("warningValue", w, pos=myEnv)
            invokeRestart("muffleWarning")
        }

        retval <- withCallingHandlers(expr = expr,
                                      warning = catchWarning)


        if( !exists("warningValue", where=myEnv, inherits=FALSE) )
            stop("Expected a warning message")
    }

checkWarning( s <- binsearch(fun = function(x) x-10, range=c(1,9) ) )
with(s, stopifnot(where==9, value==-1, flag=="Upper Boundary" ) )

checkWarning( s <- binsearch(fun = function(x) x-1, range=c(2,10) ) )
with(s, stopifnot(where==2, value==1, flag=="Lower Boundary" ) )






