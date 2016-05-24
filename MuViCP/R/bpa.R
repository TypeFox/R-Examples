##Version May 11, 2015
##This is a further optimized version of DS Calculus for singleton focal elements
##This version represents the set of focal elements as a vector instead of a list
##For the entire set (Omega), the special symbol 'Inf' is used

## WARNING ##
##The current implementation specifically DOES NOT work with belief functions
##That is, belief functions cannot be obtained
##THERE IS AN IMPLICIT ASSUMPTION THAT WE WORK WITH ONLY ATOMIC SETS IN THE BPA'S
##This assumption is used to simplify the programs
##This assumption is enough for the purposes of classification that we wish to use DS for
##It also makes the programs faster, since combination is then an O(n) problem
##Atomicity implies bpa = belief
##If your sets are not atomic, you should be wary of using these programs

##This is a set of 2 functions

##The first function, bpa, is a class representing a bpa
##A bpa operates on all possible subsets of a set S
##However, it assigns a positive probability to only some of them, called focal elements
##The setlist is the list of focal elements
##The mlist is the vector of assigned probabilities
##S itself is always a focal element
##A bpa must sum to 1
##In the current implementation, a bpa is either
## --automatically initialized (complete ignorance)
## --completely specified via assign.bpa
##When specifying focal elements via assign.bpa, any residual mass is automatically given to S

##The second function, combine.bpa, works on a pair of objects of class bpa
##It combines the two bpa over their focal elements
##It returns an object of type bpa
##Combination functions removed to separate file (June 11)

bpa <- function(n = 1, setlist = c(1:n, Inf), mlist = c(rep(0,n), 1))
    {
        stopifnot(length(setlist) == length(mlist))
        N <- n
        mlist <- mlist
        names(mlist) <- setlist
        setlist <- setlist
        name <- NULL
        info <- NULL

        assign.bpa <- function(s, m)  
            {
                ##Check for Validity
                stopifnot(length(s) == length(m),
                          all(s %in% setlist),
                          sum(m) <= 1)
                ##Assign
                mlist[setlist %in% s] <<- m
                mlist[N+1] <<- (1 - sum(m))
                names(mlist) <- setlist
            }

        get.mass <- function(s)
            {
                k <- which(setlist %in% s)
                if(length(k) == 0) ##That is element does not occur
                    return(0)
                else
                    return(mlist[k])
            }

        get.assigned.class <- function()
            {
                x <- order(mlist, decreasing = TRUE)
                return(setlist[x])
            }

        ret <- list(get.N = function() N, 
                    get.setlist = function() setlist,
                    get.full.m = function() mlist,
                    get.focal.elements = function() setlist[mlist > 0],
                    get.m = function() mlist[mlist > 0],
                    get.mass = get.mass,
                    assign.bpa = assign.bpa,
                    get.assigned.class = get.assigned.class,
                    get.assigned.ratios = function() return(mlist / mlist['Inf']),
                    set.name = function(name) name <<- name,
                    get.name = function() name,
                    set.info = function(info) info <<- info,
                    get.info = function() info)
        class(ret) <- 'bpa'
        return(ret)
    }

