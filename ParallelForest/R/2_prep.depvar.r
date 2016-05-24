#------------------------------------------------------------------------------
#   Copyright (C) 2014  Bertram Ieong
#   No warranty provided.
#------------------------------------------------------------------------------

prep.depvar.in = function(y) {
    # Given the dependent variable vector with two classes, this function 
    # returns a list that contains the dependent variable vector expressed
    # with only integer 0s and 1s and information needed to convert it back.
    # This will be a helper function to any grow functions.


    if(length(class(y))==1 && (class(y)=="numeric" || class(y)=="integer" || class(y)=="character" || class(y)=="factor")){
        if(class(y)=="factor"){
            y.underlying = sort(unique(as.integer(y)))
        } else {
            y.underlying = sort(unique(y))
        }

        if(any(is.na(y.underlying))) stop("NA cannot be a class in the dependent variable")
        if(length(y.underlying)!=2) stop("Dependent variable must have exactly two classes")

        if(class(y)=="integer"){
            if(all(y.underlying==c(0,1))){
                return(list(y.prepped=as.integer(y)))
            }
        }

        y.prepped = vector(mode="integer", length(y))

        if(class(y)=="factor"){
            y.prepped[as.integer(y)==y.underlying[1]] = 0
            y.prepped[as.integer(y)==y.underlying[2]] = 1

            return(list(y.prepped=as.integer(y.prepped), y.underlying=y.underlying, factor.levels=levels(y)))

        } else {
            y.prepped[y==y.underlying[1]] = 0
            y.prepped[y==y.underlying[2]] = 1

            return(list(y.prepped=as.integer(y.prepped), y.underlying=y.underlying))
        }

    } else if(length(class(y))==2 && any(class(y)=="ordered") && any(class(y)=="factor")) {
        y.unorderedfac = y
        class(y.unorderedfac) = "factor"

        retlist = prep.depvar.in(y.unorderedfac)
        retlist = append(retlist,list(is.ordered.factor=TRUE))
        return(retlist)

    } else {
        stop("The dependent variable must be either of class numeric, integer, character, or factor.")
    }

}


prep.depvar.out = function(y.prepped, y.underlying, factor.levels, is.ordered.factor) {
    # Undos what prep.depvar.in would do. This will be a helper function to any
    # predict functions.

    # if original had only 0s and 1s
    if(missing(y.underlying)) return(y.prepped)

    # if original had numerics other than just 0s and 1s, and strings, and factors
    y = vector(mode=class(y.underlying), length(y.prepped))
    y[y.prepped==0] = y.underlying[1]
    y[y.prepped==1] = y.underlying[2]

    if(missing(factor.levels)) return(y)

    # if original was factors
    y = factor(y, levels=seq(length(factor.levels)), labels=factor.levels)

    if(missing(is.ordered.factor)) return(y)

    # if original was ordered factors
    class(y) = c("ordered","factor")
    if(is.ordered.factor){
        return(y)
    } else {
        stop("Unexpected error")
    }

}
