"randtest.enfa" <- function(xtest, nrepet=999, ...)
{
    ## Verifications
    if (!inherits(xtest,"enfa"))
        stop("should be an object of class \"enfa\"")
    if (!isTRUE(all.equal(xtest$cw, rep(1,length(xtest$cw)))))
        warning("not yet implemented for unequal column weightsw: \n column weights not taken into account")

    ## External call to the C function "randenfar": randomizes the weights
    ## to test the significance of the first axis of specialization
    tab<-as.matrix(xtest$tab)
    pr<-xtest$pr
    res<-.C("randenfar", as.double(t(tab)), as.double(pr),
            as.integer(ncol(tab)), as.integer(nrow(tab)),
            as.integer(nrepet), double(nrepet), PACKAGE="adehabitatHS")[[6]]
    return(as.randtest(res, xtest$s[1], call = match.call()))
  }

