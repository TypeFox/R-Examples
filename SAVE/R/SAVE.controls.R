##########################################################################
## SAVE main Method
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##                                                        and Rui Paulo
##
##########################################################################

SAVE.controls <-
function (lower = NULL, upper = NULL, optim.method = "BFGS", parinit = NULL,...) {
    #a<- match.call(expand.dots=T)
    a <- as.list(.expand.call(call=sys.call(sys.parent(0))))
    #print(a[-1])
    aa <- a[-1]
    frmls <- formals(deparse(a[[1]]))
    # To remove the '...'
    frmls <- frmls[-length(frmls)]
    par.deprecated <- which(!(names(aa) %in% names(frmls)))
    if (length(par.deprecated) != 0) {
        cat ("These parameters will be deprecated in the kriging at the stage I parameters (the ones invoved in the construction of the emulator) estimation step.\n")
        print (as.character(names(aa)[par.deprecated]))
        a <- a[-(par.deprecated+1)]
        cat ("These parameters will be passed to the kriging process at the stage I parameters (the ones invoved in the construction of the emulator) estimation step.\n")
        print (as.list(a[-1]))
    }
    #    kriging.controls <- list(lower = lower, upper = upper, optim.method = optim.method, parinit = parinit)
    kriging.controls <- a[-1]
    return(kriging.controls)
}
