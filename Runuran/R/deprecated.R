#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   DEPRECATED functions!                                                 ##
##                                                                         ##
#############################################################################

#############################################################################
##                                                                          #
## Special sampling methods                                                 #
##                                                                          #
#############################################################################

## -- DAU: Alias-Urn Method ------------------------------------------------
## DEPRECATED!
## use function 'ur(dau.new(...),n)' instead

urdau <- function (n, probvector, from = 0, by = 1) {
        ## create UNU.RAN object
        unr <- dau.new(pv=probvector, from=0)
        ## draw sample
        if (from==0 && by==1)
                unuran.sample(unr,n)
        else
                from + by * unuran.sample(unr,n)
}

## -- DGT: Guide Table Method -----------------------------------------------
## DEPRECATED!
## use function 'ur(dgt.new(...),n)' instead

urdgt <- function (n, probvector, from = 0, by = 1) {
        ## create UNU.RAN object
        unr <- dgt.new(pv=probvector, from=0)
        ## draw sample
        if (from==0 && by==1)
                unuran.sample(unr,n)
        else
                from + by * unuran.sample(unr,n)
}

## -- End -------------------------------------------------------------------
