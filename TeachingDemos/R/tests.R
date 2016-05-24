SnowsPenultimateNormalityTest <- function(x){


        # the following function works for current implementations of R
        # to my knowledge, eventually it may need to be expanded
        is.rational <- function(x){
                rep( TRUE, length(x) )
        }


        tmp.p <- if( any(is.rational(x))) {
                0
        } else {
                # current implementation will not get here if length
                # of x is positive.  This part is reserved for the
                # ultimate test
                1
        }


        out <- list(
                p.value = tmp.p,
                alternative = strwrap(paste('The data does not come from a',
        'strict normal distribution (but may represent a distribution',
        'that is close enough)'), prefix="\n\t"),
                method = "Snow's Penultimate Normality Test",
                data.name = deparse(substitute(x))
        )


        class(out) <- 'htest'
        out
}

SnowsCorrectlySizedButOtherwiseUselessTestOfAnything <- function(x,
            data.name=deparse(substitute(x)), alternative='You Are Lucky',
                                                                 ...,
                                                                 seed) {
    if( !missing(seed) ) {
        if( is.numeric(seed) ) {
            set.seed(seed)
        } else {
            char2seed(seed)
        }
    }

    tmp.p <- runif(1)

    out <- list(
                p.value = tmp.p,
                data.name=data.name,
                method = "Snow's Correctly Sized But Otherwise Useless Test of Anything",
                alternative=alternative)
    if( !missing(seed) ) out$seed <- seed
    names(tmp.p) <- 'Random Uniform Value'
    out$statistic <- tmp.p

    class(out) <- 'htest'
    return(out)
}

