setMethod("show", "MultivariateDistribution",
    function(object){
        txt <- gettextf("Distribution object of class: %s\n", class(object)[1])
        parameter <- param(object)
        Names <- slotNames(parameter)
        if(length(Names) > 1){
            for(i in Names[Names != "name"])
                txt <- c(txt,
                    gettextf("%s: %s\n", i, slot(parameter, i)))
        }
        return(txt)
    })
setMethod("show", "EuclCondition",
    function(object){
        cat(gettextf("name:\t%s\n", object@name))
        cat(gettextf("Range:\t%s with dimension %s\n", object@Range@name, object@Range@dimension))
    })
setMethod("show", "LMParameter",
    function(object){
        cat(gettextf("name:\t%s\n", object@name))
        cat(gettextf("theta:\t%s\n", object@theta))
        cat(gettextf("intercept:\t%s\n", object@intercept))
        cat(gettextf("scale:\t%s\n", object@scale))
    })
setMethod("show", "UnivariateCondDistribution",
    function(object){
        txt <- gettextf("Distribution object of class: %s\n", class(object)[1])
        parameter <- param(object)
        Names <- slotNames(parameter)
        if(length(Names) > 1){
            for(i in Names[Names != "name"])
                txt <- c(txt,
                    gettextf("%s: %s\n", i, slot(parameter, i)))
        }
        cat(txt)
        cat(gettext("## cond:\n"))
        show(object@cond)
    })
