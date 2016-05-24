detrenderEnv <-
function ()

{

    pos <- match("detrenderEnv", search())

    if (is.na(pos)) {

        detrenderEnv <- list()

        attach(detrenderEnv, length(search())-1 )

        rm(detrenderEnv)

        pos <- match("detrenderEnv", search())

    }

    return(pos.to.env(pos))

}

