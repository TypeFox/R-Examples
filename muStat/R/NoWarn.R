`NoWarn` <-
function(x){
    PrevWarn <- options(warn=(-1))$warn
    x
    options(warn=PrevWarn)
    invisible(endfunction())
}

