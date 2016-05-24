setClass("NBMinerParameter",
    representation(
        pi      = "numeric",
        theta   = "numeric",
        n       = "integer",
        k       = "numeric",
        a       = "numeric",
        minlen  = "integer",
        maxlen  = "integer",
        rules   = "logical"
    )
)

setClass("NBMinerControl",
    representation(
        verbose = "logical",
        debug = "logical"
    ),
    prototype(verbose = FALSE, debug = FALSE)
)


