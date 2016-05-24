.onLoad <- function(libname, pkgname) {
    # allow file suffixes of .Rnw, .Snw and .nw
    if (R.version$major >=3)
        tools::vignetteEngine("noweb", weave=noweave, tangle=notangle,
                              pattern="[.][RS]?nw", package="noweb")

    # The "data" directory is a more logical place for this, but when
    #  noweb is building itself, this is needed earlier in the process.
    assign("nowebSyntax",
           list(
             doc =      "^@",
             code=     "^<<(.*)>>=.*",
             coderef = "^ *<<([^>]*)>> *",
             docexpr = "\\\\Sexpr\\{([^\\}]*)\\}",
             sqexpr  = "\\[\\[([^\\]]*)]]",
             escapeat= "@<<[^>]*>>",
             verbatim= c("verbatim", "semiverbatim", "Verbatim"),
             verb = "\\\\[Vv]erb"),
           envir = asNamespace("noweb"))
}
