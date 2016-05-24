
".phylobase.Options" <-
    list(retic = "fail",
         singleton = "warn",
         multiroot = "warn",
         poly = "ok",
         allow.duplicated.labels = "warn")


.onAttach <- function(library, pkg)
{
    ## we can't do this in .onLoad
    unlockBinding(".phylobase.Options", asNamespace("phylobase"))
}


