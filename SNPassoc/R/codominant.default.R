`codominant.default` <-
function(o)
{
 if (length(unique(o[!is.na(o)]))>3)
    stop("variable should have 3 levels max")
 else factor(o)
}

