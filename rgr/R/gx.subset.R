gx.subset <-
function (dfname, subset = TRUE) 
{
    subset <- eval(substitute(subset), dfname)
    data.frame(lapply(dfname[subset, ], function(x) if (is.factor(x)) 
        x[, drop = TRUE]
    else x), row.names = row.names(dfname)[subset])
}
