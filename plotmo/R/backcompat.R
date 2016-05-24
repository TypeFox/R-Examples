# backcompat.R

# externs for plotmo 2.2.1
get.plotmo.clip.limits <- function(object, env, type, y, trace)
{
    warning0("get.plotmo.clip.limits is deprecated")
    NA
}
get.plotmo.default.type <- function(obj, env)
{
    warning0("get.plotmo.default.type is deprecated, please use plotmo_type instead")
    plotmo_type(object=obj, trace=0,
                fname="deprecated.get.plotmo.default.type", type=NULL)
}
get.plotmo.nresponse <- function(object, nresponse, y)
{
    warning0("get.plotmo.nresponse is deprecated, please use plotmo_nresponse instead")
    plotmo_nresponse(y=y, object=object, nresponse=nresponse,
                     trace=0, fname="deprecated.get.plotmo.nresponse")
}
get.plotmo.pairs <- function(object, env, x, trace, all2)
{
    warning0("get.plotmo.pairs is deprecated, please use plotmo.pairs instead")
    plotmo_pairs(object=object, x=x, nresponse=1, trace=trace,
                 all2=all2, degree2=TRUE)
}
get.plotmo.singles <- function(object, env, x, trace, all1)
{
    warning0("get.plotmo.singles is deprecated, please use plotmo.singles instead")
    plotmo_singles(object=object, x=x, nresponse=1, trace=trace,
        degree1=1, all1=all1)
}
get.plotmo.x <- function(object, env, trace)
{
    warning0("get.plotmo.x is deprecated, please use plotmo_x instead")
    plotmo_x(object=object, trace=trace, stringsAsFactors=TRUE)
}
get.plotmo.y <- function(object, env, y.column, expected.len, trace)
{
    warning0("get.plotmo.y is deprecated, please use plotmo_y instead")
    plotmo_y(object=object, nresponse=1, trace=trace)
}
get.plotmo.ylim <- function(object, env, type, trace)
{
    warning0("get.plotmo.ylim is deprecated")
    NA
}
# called by earth 4.2.0 TODO remove these when earth is updated
get.plotmo.type.wrapper <- function(object, env, type, func.name="plotmo")
    plotmo_type(object=object, trace=0, fname=func.name, type=type)

get.plotmo.x.wrapper <- function(object, env, trace=0)
    plotmo_x(object=object, trace=trace)

get.plotmo.y.wrapper <- function(object, env, y.column, expected.len, trace)
    plotmo_y(object=object, nresponse=y.column, trace=trace,
             expected.len=expected.len)$y

get.plotmo.y.default <- function(object, env, y.column, expected.len, trace)
    plotmo_y(object=object, nresponse=y.column, trace=trace,
             expected.len=expected.len)$y
