topmodels.bma <-
function (bmao) 
{
    if (!is.bma(bmao)) {
        stop("you need to provide a bma object")
    }
    return(.post.topmod.bma(bmao))
}
