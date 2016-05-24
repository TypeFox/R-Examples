"specModel" <-
function (theta, model) 
{
    if (model@specfun == "bspline") 
        E <- calcEbspline(unlist(theta), model)
    if (model@specfun == "gaus") 
        E <- calcEhiergaus(theta, model@x2, model@nupow)
    E
}
