tofKpredict.nohull <-
function (pars, fixedlist, localmaxrange)
{
    fK <- tofullKrigingspace(fittedlist = pars, fixedlist = fixedlist)
    return(purefn(fK, testhull = F, constraints = NULL))
}
