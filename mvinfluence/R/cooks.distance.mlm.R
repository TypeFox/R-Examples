cooks.distance.mlm <-
function (model, infl = mlm.influence(model, do.coef = FALSE), ...) 
{
    cookd <- infl$CookD
    m <- infl$m
    names(cookd) <- if(m==1) infl$subsets else apply(infl$subsets,1, paste, collapse=',')
    cookd
}
