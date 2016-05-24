# missing values are not allowed.
`calculateMukv` <-
function (observations, predGrid, model, formulaString, ...)
{
    prG = predGrid
    obs = observations
    if (missing(formulaString) || is.null(formulaString)) {
         eq = dum ~ 1
    } else eq = formulaString
############ changed following 'if statement' 
    if (!"data.frame" %in% getSlots(class(obs)) & (terms(eq)[[3]] ==  
        1 || all(all.vars(eq)[-1] %in% dimnames(coordinates(obs))[[2]]))) {
      obs = SpatialPointsDataFrame(obs, data = data.frame(dum = rep(1,
            dim(coordinates(obs))[1])))
      names(obs) = as.character(eq[[2]])
    }
    red_ann_gam_krig = krige(eq, obs, prG, model)
    return(mean(red_ann_gam_krig$var1.var))
}

