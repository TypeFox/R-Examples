effectsSelect <-
function (nloc, max.level = NULL, max.dom = NULL, effects = NULL) 
{
    if (is.null(effects)) 
        return(effectsNamesGeneral(nloc, max.level, max.dom))
    effects <- effects[sapply(effects, statusMaxLevel, max.level)]
    effects <- effects[sapply(effects, statusMaxDom, max.dom)]
    return(effects)
}
