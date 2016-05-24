##' A function to bind the different entity level.
##'
##' A data frame is chosen over the list is solely for the purpose of
##' transition to ggplot2.
##'
##' @param territory The data frame which contains the
##' territory/country level data
##' @param subregion The sub aggregated region aggregate
##' @param region The macro region aggregate
##' @param world The world aggregate
##' @export

ebind = function(territory = NULL, subregion = NULL,
                    region = NULL, world = NULL){
    nTerritory = NROW(territory)
    nSubregion = NROW(subregion)
    nRegion = NROW(region)
    nWorld = NROW(world)
    tmp = rbind(territory, subregion, region, world)
    area.df = data.frame(tmp, Area = c(rep("Territory", nTerritory),
                                       rep("subRegion", nSubregion),
                                       rep("Region", nRegion),
                                       rep("World", nWorld)),
        stringsAsFactors = FALSE)
    area.df
}
