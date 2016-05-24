################################################################################
################################################################################
## author Mathieu Cossutta <mcossutta@gmail.com>   
## author Davide Morselli <dmorsell@gunil.ch>   
## author Till Junge <till.junge@gmail.com>                                   ##
##                                                                            ##
## Copyright (c) UNIL (Universite de Lausanne)                                ##
## NCCR - LIVES (National Centre of Competence in Research LIVES              ##
## Overcoming vulnerability: life course perspectives",                       ##
## <http://www.lives-nccr.ch/>)                                               ##
##                                                                            ##
## geospacom and spacom are free software: you can redistribute it and/or modify it under    ##
## the terms of the GNU General Public License as published by the Free       ##
## Software Foundation, either version 2 of the License or any later version. ##
##                                                                            ##
## geospacom and spacom are distributed in the hope that they will be useful, but WITHOUT ANY  ##
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  ##
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      ##
## details, see <http://www.gnu.org/licenses/>.                               ##
################################################################################
################################################################################


fix.holes <-
  function(poly.dat){
    n.poly.all<-numeric()
    for (k in 1:nrow(poly.dat@data)){
      n.poly.all[k]<-length(poly.dat@polygons[[k]]@Polygons)
    }
    has.hole<-which(n.poly.all>1)
    n.poly<-n.poly.all[has.hole]
    for (k in 1:length(has.hole)){
      print(paste("Fixing holes in the shapefile. Please wait. Polygon", k,"of",length(has.hole)))
      poly.dat@polygons[[has.hole[k]]]@Polygons[[n.poly[k]]]@hole<-F
    }
    return(poly.dat)
  }




## Warning the shp file should contain in the data frame a field with the id
##  of the area. If not, you should add it by end or using the following
## function:
## so you should produce a correspondance vector v such that v[i] is the
##  label of the i area in the field v[i]


addID <- function(poly,correspondance,area,name){
  l <- nrow(poly@data)
  poly@data[[name]] <- rep(0,nrow(poly@data))
  for (i in 1:l){
    poly@data[[name]][poly@data[[area]]==correspondance[i]] <-i
  }
  return(poly)
}

performAddFields <-  function(poly,dataframe,context.id,names=NULL) {
    y <- as.character(names(dataframe))
    y <- y[!y==context.id]
    # Check that the names are in y
    for (name in names) {
      if (!name %in% y) {
        stop("The variable '", name, "' is not in the context data")}}
    if (is.null(names)){names <-y}
    l <-length(poly)
    m <-length(names)
    for (k in 1:m){
      name <- names[[k]]
      poly@data[[name]]<- rep(NA,l)
      
      for (i in 1:l)
      {
        #poly@data[[name]][poly@data[[context.id]]==i]<- dataframe[[name]][i]
        poly@data[[name]][poly@data[[context.id]]==dataframe[[context.id]][i]]<- dataframe[[name]][i]
      }
    }
    return(list(poly,names))
  }

performSinglePlot <- function(poly,name,main=name,method="equal",nbr=10,...)
  {
    at = classIntervals(unique(poly@data[[name]][!is.na(poly@data[[name]])]), n = nbr, style = method, ...)$brks
    at[length(at)]<- at[length(at)]+0.01
    #col = brewer.pal(nbr,"YlOrRd")
    col=rev(heat.colors(nbr))
    spplot(
      poly,
      name,
      at=at,
      col.regions=col,
      main=main,
      mar=c(0,0,0,0),...
    )
  }
