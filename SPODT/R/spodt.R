spodt <- function(formula ,data,
                  weight=FALSE, graft=0,
                  level.max=5, min.parent=10, min.child=5, rtwo.min=0.001)
{
    if (class(data)!="SpatialPointsDataFrame") stop("use a SpatialPointsDataFrame")
    if (is.na((is.projected(data)))|(! is.projected(data)) ) warning("the coordinates are not projected. Please, provide projected coordinates or be sure to use euclidian coordinates!")

    coord.x    <- coordinates(data)[,1]
    coord.y    <- coordinates(data)[,2]
    loc.data   <- as.numeric(row.names(data@data))
    data.temp  <- data@data

    Call <- match.call()
    indx <- match(c("formula", "data"), names(Call), nomatch = 0L)
    if (indx[1] == 0L) stop("a 'formula' with the cofactors is required\n for single spatial analysis (with no cofactor) the right hand side should be z~1")

    dataset.prep <- model.frame(Call$formula, data=data.temp)

    dataset <- cbind(loc.data, coord.x, coord.y, dataset.prep)
    colnames(dataset)[1:4] <- c("loc", "x", "y", "z")

    spodt.fct(dataset, weight=weight, graft=graft, level.max=level.max,
    min.parent=min.parent, min.child=min.child, rtwo.min=rtwo.min)
}