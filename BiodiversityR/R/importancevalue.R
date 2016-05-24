`importancevalue` <-
function(x, site="plotID", species="species", 
    count="count", basal="basal", 
    factor="forest", level="") 
{

    if (any(names(x) == site) == F) {stop ("site variable not defined")}
    if (any(names(x) == species) == F) {stop ("species variable not defined")}
    if (any(names(x) == count) == F) {stop ("count variable not defined")}
    if (any(names(x) == basal) == F) {stop ("basal area variable not defined")}
    if (factor != "") {
        if (any(levels(droplevels(factor(x[, factor]))) == level) == F) {stop ("specified level not among factor levels")}
        subs <- x[, factor]==level
        x <- x[subs, ,drop=F]
    }
    species.names <- levels(droplevels(factor(x[, species])))
    p <- length(species.names)
    result <- array(dim=c(p, 7))
    colnames(result) <- c("frequency", "density", "dominance", "frequency.percent", "density.percent", "dominance.percent", "importance.value")
    rownames(result) <- species.names
    total.plots <- sum(table(x[, site]) > 0)
    for (j in 1:p) {
        subs <- x[, species] == rownames(result)[j]
        spec.data <- x[subs, , drop=F]
        spec.data <- data.frame(spec.data)
        result[j, "frequency"] <- sum(table(spec.data[, site]) > 0) / total.plots
        result[j, "density"] <- sum(spec.data[, count])
        result[j, "dominance"]  <- sum(spec.data[, basal])             
    }
    total.freq <- sum(result[, "frequency"])
    total.density <- sum(x[, count])
    total.dominance <- sum(x[, basal])
    for (j in 1:p) {
        result[j, "frequency.percent"] <- result[j, "frequency"] / total.freq * 100
        result[j, "density.percent"] <- result[j, "density"] / total.density * 100
        result[j, "dominance.percent"] <- result[j, "dominance"] / total.dominance * 100
        result[j, "importance.value"] <- sum(result[j, c("frequency.percent", "density.percent", "dominance.percent")])
    }    
    result <- result[order(result[, "importance.value"], decreasing=T),]
    return(result)
}


`importancevalue.comp` <-
function(x, site="plotID", species="species", 
    count="count", basal="basal", 
    factor="forest") 
{

    groups <- table(x[, factor])
    m <- length(groups)
    levels <- names(groups)
    result <- list(values=levels)
    for (i in 1:m) {
        resultx <- importancevalue(x=x, site=site, species=species, count=count, basal=basal, 
            factor=factor, level=levels[i])
        result[[levels[i]]] <- resultx
    }
    return(result)
}


