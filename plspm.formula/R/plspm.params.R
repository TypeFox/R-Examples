plspm.params <-
function (Formula, Data) 
{
    coupe <- strsplit(Formula, split = "\n")[[1]]
    coupe <- gsub(" ", "", coupe)
    ind.def <- grep("=~", coupe)
    if (length(ind.def) == 0) {
        stop("formules de relations externes non definies")
    }
    else {
    }
    ind.inter <- grep("~~", coupe)
    if (length(ind.inter) == 0) {
        stop("formules de relations structurelles non definies")
    }
    else {
    }
    res.def <- fraction.formula(coupe, ind.def)
    res.inter <- fraction.formula(coupe, ind.inter)
    ovect.mu <- omatrix.inner(res.def[[1]], res.inter, mat = TRUE, 
        iplot = FALSE)
    M <- ovect.mu$matrice
    ovect <- ovect.mu$ordre
    outer.l <- outer.list(res.def, res.inter, Data, ovect)
    result <- list(inner.mat = M, outer.list = outer.l)
    return(result)
}
