comps <-
function (x, bonds = c("entire", "strong", "weak")) 
{
    ifelse(isTRUE(is.null(dimnames(x)[1]) == TRUE | is.null(dimnames(x)[1][[1]]) == 
        TRUE) == TRUE, lbs <- 1:nrow(x), lbs <- dimnames(x)[[1]])
    bd <- bundles(x, lb2lb = FALSE, collapse = TRUE)
    switch(match.arg(bonds), entire = lbd <- bd, strong = lbd <- list(bd$recp, 
        bd$txch, bd$mixd, bd$full), weak = lbd <- list(bd$asym, 
        bd$tent))
    tx <- transl(unlist(lbd))
    while (length(tx) > length(transl(tx))) {
        tx <- transl(tx)
    }
    com <- list()
    for (i in 1:length(tx)) {
        com[[i]] <- lbs[as.numeric(dhc(tx[i]))]
    }
    rm(i)
    lt <- list(com = com, isol = lbs[which(!(dimnames(x)[[1]] %in% 
        dimnames(suppressWarnings(rm.isol(x, diag.incl = FALSE)))[[1]]))])
    return(lt)
}
