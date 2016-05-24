genZ2X <-
function (genZ, reference = "F2", max.level = NULL, max.dom = NULL) 
{
    ans <- NULL
    eff <- effectsNamesGeneral(nloc = ncol(genZ)/3, max.level = max.level, 
        max.dom = max.dom)
    for (e in rev(eff)) {
        if (!(e %in% colnames(ans))) {
            partial <- partialX(genZ = genZ, reference = reference, 
                effect = e)
            partial <- partial[, !(colnames(partial) %in% colnames(ans)), 
                drop = FALSE]
            neweffects <- colnames(partial) %in% eff
            ans <- cbind(ans, partial[, neweffects, drop = FALSE])
        }
    }
    ans <- ans[, eff]
    rownames(ans) <- rownames(genZ)
    return(ans)
}
