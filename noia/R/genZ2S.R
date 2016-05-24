genZ2S <-
function (genZ = NULL, reference = "F2", nloc = NULL, max.level = NULL, 
    max.dom = NULL) 
{
    "strrev" <- function(ss) {
        sapply(lapply(strsplit(ss, character(0)), rev), paste, 
            collapse = "")
    }
    if (is.null(nloc)) {
        if (is.null(genZ)) {
            stop("Function Z2S: number of loci unknown; either zmat or nloc must be provided")
        }
        nloc <- ncol(genZ)/3
    }
    ans <- 1
    for (i in 1:nloc) {
        eff <- colnames(ans)
        ans <- kronecker(Sloc(reference = reference, i, genZ), 
            ans)
        if (is.null(eff)) {
            colnames(ans) <- effectsNamesGeneral(1)
        }
        else {
            colnames(ans) <- strrev(kronecker(effectsNamesGeneral(1), 
                strrev(eff), "paste", sep = ""))
        }
        if (!(is.null(max.level) && is.null(max.dom))) 
            ans <- ans[, effectsSelect(nloc = nloc, max.level = max.level, 
                max.dom = max.dom, effects = colnames(ans))]
    }
    rownames(ans) <- genNames(nloc)
    return(ans)
}
