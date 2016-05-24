logbin.reparameterise <- function(np.coefs, terms, data, allref, design.ref)
{
    termlabels <- attr(terms, "term.labels")
    design.type <- sapply(allref, attr, "type")
    ref.vector <- as.vector(design.ref, mode = "integer")
    
    coefs.int <- coefs.int.reparam <- np.coefs[1]
    coefs.model <- coefs.model.reparam <- np.coefs[-1L]
    coef.count <- 0L
    
    coef.names <- "(Intercept)"
    
    for (i in seq_len(length(design.ref))) {
        varname <- gsub("`", "", termlabels[i])
        varref <- allref[[termlabels[i]]][[ref.vector[i]]]
        if (design.type[i] == 1) {
            thiscoef <- coefs.model[(coef.count + 1L)]
            cont.min <- min(data[[varname]])
            cont.max <- max(data[[varname]])
            if (varref == 1) {
                coefs.int.reparam <- coefs.int.reparam - thiscoef * cont.min
                coefs.model.reparam[(coef.count + 1L)] <- thiscoef
            } else {
                coefs.int.reparam <- coefs.int.reparam + thiscoef * cont.max
                coefs.model.reparam[(coef.count + 1L)] <- -thiscoef
            }
            coef.names <- append(coef.names, varname)
            coef.count <- coef.count + 1L
        } else if (design.type[i] == 2) {
            lev <- levels(factor(data[[varname]]))
            nlev <- nlevels(factor(data[[varname]]))
            ref.orig <- varref
            ref.new <- lev[1L]
            if (ref.orig != ref.new) {
                thiscoef <- coefs.model[(coef.count + 1L):(coef.count + nlev - 1L)]
                thiscoef.new <- append(thiscoef, 0, after = which(lev == ref.orig) - 1L)[-1L] - thiscoef[1L]
                coefs.int.reparam <- coefs.int.reparam + thiscoef[1L]
                coefs.model.reparam[(coef.count + 1L):(coef.count + nlev - 1L)] <- thiscoef.new
            }
            coef.names <- append(coef.names, paste(varname, lev[-1L], sep = ""))
            coef.count <- coef.count + nlev - 1L
        } else if (design.type[i] == 3) {
            lev <- levels(factor(data[[varname]]))
            nlev <- nlevels(factor(data[[varname]]))
            ref.orig <- varref
            thiscoef <- coefs.model[(coef.count + 1L):(coef.count + nlev - 1L)]
            thiscoef.sum <- append(rev(cumsum(rev(thiscoef))), 0)
            thiscoef.ord <- thiscoef.sum[match(lev, ref.orig)]
            thiscoef.new <- thiscoef.ord[-1L] - thiscoef.ord[1L]
            coefs.int.reparam <- coefs.int.reparam + thiscoef.ord[1L]
            coefs.model.reparam[(coef.count + 1L):(coef.count + nlev - 1L)] <- thiscoef.new
            coef.names <- append(coef.names, paste(varname, lev[-1L], sep = ""))
            coef.count <- coef.count + nlev - 1L
        }
    }
    coefs <- c(coefs.int.reparam, coefs.model.reparam)
    names(coefs) <- coef.names
    
    design <- model.matrix(terms, data)
    
    list(coefs = coefs, design = design)
}