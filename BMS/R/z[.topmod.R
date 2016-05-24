`[.topmod` <-
function (x, i, ...) 
{
    tm = x
    idx = i
    if (any(is.na(suppressWarnings(as.integer(idx))))) 
        idx = 1:length(tm$lik())
    if (length(tm$betas_raw()) > 1) {
        bbeta = TRUE
        bet = as.vector(tm$betas()[, idx])
        bet = bet[bet != 0]
    }
    else {
        bbeta = FALSE
        bet = numeric(0)
    }
    if (length(tm$betas2_raw()) > 1) {
        bbeta2 = TRUE
        bet2 = as.vector(tm$betas2()[, idx])
        bet2 = bet2[bet2 != 0]
    }
    else {
        bbeta2 = FALSE
        bet2 = numeric(0)
    }
    fixvec = tm$fixed_vector()
    if (!length(as.vector(fixvec))) 
        fixvec = numeric(0)
    else fixvec = as.vector(t(fixvec[, idx]))
    .top10(nmaxregressors = tm$nregs, nbmodels = tm$nbmodels, 
        bbeta = bbeta, lengthfixedvec = nrow(tm$fixed_vector()), 
        bbeta2 = bbeta2, inivec_lik = tm$lik()[idx], inivec_bool = tm$bool()[idx], 
        inivec_count = tm$ncount()[idx], inivec_vbeta = bet, 
        inivec_vbeta2 = bet2, inivec_veck = tm$kvec_raw()[idx], 
        inivec_fixvec = fixvec)
}
