pmp.bma <-
function (bmao, oldstyle = FALSE) 
{
    if (!(is.bma(bmao) || is.topmod(bmao))) 
        stop("bmao needs to be a 'bma' object!")
    if (is.topmod(bmao)) {
        topmods = bmao
        was.enum = FALSE
        cumsumweights = sum(topmods$ncount())
        log.null.lik = 0
    }
    else {
        topmods = bmao$topmod
        log.null.lik = (1 - bmao$info$N)/2 * log(as.vector(crossprod(bmao$X.data[, 
            1] - mean(bmao$X.data[, 1]))))
        cumsumweights = bmao$info$cumsumweights
        was.enum = (bmao$arguments$mcmc == "enum")
    }
    lt1 = suppressWarnings(topmods$lik() - max(topmods$lik()))
    lt1 = exp(lt1)/sum(exp(lt1))
    if (was.enum) {
        lt2 = exp(topmods$lik() - log.null.lik)/cumsumweights
    }
    else {
        lt2 = topmods$ncount()/cumsumweights
    }
    cpoint = min(length(lt1), length(lt2))
    lt1 = lt1[1:cpoint]
    lt2 = lt2[1:cpoint]
    if (!oldstyle) 
        lt1 <- lt1 * sum(lt2)
    topmodout = rbind(lt1, lt2)
    rownames(topmodout) = c("PMP (Exact)", "PMP (MCMC)")
    colnames(topmodout) = topmods$bool()
    return(t(topmodout))
}
