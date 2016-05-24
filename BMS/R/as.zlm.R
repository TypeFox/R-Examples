as.zlm <-
function (bmao, model = 1) 
{
    thiscall = match.call()
    if (!is.bma(bmao)) 
        stop("bmao needs to be a bma object")
    bools = bmao$topmod$bool()
    if (all(is.character(model)) && length(model) == 1) {
        model = (1:length(bools))[bools == model[[1]]]
        if (length(model) == 0) 
            stop("Provided model hex-index was not found in bmao object topmodels")
    }
    else if (all(is.character(model)) && (length(model) > 1)) {
        mix = match(model, bmao$reg.names)
        if (any(is.na(mix))) 
            stop("Provided variable names do not conform to bma object")
        ll = logical(bmao$info$K)
        ll[mix] = TRUE
        model = (1:length(bools))[bools == bin2hex(ll)]
        rm(ll, mix)
        if (length(model) == 0) 
            stop("Model conforming to provided variable names was not found in bmao object topmodels")
    }
    else if ((length(model) == bmao$info$K) && (is.numeric(model) || 
        is.logical(model))) {
        model = (1:length(bools))[bools == bin2hex(model)]
        if (length(model) == 0) 
            stop("Provided binary model index was not found in bmao object topmodels")
    }
    else if ((length(model) == 1) && (is.numeric(model) || is.logical(model))) {
        if (model < 1 | model > length(bools)) 
            stop("Provided numeric model index was not found in bmao object topmodels")
    }
    else stop("model needs to be an integer, logical or character model index representation (hexcode or variable names)")
    inclvbls = as.logical(bmao$topmod$bool_binary()[, model, 
        drop = TRUE])
    yXdf = as.data.frame(bmao$X.data)
    zlmres = zlm(as.formula(yXdf[, c(TRUE, inclvbls)]), data = yXdf, 
        g = bmao$gprior.info)
    zlmres$call <- thiscall
    return(zlmres)
}
