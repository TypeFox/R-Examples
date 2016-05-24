"anova.grouped" <-
function(object, object2, ...){
    if(!inherits(object, "grouped")) 
        stop("Use only with 'grouped' objects.\n")
    if(missing(object2))
        stop("'anova.grouped()' performs only an LRT between two nested models; for \n\t\t\t\ta single model use 'summary.grouped()'.\n")
    if(!inherits(object2, "grouped"))
        stop(deparse(substitute(object2)), " must be a 'grouped' object.\n")
    dets0 <- object$details
    dets1 <- object2$details
    y0 <- dets0$y
    y1 <- dets1$y
    if(!all(y0 == y1))
        stop("either the two 'grouped' objects represent models fitted in different responses \nor models fitted under different coarsening mechanisms.\n")
    if(dets0$link != dets1$link)
        stop("the two models assume different link functions.\n")
    if(dets0$distr != dets1$distr)
        stop("the two models assume different distribution for the unerlying true response variable.\n")
    df0 <- length(object$coef)
    df1 <- length(object2$coef)
    if(df0 >= df1)
        stop(deparse(substitute(object)), " should be nested in ", deparse(substitute(object2)))
    L0 <- dets0$logLik
    L1 <- dets1$logLik
    L01 <- - 2 * (L0 - L1)
    if(L01 < 0 || !all(colnames(dets0$X) %in% colnames(dets1$X)))
        warning("it appears that the two models are not nested.\n")
    p.val <- 1 - pchisq(L01, df1 - df0)
    out <- list(name0 = deparse(substitute(object)), L0 = L0, df0 = df0, AIC0 = AIC(object), BIC0 = AIC(object, k = log(dets0$n)),
                name1 = deparse(substitute(object2)), L1 = L1, df1 = df1, AIC1 = AIC(object2), BIC1 = AIC(object, k = log(dets1$n)),
                L01 = L01, p.value = p.val)
    class(out) <- "aov.grouped"
    out
}

