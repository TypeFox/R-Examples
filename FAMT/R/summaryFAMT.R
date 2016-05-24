summaryFAMT <-
function (obj, pi0 = NULL, alpha = 0.15, info = c("ID", "Name")) 
{
    if (!any(is.element(c("FAMTdata", "FAMTmodel"), class(obj)[1]))) 
        stop("Class of obj should be either FAMTdata or FAMTmodel")
    if (class(obj)[1] == "FAMTdata") {
        data = obj
        m = nrow(data$expression)
        n = ncol(data$expression)
        nbnarow = (1:m)[apply(is.na(data$expression), 1, sum) > 
            0]
        nbnacol = (1:n)[apply(is.na(data$expression), 2, sum) > 
            0]
        s.expr = vector(length = 2, "list")
        names(s.expr) = c("Number of tests", "Sample size")
        s.expr[[1]] = m
        s.expr[[2]] = n
        res = list(expression = s.expr, covariates = summary(data$covariates), 
            annotations = summary(data$annotations))
    }
    if (class(obj)[1] == "FAMTmodel") {
        model = obj
        nbreject = matrix(0, nrow = length(alpha), ncol = 3)
        dimnames(nbreject) = list(1:length(alpha), c("alpha", 
            "Raw analysis", "FA analysis"))
        nbreject[, 1] = alpha
        if (is.null(pi0)) 
            pi0 = pi0FAMT(model, method = "smoother", diagnostic.plot = FALSE)
        m = length(model$adjpval)
        ordadj = order(model$adjpval)
        ordadjpval = model$adjpval[ordadj]
        fafdr = m * pi0 * ordadjpval/(1:m)
        nbreject[, 3] = unlist(lapply(alpha, function(t, fdr, 
            ordpval) {
            bool = fdr <= t
            if (all(bool) == FALSE) 
                res = 0
            if (any(bool) == TRUE) {
                kmax = max((1:length(fdr))[fdr <= t])
                threshold = ordpval[kmax]
                reject = (1:length(fdr))[ordpval <= threshold]
                res = length(reject)
            }
            return(res)
        }, fdr = fafdr, ordpval = ordadjpval))
        ord = order(model$pval)
        ordpval = model$pval[ord]
        rawfdr = m * pi0 * ordpval/(1:m)
        nbreject[, 2] = unlist(lapply(alpha, function(t, fdr, 
            ordpval) {
            bool = fdr <= t
            if (all(bool) == FALSE) 
                res = 0
            if (any(bool) == TRUE) {
                kmax = max((1:length(fdr))[fdr <= t])
                threshold = ordpval[kmax]
                reject = (1:length(fdr))[ordpval <= threshold]
                res = length(reject)
            }
            return(res)
        }, fdr = rawfdr, ordpval = ordpval))
        nbde = nbreject[nrow(nbreject), 3]
        if (nbde == 0) 
            DE = integer(0)
        if (nbde > 0) {
            de = order(model$adjpval)[1:nbde]
            columns = is.element(colnames(model$adjdata$annotations), 
                info)
            DE = model$adjdata$annotations[de, columns]
        }
        res = list(nbreject = nbreject, DE = DE, pi0 = pi0)
    }
    return(res)
}
