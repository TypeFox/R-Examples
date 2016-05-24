`display.tests` <-
function (sig.rows = "all", summ = "anova", tests, form = parameter.list$form,
        use.model = parameter.list$use.model, ...){
    if (missing(tests)) {
        if (identical(sig.rows, "all")) {
            sig.rows <- 1:dim(sigs)[1]
        }
        tests <- match(rownames(sigs)[sig.rows], rownames(clust.mat))
        if (any(is.na(tests))) {
            tests <- tests[!is.na(tests)]
            warning("Nonexistent row requested from clust.mat")
        }
    }
    if(class(use.model)=="character"){
        use.model <- get(use.model)
    }
    ret <- lapply(tests, function(x) {
        dat <- data.frame(Y=t(clust.mat[x,,drop=FALSE]), unique(parameter.list$covariates))
        colnames(dat)[1] <- "Y"
        use.model(form, dat, ...)
    })
    if(!identical(summ,"none")){
        if(class(summ)=="character"){
            summ <- get(summ)
        }
        ret <- lapply(ret, summ)
    }
    names(ret) <- rownames(clust.mat)[tests]
    ret
}
