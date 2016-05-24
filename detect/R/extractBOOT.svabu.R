extractBOOT.svabu <-
function(object, model = c("full", "sta", "det", "phi", "disp"), ...)
{
    tmp <- attr(object, "bootstrap")
    if (is.null(tmp))
        return(NULL)
    model <- match.arg(model)
    if (!inherits(object, "svabu_nb") && model == "disp")
        stop("model = 'disp' available for svabu_nb models only.")
    if (!object$zeroinfl && model == "phi")
        stop("not ZI model, can't provide values for 'phi'")
    if (inherits(object, "svabu_nb")) {
        log_sigma <- tmp["log.sigma",]
        att <- attributes(tmp)
        att$dim <- att$dimnames <- NULL
        tmp <- tmp[rownames(tmp) != "log.sigma",]
        att$dim <- dim(tmp)
        att$dimnames <- dimnames(tmp)
        attributes(tmp) <- att
    }
    if (model == "disp") {
        cfs <- mean(log_sigma)
        ses <- sd(log_sigma)
        vcv <- matrix(ses^2, 1, 1)
        cors <- matrix(1,1,1)
        names(cfs) <- names(ses) <- "log.sigma"
        colnames(vcv) <- rownames(vcv) <- "log.sigma"
        colnames(cors) <- rownames(cors) <- "log.sigma"
    } else {
        cf <- coef(object, "full")
        cfs <- rowMeans(tmp)
        ses <- apply(tmp, 1, sd)
        vcv <- cov(t(tmp))
        cors <- cor(t(tmp))
        names(cfs) <- names(ses) <- names(cf)
        colnames(vcv) <- rownames(vcv) <- names(cf)
        colnames(cors) <- rownames(cors) <- names(cf)
        wi <- switch(model,
            "full"=rep(TRUE, length(cf)),
            "sta"=grepl("sta_", rownames(tmp)),
            "det"=grepl("det_", rownames(tmp)),
            "zif"=grepl("zif_", rownames(tmp)))
        cfs <- cfs[wi]
        ses <- ses[wi]
        vcv <- vcv[wi, wi, drop=FALSE]
        cors <- cors[wi, wi, drop=FALSE]
    }
    list(coefficients=cfs, std.error=ses, vcov=vcv, cor=cors,
        B = ncol(tmp) - 1, type = attr(tmp, "type"), model = model)
}

