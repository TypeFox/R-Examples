mmm2 <-
function (formula, id, data = NULL, rtype = TRUE, interaction = NULL, 
    R = NULL, b = NULL, tol = 0.001, maxiter = 25, family = "gaussian", 
    corstr = "independence", Mv = 1, silent = TRUE, scale.fix = FALSE, 
    scale.value = 1) 
{
    mf <- model.frame(formula = formula, data = data)
    x <- as.matrix(model.matrix(attr(mf, "terms"), data = mf)[, -1])
    if (ncol(x) == 1) 
        colnames(x) <- colnames(model.matrix(attr(mf, "terms"), 
            data = mf))[2:(ncol(x) + 1)]
    colnames(x) <- gsub("\\s", "", gsub("asfactor", "", gsub("[[:punct:]]", "", colnames(x))))
    y <- model.response(mf)
    nresp <- ncol(y)
    nrowdata <- nrow(y)
    mresp <- matrix(as.numeric(t(y)))
    rep.nresp <- rep(nresp, nrow(x))
    if (ncol(x) == 1){
    covmat <- as.matrix(as.numeric(x)[rep(1:nrow(x), rep.nresp)])
    colnames(covmat) <- colnames(x)
    } else {
    covmat <- x[rep(1:nrow(x), rep.nresp), ]
    }
    if(is.vector(id) == TRUE){id <- as.matrix(id)
    } else {
    id <- as.matrix(model.frame(id, data))
    }
    id <- id[rep(1:nrow(id), rep.nresp), ]
    if (rtype == TRUE) {
        r <- matrix(rep(0, nresp * (nresp - 1)), nrow = nresp)
        for (i in 1:(nresp - 1)) {
            r[(nrow(r) - (i - 1)), i] <- 1
        }
        resptype <- NULL
        for (i in 1:(nrow(covmat)/nrow(r))) {
            resptype <- rbind(resptype, r)
        }
        if (ncol(resptype) == 1) {
            colnames(resptype) <- c("rtype")
        }
        else {
            colnames(resptype) <- paste("rtype", seq(1:(nresp - 
                1)), sep = "")
        }
        interact <- NULL
        if (length(interaction) != 0) {
            for (i in 1:length(interaction)) {
                interact <- cbind(interact, resptype * covmat[, 
                  interaction[i]])
            }
            covmat <- cbind(covmat, resptype, interact)
        }
        else {
            covmat <- cbind(covmat, resptype)
        }
    }
    covmat <- as.data.frame(covmat)
    if (length(interaction != 0)) {
        colnames(covmat)[(ncol(x) + nresp):ncol(covmat)] <- paste(colnames(resptype), 
            rep(colnames(x)[interaction], each = length(colnames(resptype))), 
            sep = "*")
    }
    covn1 <- colnames(covmat)[1]
    if (substr(colnames(x)[1], start = nchar(colnames(x)[1]), 
        stop = nchar(colnames(x)[1])) == "]") {
        vn1 <- paste("covariate", seq(ncol(x)), sep = "")
        vn2 <- paste(colnames(resptype), rep(vn1[interaction], 
            each = length(colnames(resptype))), sep = "*")
        colnames(covmat) <- c(vn1, colnames(resptype), vn2)
    }
    formula2 <- as.formula(paste("mresp ~ ", paste(colnames(covmat), collapse = "+")))
    #library(gee)
    fit <- gee(formula2, id = id, data = covmat, R = R, b = b, 
        tol = tol, maxiter = maxiter, family = family, corstr = corstr, 
        Mv = Mv, silent = silent, scale.fix = scale.fix, scale.value = scale.value)
    fit$title <- "Multivariate Marginal Models with Shared Regression Parameters"
    fit$version <- "Version 1.2 (12/2013)"
    fit
}
