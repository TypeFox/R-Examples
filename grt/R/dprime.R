dprime <- function(x,
    category,
    response,
    par = list(),
    zlimit = Inf,
    type = c("SampleIdeal", "Observer"))
{
    type <- match.arg(type)
    if(type == "SampleIdeal"){
        if(nrow(x) != length(category))
            stop("nrow(x) and length(category) are not the same")
        g <- as.factor(category)
        lev <- levels(g)
        if(length(lev) > 2) stop("more than 2 categories found")
        nc <- as.vector(table(g))
        if(any(sapply(par[c("noise","coeffs","bias")],is.null))){
            mc <- mcovs.default(x = x, grouping = category, pooled = TRUE)
            dprime <- dprimef(mc$means, mc$covs, noise=par$noise)
        } else {
            if(par$noise < 0) stop("par$noise must be positive non-zero value")
            h_coefs <- c(par$coeffs, par$bias)/par$noise
            h <- -as.matrix(cbind(x,1)) %*% as.matrix(h_coefs)
            # Truncate the large z-scores
            h[abs(h) > zlimit] <- sign(h[abs(h) > zlimit]) * zlimit
            rates <- tapply(h, g, function(x){sum(pnorm(x))}) / nc
            tmp <- qnorm(1 - rates)
            dprime <- as.vector(tmp[2] - tmp[1])
        }
    }else if(type == "Observer"){
        if(length(category) != length(response))
            stop("length(category) and length(response) are not the same")
        rates <- table(response,category)
        tmp <- qnorm(1-(rates[1,]/colSums(rates)))
        dprime <- as.vector(tmp[2] - tmp[1])
    }
    dprime
}

dprimef <- function(means, covs, noise=NULL)
{
    if(!is.list(means)) stop("means is not a list")
    if(is.list(covs)) covs <- as.matrix(covs[[1L]])
    if(!is.null(noise)){
        if(noise < 0) stop("noise must be positive non-zero value")
        covs <- covs + diag((noise^2), nrow=nrow(covs), ncol=ncol(covs))
    }
    tmp <- (means[[1]] - means[[2]]) %*% qr.solve(covs) %*% (means[[1]] - means[[2]])
    dprime <- as.vector(sqrt(tmp))
    dprime
}