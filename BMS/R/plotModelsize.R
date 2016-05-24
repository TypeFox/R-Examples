plotModelsize <-
function (bmao, exact = FALSE, ksubset = NULL, include.legend = TRUE, 
    do.grid = TRUE, ...) 
{
    dotargs = match.call(expand.dots = FALSE)$...
    if (length(exact) > 1) {
        topmodidx = exact
        exact = TRUE
    }
    else {
        topmodidx = NA
    }
    K = bmao$info$K
    if (is.element("mprior.info", names(bmao))) 
        m = bmao$mprior.info$mp.msize
    else m = bmao$arguments$prior.msize
    pmp.10 = pmp.bma(bmao$topmod[topmodidx], oldstyle = TRUE)
    if (exact) {
        modelSmean = sum(apply(.post.topmod.bma(bmao$topmod[topmodidx]), 
            2, function(x) length(which(x == 1))) * pmp.10[, 
            1])
        modelS.var = sum(apply(.post.topmod.bma(bmao$topmod[topmodidx]), 
            2, function(x) length(which(x == 1)))^2 * pmp.10[, 
            1]) - modelSmean^2
        x = apply(.post.topmod.bma(bmao$topmod[topmodidx]), 2, 
            function(x) length(which(x == 1)))
        y = pmp.10[, 1]
        result = c()
        for (i in sort(unique(x))) result = c(result, sum(y[which(x == 
            i)]))
        names(result) = sort(unique(x))
        kvec = rep(0, (K + 1))
        kvec[(as.numeric(names(result)) + 1)] = result
    }
    else {
        k.vec = bmao$info$k.vec
        summi = sum(k.vec)
        modelSmean = sum((1:length(k.vec)) * (k.vec/summi)) - 
            1
        kvec = k.vec/sum(k.vec)
        modelSmean.sq = sum(((1:length(k.vec))^2) * (k.vec/summi))
        modelS.var = modelSmean.sq - modelSmean^2
    }
    upper = min(ceiling(modelSmean + 5 * modelS.var), K)
    lower = max(floor(modelSmean - 5 * modelS.var), 0)
    if (is.element("mp.Kdist", names(bmao$mprior.info))) {
        prior = bmao$mprior.info$mp.Kdist
    }
    else if (is.element("theta", names(bmao$arguments))) {
        theta = bmao$arguments$theta
        if (theta == "random") {
            beta.bin = function(a = 1, b = (K - m)/m, K = K, 
                w = 0:K) {
                return(lgamma(a + b) - (lgamma(a) + lgamma(b) + 
                  lgamma(a + b + K)) + log(choose(K, w)) + lgamma(a + 
                  w) + lgamma(b + K - w))
            }
            prior = exp(beta.bin(a = 1, b = (K - m)/m, K = K, 
                w = 0:K))
        }
        if (theta != "random") {
            prior = dbinom(x = 0:K, size = K, prob = m/K, log = FALSE)
        }
    }
    else {
        prior = rep(NA, length(kvec))
    }
    mat = cbind(kvec, prior)
    upper.ylim = max(kvec, prior, na.rm = TRUE)
    if (is.null(ksubset)) {
        ksubset = (lower:upper)
    }
    dotargs = .adjustdots(dotargs, type = "l", ylim = c(0, 1.1 * 
        upper.ylim), lwd = 1.5, xaxt = "n", col = c("steelblue3", 
        "tomato"), main = paste("Posterior Model Size Distribution", 
        "\n", "Mean:", round(modelSmean, 4)), cex.main = 0.8, 
        xlab = "Model Size", ylab = "", lty = 1:2, pch = 4, cex.axis = 0.9)
    matsubset = mat[ksubset + 1, ]
    eval(as.call(c(list(as.name("matplot"), as.name("matsubset")), 
        as.list(dotargs))))
    if (as.logical(do.grid)) 
        grid()
    points(kvec[ksubset + 1], cex = 0.8, pch = eval(dotargs$pch))
    axis(1, las = 1, at = 1:length(ksubset), labels = ksubset, 
        cex.axis = eval(dotargs$cex.axis))
    if (include.legend) {
        if (is.null(prior) || all(is.na(prior))) {
            legend(x = "topright", lty = eval(dotargs$lty), legend = c("Posterior"), 
                col = eval(dotargs$col), ncol = 1, bty = "n", 
                lwd = eval(dotargs$lwd))
        }
        else {
            legend(x = "topright", lty = eval(dotargs$lty), legend = c("Posterior", 
                "Prior"), col = eval(dotargs$col), ncol = 2, 
                bty = "n", lwd = eval(dotargs$lwd))
        }
    }
    return(invisible(list(mean = modelSmean, var = modelS.var, 
        dens = kvec)))
}
