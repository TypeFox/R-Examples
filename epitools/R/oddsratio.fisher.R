oddsratio.fisher <-
function (x, y = NULL, conf.level = 0.95, rev = c("neither", 
    "rows", "columns", "both"), correction = FALSE, verbose = FALSE) 
{
    if (is.matrix(x) && !is.null(y)) {
        stop("y argument should be NULL")
    }
    if (is.null(y)) {
        x <- epitable(x, rev = rev)
    }
    else {
        x <- epitable(x, y, rev = rev)
    }
    tmx <- table.margins(x)
    p.exposed <- sweep(tmx, 2, tmx["Total", ], "/")
    p.outcome <- sweep(tmx, 1, tmx[, "Total"], "/")
    nr <- nrow(x)
    fisher <- matrix(NA, nr, 3)
    fisher[1, 1] <- 1
    for (i in 2:nr) {
        xx <- rbind(x[1, ], x[i, ])
        ftestxx <- fisher.test(xx, conf.level = conf.level)
        est <- ftestxx$estimate 
        ci <- ftestxx$conf.int 
        fisher[i, ] <- c(est, ci)
    }
    pv <- tab2by2.test(x, correction = correction)
    colnames(fisher) <- c("estimate", "lower", "upper")
    rownames(fisher) <- rownames(x)
    cn2 <- paste("odds ratio with", paste(100 * conf.level, "%", 
        sep = ""), "C.I.")
    names(dimnames(fisher)) <- c(names(dimnames(x))[1], cn2)
    rr <- list(x = x, data = tmx, p.exposed = p.exposed, p.outcome = p.outcome, 
        measure = fisher, conf.level = conf.level, p.value = pv$p.value, 
        correction = pv$correction)
    rrs <- list(data = tmx, measure = fisher, p.value = pv$p.value, 
        correction = pv$correction)
    attr(rr, "method") <- "Conditional MLE & exact CI from 'fisher.test'"
    attr(rrs, "method") <- "Conditional MLE & exact CI from 'fisher.test'"
    if (verbose == FALSE) {
        rrs
    }
    else rr
}
