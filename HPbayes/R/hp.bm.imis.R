hp.bm.imis <- function (prior, nrisk, ndeath, K, d = 10, B = 400, age = c(1e-05, 
    1, seq(5, 100, 5)), CI = 95) 
{
    low <- apply(prior, 2, min)
	high <- apply(prior, 2, min)   

opt.result <- loop.optim(prior = prior, nrisk = nrisk, ndeath = ndeath, d=d, theta.dim=ncol(prior), age=age)
    opt.mu.d <- opt.result$opt.mu.d
    opt.cov.d <- opt.result$opt.cov.d
    theta.new <- opt.result$theta.new
    d.keep <- opt.result$d.keep
    log.like.0 <- opt.result$log.like.0
    wts.0 <- opt.result$wts.0
    
samp.po <- samp.postopt(opt.cov.d = opt.cov.d, opt.mu.d = opt.mu.d, 
        prior = prior, d.keep = d.keep, B=B, d=d, B0=nrow(prior))
    H.k <- samp.po$H.k
    H.new <- samp.po$H.new
    B1 <- samp.po$B1

ll.postopt <- like.resamp(K = K, log.like.0 = log.like.0, 
opt.cov.d = opt.cov.d, opt.mu.d = opt.mu.d, d.keep = d.keep, d=d, theta.dim=ncol(prior))
    h.mu <- ll.postopt$h.mu
    h.sig <- ll.postopt$h.sig
    log.like <- ll.postopt$log.like
    K <- ll.postopt$K
    
result <- final.resamp(K = K, B1 = B1, B=B, H.new = H.new, H.k = H.k, 
        log.like = log.like, d.keep = d.keep, prior = prior, 
        h.mu = h.mu, h.sig = h.sig, nrisk = nrisk, ndeath = ndeath, 
        theta.dim=ncol(prior), age = age)
    H.final <- result$H.new
    nup <- result$nup
    vwts <- result$vwts
    ewts <- result$ewts
    frac.nup <- result$frac.nup
    mwts <- result$mwts
    wts.k <- result$wts.k
    mwt.case <- result$mwt.case
    
    loCI <- ((100 - CI)/2)/100
    hiCI <- 1 - (((100 - CI)/2)/100)
    loci <- round(rbind(quantile(H.final[, 1], loCI), quantile(H.final[, 
        2], loCI), quantile(H.final[, 3], loCI), quantile(H.final[, 
        4], loCI), quantile(H.final[, 5], loCI), quantile(H.final[, 
        6], loCI), quantile(H.final[, 7], loCI), quantile(H.final[, 
        8], loCI)), digits = 3)
    Median <- round(rbind(median(H.final[, 1]), median(H.final[, 
        2]), median(H.final[, 3]), median(H.final[, 4]), median(H.final[, 
        5]), median(H.final[, 6]), median(H.final[, 7]), median(H.final[, 
        8])), digits = 3)
    hici <- round(rbind(quantile(H.final[, 1], hiCI), quantile(H.final[, 
        2], hiCI), quantile(H.final[, 3], hiCI), quantile(H.final[, 
        4], hiCI), quantile(H.final[, 5], hiCI), quantile(H.final[, 
        6], hiCI), quantile(H.final[, 7], hiCI), quantile(H.final[, 
        8], hiCI)), digits = 3)
    pout <- data.frame(loci, Median, hici)
    names(pout) <- c("Low CI", "Median", "High CI")
    print(pout)
    return(list(H.final = H.final, h.mu = h.mu, h.sig = h.sig, 
        log.like = log.like, log.like.0 = log.like.0, wts.0 = wts.0, 
        d.keep = d.keep, nup = nup, frac.nup = frac.nup, vwts = vwts, 
        ewts = ewts, mwts = mwts, wts.k = wts.k, mwt.case = mwt.case, out = pout))
}