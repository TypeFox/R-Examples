calculate.CN <-
function(coverage, method.pooled="amplicon", method.mt="Bonferroni", thres.cov=100, thres.p=0.05) {
  method.pooled <- method.pooled[1]
  method.mt <- method.mt[1]
  m.coverage.amplicon <- apply(coverage, 1, mean)
  m.coverage.sample <- apply(coverage, 2, median)
  C <- scale(coverage, center=FALSE, scale=m.coverage.sample/mean(m.coverage.sample))
  m <- apply(C, 1, median)
  CN <- 2 * t(scale(t(C), center=FALSE, scale=m))
  nsamples <- ncol(CN)
  namplicons <- nrow(CN)
  index.ok <- which(m.coverage.amplicon >= thres.cov)
  S <- matrix(nrow=namplicons+1, ncol=21)  
  rownames(S) <- c(rownames(CN), "pooled")
  colnames(S) <- c("mean_cov", "W_CN", "p_CN", "mean_CN", "sd_CN", "median_CN", "mad_CN", "loss_lsig_CN", "loss_sig_CN", "loss_hsig_CN", "gain_lsig_CN", "gain_sig_CN", "gain_hsig_CN", "loss_sig_n", "loss_hsig_n", "gain_sig_n", "gain_hsig_n", "loss_sig_pooled_n", "loss_hsig_pooled_n", "gain_sig_pooled_n", "gain_hsig_pooled_n")
  S[, "mean_cov"] <- c(m.coverage.amplicon, mean(coverage[index.ok, ]))
  S[, "mean_CN"] <- c(apply(CN, 1, mean), mean(CN[index.ok, ]))
  S[, "sd_CN"] <- c(apply(CN, 1, sd), sd(CN[index.ok, ]))
  S[, "median_CN"] <- c(apply(CN, 1, median), median(CN[index.ok, ]))
  s <- c(apply(CN, 1, mad), mad(CN[index.ok, ]))
  S[, "mad_CN"] <- s
  S[, "gain_lsig_CN"] <- sapply(s, function(sigma){qnorm(thres.p, mean=2, sd=sigma, lower.tail=FALSE)})
  S[, "gain_sig_CN"] <- sapply(s, function(sigma){qnorm(thres.p/nsamples, mean=2, sd=sigma, lower.tail=FALSE)})
  S[, "gain_hsig_CN"] <- sapply(s, function(sigma){qnorm(thres.p/nsamples/namplicons, mean=2, sd=sigma, lower.tail=FALSE)})
  S[, "loss_lsig_CN"] <- sapply(s, function(sigma){qnorm(thres.p, mean=2, sd=sigma, lower.tail=TRUE)})
  S[, "loss_sig_CN"] <- sapply(s, function(sigma){qnorm(thres.p/nsamples, mean=2, sd=sigma, lower.tail=TRUE)})
  S[, "loss_hsig_CN"] <- sapply(s, function(sigma){qnorm(thres.p/nsamples/namplicons, mean=2, sd=sigma, lower.tail=TRUE)})
  for (i in 1:namplicons) {
    model <- shapiro.test(CN[i, ])
    S[i, "W_CN"] <- model$statistic
    S[i, "p_CN"] <- model$p.value
    cn <- CN[i, ]
    S[i, "gain_sig_n"] <- length(which(cn >  S[i, "gain_sig_CN"]))
    S[i, "gain_sig_pooled_n"] <- length(which(cn >  S["pooled", "gain_sig_CN"]))
    S[i, "gain_hsig_n"] <- length(which(cn >  S[i, "gain_hsig_CN"]))
    S[i, "gain_hsig_pooled_n"] <- length(which(cn >  S["pooled", "gain_hsig_CN"]))
    S[i, "loss_sig_n"] <- length(which(cn <  S[i, "loss_sig_CN"]))
    S[i, "loss_sig_pooled_n"] <- length(which(cn <  S["pooled", "loss_sig_CN"]))
    S[i, "loss_hsig_n"] <- length(which(cn <  S[i, "loss_hsig_CN"]))
    S[i, "loss_hsig_pooled_n"] <- length(which(cn <  S["pooled", "loss_hsig_CN"]))
  }
  amplicons <- rownames(CN)[index.ok]
  CN <- CN[amplicons, ]
  samples <- colnames(CN)
  namplicons <- length(amplicons)
  nsamples <- length(samples)
  P <- matrix(nrow=namplicons*nsamples, ncol=6)
  colnames(P) <- c("amplicon", "sample", "CN", "p", "p_samples", "p_samples_amplicons") 
  rownames(P) <- paste(rep(amplicons, each=nsamples), rep(samples, namplicons), sep="_")
  cat("Fitting null models for", length(amplicons), "amplicons (each dot represents 10 amplicons)") 
  i <- 0
  for (x in amplicons) {
    i <- i+1 
    if (i %% 10 == 0) cat(".")
    for (y in samples) {
      z <- paste(x, y, sep="_")
      P[z, "amplicon"] <- x
      P[z, "sample"] <- y
      cn <- CN[x, y]
      P[z, "CN"] <- cn
      if (method.pooled == "pooled") SD <- S["pooled", "mad_CN"]
      if (method.pooled == "amplicon") SD <- S[x, "mad_CN"]
      if (cn < 2) P[z, "p"] <- pnorm(cn, mean=2, sd=SD, lower.tail=TRUE)
      if (cn >= 2) P[z, "p"] <- pnorm(cn, mean=2, sd=SD, lower.tail=FALSE)
    }
    index <- which(P[, "amplicon"] == x)
    p <- as.numeric(P[index, "p"])
    M <- multtest::mt.rawp2adjp(p, proc=method.mt)
    P[index, "p_samples"] <- M$adjp[order(M$index), method.mt]
  }
  cat("\n")
  p <- as.numeric(P[, "p"])
  M <- multtest::mt.rawp2adjp(p, proc=method.mt)
  P[, "p_samples_amplicons"] <- M$adjp[order(M$index), method.mt]
  attr(CN, "model") <- S
  attr(CN, "P") <- P
  attr(CN, "method.pooled") <- method.pooled
  attr(CN, "method.mt")  <- method.mt
  attr(CN, "thres.p") <- thres.p
  attr(CN, "thres.cov") <- thres.cov
  attr(CN, "amplicons.removed") <- setdiff(rownames(coverage), amplicons) 
  return(CN)
}
