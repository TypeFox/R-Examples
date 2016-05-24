## this is a hack and could be made more object oriented!
coeff.extract <- function (object) {
  if ("lm" %in% class(object)) { # includes glm
    return(coefficients(summary(object))[ , 1:2]) # always Estimate and Std. Error
  }
  if (class(object) == "coxph") {
    tmp <- data.frame(beta = object$coefficients, se = sqrt(diag(object$var)))
    rownames(tmp) <- names(object$coefficients)
    colnames(tmp) <- c("Estimate", "Std. Error")
    return(tmp)
  }
  stop("no way to extract coefficients from object of class :", class(object))
}
  
grs.onesnp.apply <- function (params, object, coeff.extract.fun = coeff.extract) {
  env = parent.frame() # to eval update in
  stopifnot(is.data.frame(params))
  stopifnot(all(c("snp", "coded.allele") %in% names(params)))
  params$codeB <- paste(params$snp, params$coded.allele, sep = "_")
  if ("beta" %in% names(params)) warning("overwriting existing beta column")
  params$beta <- NA
  if ("se" %in% names(params)) warning("overwriting existing se column")
  params$se <- NA
  for (idx in 1:nrow(params)) {
    ## update(..., evaluate = TRUE) uses eval(..., parent.frame()) where
    ## parent.frame() is the calling environment of update().  We need to go
    ## at least one frame further up the calling stack, hence:
    newcall <- update(object, formula = update.formula(formula(object), paste("~ . +", params$codeB[idx])), evaluate = FALSE)
    ## there is no obvious way to check that all params$codeB are accessible
    ## to update() without clumsily dissecting object$call, hence just catch:
    newcoef <- try(coeff.extract.fun(eval(newcall, env))) # silent = TRUE not helpful
    if (class(newcoef) != "try-error") {
      params$beta[idx] <- newcoef[params$codeB[idx], "Estimate"]
      params$se[idx] <- newcoef[params$codeB[idx], "Std. Error"]
    }
  }
  return(params)
}

grs.make.scores <- function (params, snpdata, appendage = ".score") {
  stopifnot(is.data.frame(params))
  stopifnot(all(c("snp", "coded.allele", "coef") %in% names(params)))
  params$codeB <- paste(params$snp, params$coded.allele, sep = "_")
  stopifnot(is.snpdata(snpdata))
  if(!all(params$codeB %in% names(snpdata$data))) snpdata <- align.snpdata.coding(params, snpdata)$snpdata
  stopifnot(all(params$codeB %in% names(snpdata$data))) # this should not happen
  for (score in unique(params$score)) {
    if (paste(score, appendage, sep = "") %in% names(snpdata$data)) warning(paste("overwriting column : ", score, appendage, sep = ""))
    snpdata$data[[paste(score, appendage, sep = "")]] <- apply(subset(snpdata$data, select = params$codeB[params$score == score]), 1, function(xx) return(sum(xx*params$coef[params$score == score])))
  }
  return(snpdata)
}

grs.summary <- function(w, b, s, n) {
  stopifnot(length(b) == length(w))
  stopifnot(length(s) == length(w))
  f <- !is.na(w) & !is.na(b) & !is.na(s) # filter
  m <- sum(f) # number of snps used
  X2m <- sum((b[f]*s[f]^-1)^2)
  R2m <- 1-exp(-X2m/n)
  ahat <- sum(w[f]*b[f]*s[f]^-2)/sum(w[f]^2*s[f]^-2)
  aSE <- sqrt(1/sum(w[f]^2*s[f]^-2))
  X2rs <- (ahat/aSE)^2
  pval <- pchisq(X2rs, df = 1, lower.tail = FALSE)
  R2rs <- 1-exp(-X2rs/n)
  Qrs <- X2m-X2rs
  phet <- if (m>1) pchisq(Qrs, df = m-1, lower.tail = FALSE) else NA
  return(list(m=m, n=n, X2m=X2m, R2m=R2m,
              ahat=ahat, aSE=aSE, X2rs=X2rs, R2rs=R2rs, pval=pval,
              Qrs=Qrs, phet=phet))
}

grs.plot <- function(w, b, s, text = NULL, textpos = NULL, textcex = 0.5, alpha = 0.05) {
  f <- !is.na(w) & !is.na(b) & !is.na(s)
  ws <- sign(w)
  plot((w/ws)[f], (b/ws)[f],
       xlim = c(0, 1.1*max((w/ws)[f])),
       ylim = range(c((b/ws)[f] + qnorm(.025)*s[f], (b/ws)[f] + qnorm(.975)*s[f])),
       type = "n", ann = FALSE, las = 1)
  abline(h = 0, lty = "dotted")
  ahat <- sum((w*b/s^2)[f])/sum((w^2/s^2)[f])
  aSE <- sqrt(1/sum((w^2*s^-2)[f]))
  abline(a = 0, b = ahat, col = "red")
  abline(a = 0, b = ahat + qnorm(alpha/2)*aSE, col = "red", lty = "dashed")
  abline(a = 0, b = ahat + qnorm(1-alpha/2)*aSE, col = "red", lty = "dashed")
  for (idx in which(f)) {
    lines(rep((w/ws)[idx], 2), (b/ws)[idx] + qnorm(c(alpha/2, 1-alpha/2))*s[idx], col = "grey")
  }
  if (!is.null(text)) {
    if (is.null(textpos)) textpos <- rep(1, length(text))
    text((w/ws)[f], (b/ws)[f], text[f], pos = textpos[f], cex = textcex)
  }
  points((w/ws)[f], (b/ws)[f])
}

grs.filter.Qrs <- function (w, b, s, p.thresh = 0.05) {
  stopifnot(length(b) == length(w))
  stopifnot(length(s) == length(w))
  f <- !is.na(w) & !is.na(b) & !is.na(s)

  while (sum(f) > 1) {
    if (pchisq(sum((b[f] * s[f]^-1)^2) - sum(w[f] * b[f] * s[f]^-2)^2/sum(w[f]^2 * s[f]^-2)^2 * sum(w[f]^2 * s[f]^-2), df = sum(f) - 1, lower.tail = FALSE) > p.thresh) break
    
    f[which.min(sapply(1:length(f), function(ii) {
      if (!f[ii]) return(NA)
      ff <- f & 1:length(f) != ii
      sum((b[ff] * s[ff]^-1)^2) - sum(w[ff] * b[ff] * s[ff]^-2)/sum(w[ff]^2 * s[ff]^-2)^2 * sum(w[ff]^2 * s[ff]^-2)}))] <- FALSE
  }
  return(f)
}


# hacky but works... check env of update calls inside grs.onesnp.apply is correct
test.subsample <- function(params, object, scorename, snpdata, nsub, permute = NULL) {
  subdata <- snpdata$data[sample(1:nrow(snpdata$data), nsub), ]
  if (!is.null(permute)) for (pp in permute) subdata[[pp]] <- sample(subdata[[pp]]) # permute these columns
  sm0 <- update(object, data = subdata)
  sm1 <- update(sm0, formula = update.formula(formula(sm0), paste("~ . +", scorename)))
  sstats <- grs.onesnp.apply(params, sm0, coeff.extract.fun = coeff.extract)
  ## hack to make this generic, length(residuals) gets back the sample size for most objects
  tmp <- grs.summary(sstats$coef, sstats$beta, sstats$se, length(sm0$residuals))
  tmp$ahat.exact <- coeff.extract(sm1)[scorename, "Estimate"]
  tmp$aSE.exact <- coeff.extract(sm1)[scorename, "Std. Error"]
  tmp$sign.concord <- sum(sign(sstats$coef)==sign(sstats$beta))
  tmp$sign.pval <- binom.test(tmp$sign.concord, nrow(sstats), alternative = "two.sided")$p.value
  return(unlist(tmp))
}

