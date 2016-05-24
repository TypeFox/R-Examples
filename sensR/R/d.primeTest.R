##################################################################
### Exported functions:
##################################################################
dprime_compare <-
  function(correct, total, protocol, conf.level=0.95,
           statistic = c("likelihood", "Pearson", "Wald.p", "Wald.d"),
           estim = c("ML", "weighted.avg"))
### Test the hypothesis of any-diff among d-primes
{
  ## Match arguments:
  stat <- match.arg(statistic)
  estim <- match.arg(estim)
  stopifnot(is.numeric(conf.level) && length(conf.level) == 1 &&
            0 < conf.level && conf.level < 1)
  ## get data:
  .data <- dprime_table(correct, total, protocol)
  ## get dExp, se(dExp), ci
  dExp.list <- dprime_estim(data=.data, estim=estim)
  res1 <- dprime_testStat(data=.data, res=dExp.list,
                          conf.level=conf.level,  dprime0=0,
                          statistic="likelihood")
  ## get test statistic:
  res2 <- dprime_compareStat(data=.data, dExp=dExp.list$dExp,
                             statistic=stat)
  ## Return
  res <- c(res2, list(data=.data, coefficients=coef(res1),
                      conf.level=res1$conf.level,
                      conf.int=res1$conf.int, estim=estim,
                      conf.method=res1$conf.method))
  class(res) <- "dprime_compare"
  res
}

dprime_test <-
  function(correct, total, protocol, conf.level = 0.95,
           dprime0 = 0,
           statistic = c("likelihood", "Wald"),
           alternative = c("difference", "similarity", "two.sided",
             "less", "greater"),
           estim = c("ML", "weighted.avg"))
### Estimate dExp with std.err, CI, and
### Test H_0: common.dp == dprime0
{
  ## Match arg:
  stat <- match.arg(statistic)
  alt <- match.arg(alternative)
  estim <- match.arg(estim)
  ### Test arg:
  stopifnot(is.numeric(conf.level) && length(conf.level) == 1 &&
            0 < conf.level && conf.level < 1)
  stopifnot(is.numeric(dprime0) && length(dprime0) == 1 &&
            dprime0 >= 0)
  z <- 1
  if(isTRUE(all.equal(dprime0, 0)) &&
     !alternative %in% c("difference", "greater"))
    stop("'alternative' has to be 'difference'/'greater' if 'dprime0' is 0")
  if(alt == "difference") alt <- "greater"
  if(alt == "similarity") alt <- "less"
  ## get dp.table
  .data <- dprime_table(correct, total, protocol)
  ## get dExp, se(dExp), ci
  dExp.list <- dprime_estim(data=.data, estim=estim)
  ## get test statistic
  res <-
    dprime_testStat(data=.data, res=dExp.list, conf.level=conf.level,
                    dprime0=dprime0, statistic=stat)
  res <- c(res, dExp.list, list(data=.data, alternative=alt))
  ## get p-value:
  res$p.value <- normalPvalue(res$stat.value, alternative=alt)
  ## Return
  class(res) <- "dprime_test"
  res
}

posthoc <- function(x, ...) {
  UseMethod("posthoc")
}

posthoc.dprime_compare <-
  function(x, alpha = .05,
           test = c("pairwise", "common", "base", "zero"),
           base = 1,
           alternative = c("two.sided", "less", "greater"),
           statistic = c("likelihood", "Wald"), # "score", "exact"
           padj.method = c("holm", "bonferroni", "none"), ...)
{
  ## Match arguments:
  stat <- match.arg(statistic)
  alt <- match.arg(alternative)
  p.meth <- match.arg(padj.method[1], p.adjust.methods)
  ## Test arguments:
  stopifnot(round(base) == base)
  base <- as.integer(round(base))
  ## 'test' may be numeric or character:
  dprime0 <- 0 ## default value
  if(is.character(test)) {
    test <- match.arg(test)
  } else if(is.numeric(test)) {
    dprime0 <- test
    test <- "value"
    stopifnot(length(dprime0) == 1 && dprime0 >= 0)
  } else stop("'test' has to be numeric or character")
  N <- NROW(x$data)
  if(test == "base") dprime0 <- x$data[base, "dhat"]
  ## Get posthoc coef-table and test statistics:
  fit <- getPosthoc(data=x$data, base=base, dprime0=dprime0,
                    type=test, statistic=stat)
  ## Get p-values:
  coef <- fit$coef.table
  pval <- normalPvalue(fit$statistic, alternative=alt)
  coef[["p-value"]] <- p.adjust(pval, method=p.meth)
  ## Return:
  res <- c(x, list(posthoc=coef, test=test, ## alternative=alt,
                   padj.method=p.meth, base=base,
                   posthoc.stat=statistic))
  res$alternative <- alt
  ## Add letter displays to return list:
  if(test == "pairwise" && alt == "two.sided") {
    signifs <- setNames(coef[["p-value"]] < alpha, rownames(coef))
    ## multcomp:::insert_absorb changed since version 1.2-12 (now 1.2-17)
    ## ia.res <- multcomp:::insert_absorb(signifs, decreasing=TRUE)
    ia.res <- insert_absorb(signifs, decreasing=TRUE)
    letterOrder <- order(as.numeric(substring(names(ia.res$Letters), 6)))
    res$Letters <- ia.res$Letters[letterOrder]
  }
### FIXME: Add "protocol" to Letter display?
  if(!test %in% c("pairwise", "common"))
    res$dprime0 <- dprime0
  class(res) <- c(paste0("posthoc.", class(x)), class(x))
  res
}

posthoc.dprime_test <- posthoc.dprime_compare

##################################################################
### Auxiliary functions:
##################################################################
dprime_table <-
  function(correct, total, protocol, restrict.above.guess = TRUE)
### Return a data.frame with columns:
### correct, total, protocol, pHat, se.pHat, dprime, se.dprime
### restrict.above.guess: if pHat < pGuess then phat := pGuess
{
  ## Test arguments:
  ## Testing argument-mode:
  stopifnot(is.numeric(correct) && is.numeric(total))
  stopifnot(all(is.finite(c(correct, total))))
  if(is.factor(protocol))
    protocol <- as.character(protocol)
  if(!is.character(protocol))
    stop("'protocol' should be a factor or a character vector")
  stopifnot(is.logical(restrict.above.guess))
  restrict.above.guess <- restrict.above.guess[1]
  ## Testing argument conformity:
  stopifnot(all(correct == round(correct)))
  correct <- as.integer(round(correct))
  stopifnot(all(total == round(total)))
  total <- as.integer(round(total))
  stopifnot(length(correct) == length(total))
  stopifnot(length(correct) == length(protocol))
  stopifnot(length(correct) >= 1)
  stopifnot(all(correct <= total) && all(correct >= 0) &&
            all(total > 0))
  protValid <- c("triangle", "duotrio", "threeAFC", "twoAFC", "tetrad")
  stopifnot(all(protocol %in% protValid))
  ## Extract data:
  .data <- data.frame("correct"=correct, "total"=total,
                      "protocol"=protocol, stringsAsFactors=FALSE)
  rownames(.data) <- paste("group", 1:NROW(.data), sep="")
  x <- correct
  n <- total
  ## Compute pHat and se.pHat:
  pHat <- correct/total
  if(restrict.above.guess) {
    pG <- ifelse(protocol %in% c("duotrio", "twoAFC"), 1/2, 1/3)
    ## if pHat < pGuess then phat := pGuess:
    OK <- pHat > pG
    pHat[!OK] <- pG[!OK]
  }
  .data$pHat <- pHat
  .data$se.pHat <- se.pHat <- sqrt(pHat * (1 - pHat) / total)
  ## Compute d-primes and se(d-primes):
  .data$dprime <- sapply(1:NROW(.data), function(i)
                       psyinv(pHat[i], protocol[i]))
  se.dprime <- sapply(1:NROW(.data), function(i) {
    se.pHat[i] / psyderiv(.data$dprime[i], protocol[i]) })
  se.dprime[!is.finite(se.dprime)] <- NA
  .data$se.dprime <- se.dprime
  ## Return:
  .data
}

dprime_estim <-
  function(data, estim = c("ML", "weighted.avg"))
### Estimate dExp and se.dExp
{
  ## Match arguments:
  estim <- match.arg(estim)
  ## Return list:
  res <- list(estim=estim)
  ## Estimate common d-prime and std.err:
  if(estim == "ML") {
    fit <- nlminb(start=1, objective=dprime_nll, df=data, lower=0)
    res$dExp <- fit$par
    res$se.dExp <- ## se(dExp) not available if dExp <~ 1e-4
      if(fit$par < 1e-4) NA
      else
        as.vector(sqrt(1/hessian(func=dprime_nll, x=fit$par, df=data)))
    res$nll.dExp <- fit$objective
  }
  if(estim == "weighted.avg") {
    dprime <- data$dprime
    se.dprime <- data$se.dprime
    if(!all(is.finite(c(dprime, se.dprime))))
      stop("Boundary cases occured: use 'estim = ML' instead")
    w <- se.dprime^2
    w.prime <- w / sum(w)
    res$dExp <- sum(dprime / w) / sum(1 / w)
    ## sum(w.prime * dprime)
    res$se.dExp <- sqrt(sum(w.prime^2 * w))
  }
  ## Return
  res
}

dprime_nll <- function(dp, df) {
  ## nll under the null hypothesis of a single common d-prime, dp.
  meth <- as.character(df[, 3])
  nll.i <- sapply(1:nrow(df), function(i) {
    - dbinom(x=df[i, 1], size=df[i, 2],
             prob=psyfun(dp, method=meth[i]), log=TRUE) } )
  sum(nll.i)
}

dprime_testStat <-
  function(data, res, conf.level=0.95, dprime0=0,
           statistic = c("likelihood", "Wald")) ## "score",
### Test statistic and CI for test of common d-prime = dprime0.
###
### res: a list with dExp, se.dExp and optionally nll.dExp
{
  ## Match arg:
  stat <- match.arg(statistic)
  dExp <- res$dExp
  se.dExp <- res$se.dExp
  ## Get test statistics and conf.int:
  if(stat == "likelihood") {
    ## message(paste("likelihood CI not yet implemented;",
    ##               "reporting Wald interval instead"))
    nll.0 <- dprime_nll(dprime0, data)
    if(is.null(nll.dExp <- res$nll.dExp))
      nll.dExp <- dprime_nll(dExp, data)
    LR <- 2 * (nll.0 - nll.dExp)
    ## signed likelihood root statistic:
    statis <- sign(dExp - dprime0) * sqrt(abs(LR))
  }
  if(stat == "Wald") {
    if(!all(is.finite(c(dExp, se.dExp))))
      stop("Boundary cases occured: use 'statistic = likelihood' instead",
           call.=FALSE)
    statis <- (dExp - dprime0) / se.dExp
  }
  ## Compute Wald CI:
  a <- (1 - conf.level)/2
  ci <- dExp + c(1, -1) * qnorm(a) * se.dExp
  ci <- delimit(ci, 0, Inf)
  .coef <- data.frame(dExp, se.dExp, ci[1], ci[2])
  colnames(.coef) <- c("Estimate", "Std. Error", "Lower", "Upper")
  rownames(.coef) <- "d-prime"
  ## Return:
  list(coefficients=.coef, stat.value=statis, conf.int=ci,
       statistic=stat, dprime0=dprime0, conf.level=conf.level,
       conf.method="Wald")
}

dprime_compareStat <-
  function(data, dExp,
           statistic = c("likelihood", "Pearson", "Wald.p", "Wald.d"))
### Get test statistics for the any.diff hypotheses
###
### data: data.frame with
###     correct, total, protocol, pHat, se.pHat, dprime, se.dprime
### dExp: value of the common d.prime
{
  stat <- match.arg(statistic)
  x <- data$correct
  n <- data$total
  O <- c(x, n-x)
  prot <- as.character(data$protocol)
  ##
  if(stat == "likelihood") {
    pExp <- sapply(prot, function(x) psyfun(dExp, method = x))
    E <- n * c(pExp, 1 - pExp)
    X <- 2 * sum(O * log(O/E))
    ## Alternative computation of LR statistic that only works if the
    ## estimation method is also likelihood:
    ## nll.alt <- sum(apply(mat, 1, function(x) {
    ##   - dbinom(x=x[1], size=x[2], prob=x[3], log=TRUE) }))
    ## X <- -2 * (nll.alt - dprime_nll(dExp, df=data))
  }
  if(stat == "Pearson") {
    pExp <- sapply(prot, function(x) psyfun(dExp, method = x))
    E <- n * c(pExp, 1 - pExp)
    X <- sum( (O - E)^2 / E)
  }
  if(stat == "Wald.p") {
    pExp <- sapply(prot, function(x) psyfun(dExp, method = x))
    X <- sum( (data$pHat - pExp)^2/(data$pHat * (1 - data$pHat) / n))
  }
  if(stat == "Wald.d") {
    dprime <- data$dprime
    se.dprime <- data$se.dprime
    if(!all(is.finite(c(dprime, se.dprime))))
      stop("Boundary cases occured: use 'likelihood' or 'Pearson' instead")
    X <- sum( ((dprime - dExp)/se.dprime)^2 )
  }
  df <- NROW(data) - 1
  pval <-
    if(df >= 1) pchisq(X, df=NROW(data)-1, lower.tail=FALSE)
    else NA
  ## Return:
  list(stat.value=X, df=df, p.value=pval, statistic=stat)
}

getPosthoc <-
  function(data, base = 1, dprime0 = 0, #K = NULL,
           type = c("pairwise", "common", "base", "zero", "value"),
           statistic = c("likelihood", "Wald"), ...)
### coef.table: data.frame with columns:
###   1: d-prime, 2: se(d-prime), 3: protocol
### data: data.frame with columns:
###   1: successes, 2: total, 3: protocol
### dprime0: value of d-prime under the null hypothesis
### K: matrix defining Tukey (pairwise) or Dunnett's (base) contrasts
###
### Result: a list with 2 components:
###   - statistic: vector of test statistics
###   - coef.table: data.frame with d-prime, se(d-prime), optionally
###     CI, p-value
{
  ## Match arguments:
  type <- match.arg(type)
  stat <- match.arg(statistic)
  N <- NROW(data)
  nll.alt0 <- sapply(1:N, function(i)
                     dprime_nll(data$dprime[i], df=data[i, 1:3]))
  ## Compute coefficient table and test statistics:
  if(type %in% c("pairwise", "base")) {
    K <- switch(EXPR = type,
                "pairwise" = contrMat(setNames(1:N, rownames(data)),
                  type="Tukey"),
                "base" = contrMat(setNames(1:N, rownames(data)),
                  type="Dunnett", base=base))
    coef <- data.frame("Estimate" = K %*% data[, "dprime"])
    coef$"Std. Error" <- sqrt(abs(K) %*% data[, "se.dprime"]^2)
    ## Should we have CI here?
    if(stat == "likelihood") {
      ## list of pairwise group indicators:
      select.list <- lapply(1:nrow(K), function(i) which(K[i, ] != 0))
      ## Compute the null-likelihood:
      nll.0 <- sapply(select.list, function(sel) {
        nlminb(start=1, objective=dprime_nll, df=data[sel, ],
               lower=0)$objective })
      ## Compute the likelihood under the alternative:
      nll.alt <- sapply(select.list, function(sel) sum(nll.alt0[sel]))
      ## likelihood root statistics:
      LR <- -2 * (nll.alt - nll.0)
      statis <- sign(K %*% data[, "dprime"]) * sqrt(abs(LR))
    }
    else  ## statistic == "Wald":
      statis <- coef[, "Estimate"] / coef[, "Std. Error"]
  } ## end type %in% c("pairwise", "base")
  else { ## common, zero, <numerical value>
    coef <- data[, c("dprime", "se.dprime")]
    colnames(coef) <- c("Estimate", "Std. Error")
### NOTE: This is needed since discrim does not allow dprime0:
    pd0 <- sapply(1:N, function(i) {
      meth <- data$protocol[i]
      pg <- ifelse(meth %in% c("duotrio", "twoAFC"), .5, 1/3)
      pc <- psyfun(dprime0, method=meth)
      pc2pd(pc, pg)
    })
    discrim.list <- lapply(1:N, function(i) {
      discrim(data[i, 1], data[i, 2], method=data[i, 3], pd0=pd0[i],
              statistic=stat) })
    coef$"Lower" <- sapply(discrim.list, function(dl) coef(dl)[3, 3] )
    coef$"Upper" <- sapply(discrim.list, function(dl) coef(dl)[3, 4] )
    ## Get test statistic:
    if(type == "common") {
      ## Fit null model:
      fit2 <- nlminb(start=1, objective=dprime_nll, df=data, lower=0)
      ## Likelihood under the null:
      nll.0c <- fit2$objective
      ## Likelihood under the alternative:
      alt <- sapply(seq_len(nrow(data)), function(sel) {
        nlminb(start=1, objective=dprime_nll, df=data[-sel, ],
               lower=0)$objective })
      nll.altc <- alt + nll.alt0
      ## likelihood root statistics and p-values:
      LR <- -2 * (nll.altc - nll.0c)
      statis <- sign(data[, "dprime"] - dprime0) * sqrt(abs(LR))
    }
    else ## type == "zero" or <numerical value>
      statis <- sapply(discrim.list, "[[", "stat.value")
  }
  ## Return:
  list(statistic=statis, coef.table=coef)
}

##################################################################
### Print methods:
##################################################################
print.dprime_compare <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\n\tTest of multiple d-primes:\n\n")
  Estim <- if(x$estim == "ML") "Maximum likelihood" else
  "Weighted average"
  cat(paste("Estimation method: ",
            Estim, "\n", sep=""))
  cat(paste(format(x$conf.level, digits=2), "% ",
            "two-sided confidence interval method: ",
            x$conf.method, "\n", sep=""))
  cat("\nEstimate of common d-prime:\n")
  print(x$coefficients, quote = FALSE, digits = digits, ...)

  cat("\nSignificance test:\n")
  cat("  Null hypothesis:  All d-primes are equal\n")
  cat("  Alternative:      At least 2 d-primes are different\n")
  stat.chr <- if(x$statistic == "likelihood")
    "Likelihood Ratio" else x$statistic
  cat(paste("  Chi-square statistic (", stat.chr,") = ",
            format(x$stat.value, digits=digits), ", df = ", x$df,
            "\n  p-value = ", format.pval(x$p.value), "\n", sep=""))
  cat("\n")
  ## Return:
  invisible(x)
}

print.dprime_test <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  ## Print initial messages:
  cat("\n\tTest of common d-prime:\n\n")
  Estim <- if(x$estim == "ML") "Maximum likelihood" else
  "Weighted average"
  cat(paste("Estimation method: ",
            Estim, "\n", sep=""))
  cat(paste(format(x$conf.level, digits=2), "% ",
            "two-sided confidence interval method: ",
            x$conf.method, "\n", sep=""))
  ## Print estimate:
  cat("\nEstimate of common d-prime:\n")
  print(x$coefficients, quote = FALSE, digits = digits, ...)
  ## Print siginificance test:
  cat("\nSignificance test:\n")
  stat.chr <- if(x$statistic == "likelihood")
    "Likelihood root" else "Wald"
  alt.chr <- switch(x$alternative,
                    "two.sided" = "different from",
                    "greater" = "greater than",
                    "less" = "less than",
                    "similarity" = "less than",
                    "difference" = "greater than")
  cat(paste("  ", stat.chr, " statistic = ",
            format(x$stat.value, digits=digits),
            ", p-value: ", format.pval(x$p.value, digits=digits),
            "\n", sep=""))
  cat(paste("  Alternative hypothesis: d-prime is", alt.chr,
            format(x$dprime0, digits=digits), "\n"))
  cat("\n")
  ## Return:
  invisible(x)
}

print.posthoc.dprime_compare <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\n\tPost-hoc comparison of d-primes:\n\n")

  if(x$test == "pairwise")
    cat("Pairwise d-prime differences:\n")
  else if(x$test == "base")
    cat(paste("Differences to group ", x$base, ":\n", sep=""))
  else
    cat("Group-wise d-primes:\n")
  ph.mat <- x$posthoc
  ph.mat["p-value"] <- format.pval(ph.mat$"p-value", digits=digits)
  print(ph.mat, quote = FALSE, digits = digits, ...)
  cat("---\n")
  if(x$padj.method == "none")
    cat("p-values are not adjusted for multiplicity\n")
  else
    cat(paste("p-values are adjusted with ", x$padj.method,
              "'s method\n", sep=""))
  cat("\n")
  ## cat("p-values pertain to the alternative hypotheses:\n")
  cat("Alternative hypotheses:\n")
  alt.mess <- switch(x$alternative,
                     "two.sided" = "different from",
                     "greater" = "greater than",
                     "less" = "less than")
  if(x$test == "pairwise")
    cat(paste("  pairwise differences are", alt.mess, "zero\n"))
  if(x$test == "common")
    cat(paste("  d-primes are", alt.mess, "common d-prime\n"))
  if(x$test == "base")
    cat(paste("  d-primes differences are", alt.mess, "zero\n"))
  if(x$test %in% c("zero", "value"))
    cat(paste("  d-primes are", alt.mess,
              format(x$dprime0, digits=digits), "\n"))
### FIXME print message about the confidence interval method?
### FIXME print method about test statistic (Wald, likelihood, ...)?
  if(x$test == "pairwise" && !is.null(x$Letters)) {
    cat("\nLetter display based on pairwise comparisons:\n  ")
    print(x$Letters)
  }
  cat("\n")
  ## Return:
  invisible(x)
}

print.posthoc.dprime_test <- print.posthoc.dprime_compare

##################################################################
### Functions to produce compact letter display adopted from the
### multcomp package version 1.2-12
##################################################################

insert_absorb <- function( x, Letters=c(letters, LETTERS), separator=".", decreasing = decreasing ){

  obj_x <- deparse(substitute(x))
  namx <- names(x)
  namx <- gsub(" ", "", names(x))
  if(length(namx) != length(x))
    stop("Names required for ", obj_x)
  split_names <- strsplit(namx, "-")
  stopifnot( sapply(split_names, length) == 2 )
  comps <- t(as.matrix(as.data.frame(split_names)))
  rownames(comps) <- names(x)
  lvls <- unique(as.vector(comps))
  n <- length(lvls)
  lmat <- array(TRUE, dim=c(n,1), dimnames=list(lvls, NULL) )

  if( sum(x) == 0 ){                                                        # no differences
    ltrs <- rep(get_letters(1, Letters=Letters, separator=separator), length(lvls) )
    names(ltrs) <- lvls
    colnames(lmat) <- ltrs[1]
    msl <- ltrs
    ret <- list(Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat)
    class(ret) <- "multcompLetters"
    return(ret)
  }
  else{
    signifs <- comps[x,,drop=FALSE]

    absorb <- function(m){
      for(j in 1:(ncol(m)-1)){
        for(k in (j+1):ncol(m)){
          if( all(m[which(m[,k]),k] & m[which(m[,k]),j]) ){                 # column k fully contained in column j
            m <- m[,-k, drop=FALSE]
            return(absorb(m))
          }
          else if( all(m[which(m[,j]),k] & m[which(m[,j]),j]) ){            # column j fully contained in column k
            m <- m[,-j, drop=FALSE]
            return(absorb(m))
          }
        }
      }
      return(m)
    }
    for( i in 1:nrow(signifs) ){                                            # insert
      tmpcomp <- signifs[i,]
      wassert <- which(lmat[tmpcomp[1],] & lmat[tmpcomp[2],])               # which columns wrongly assert nonsignificance
      if(any(wassert)){
        tmpcols <- lmat[,wassert,drop=FALSE]
        tmpcols[tmpcomp[2],] <- FALSE
        lmat[tmpcomp[1],wassert] <- FALSE
        lmat <- cbind(lmat, tmpcols)
        colnames(lmat) <- get_letters( ncol(lmat), Letters=Letters,
                                       separator=separator)
        if(ncol(lmat) > 1){                                                 # absorb columns if possible
          lmat <- absorb(lmat)
          colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                         separator=separator )
        }
      }
    }
  }
  lmat <- lmat[,order(apply(lmat, 2, sum))]
  lmat <- sweepLetters(lmat)                                                                  # 1st
  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x))))))]                # reorder columns
  colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                 separator=separator)
  lmat <- lmat[,order(apply(lmat, 2, sum))]                                                   # 2nd sweep
  lmat <- sweepLetters(lmat)
  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x)))),
                           decreasing = decreasing))]                # reorder columns
  colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                 separator=separator)
  ltrs <- apply(lmat,1,function(x) return(paste(names(x)[which(x)], sep="", collapse="") ) )
  msl <- matrix(ncol=ncol(lmat), nrow=nrow(lmat))                                             # prepare monospaced letters
  for( i in 1:nrow(lmat) ){
    msl[i,which(lmat[i,])] <- colnames(lmat)[which(lmat[i,])]
    absent <- which(!lmat[i,])
    if( length(absent) < 2 ){
      if( length(absent) == 0 )
        next
      else{
        msl[i,absent] <- paste( rep(" ", nchar(colnames(lmat)[absent])), collapse="" )
      }
    }
    else{
      msl[i,absent] <- unlist( lapply( sapply( nchar(colnames(lmat)[absent]),
                                               function(x) return(rep( " ",x)) ),
                                       paste, collapse="") )
    }
  }
  msl <- apply(msl, 1, paste, collapse="")
  names(msl) <- rownames(lmat)
  ret <- list( Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat,
               aLetters = Letters, aseparator = separator )
  class(ret) <- "multcompLetters"
  return(ret)
}


# All redundant letters are swept out without altering the information within a LetterMatrix.
#
# mat         ... a LetterMatrix as produced by function insert_absorb()
# start.col   ... either a single integer specifying the column to start with or a vector
#                 of max. length equal to ncol(mat) specifying the column order to be used.
# Letters     ... a set of user defined letters { default is Letters=c(letters, LETTERS) }
# separator   ... a separating character used to produce a sufficiently large set of
#                 characters for a compact letter display (default is separator=".") in case
#                 the number of letters required exceeds the number of letters available

sweepLetters <- function(mat, start.col=1, Letters=c(letters, LETTERS), separator="."){

  stopifnot( all(start.col %in% 1:ncol(mat)) )
  locked <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))          # 1 indicates that another letter dependes on this entry
  cols <- 1:ncol(mat)
  cols <- cols[c( start.col, cols[-start.col] )]
  if( any(is.na(cols) ) )
    cols <- cols[-which(is.na(cols))]

  for( i in cols){
    tmp <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))
    tmp[which(mat[,i]),] <- mat[which(mat[,i]),]                        # get items of those rows which are TRUE in col "i"
    one <- which(tmp[,i]==1)

    if( all(apply(tmp[,-i,drop=FALSE], 1, function(x) return( any(x==1) ))) ){     # there is at least one row "l" where mat[l,i] is the only item which is TRUE i.e. no item can be removed in this column
      next
    }
    for( j in one ){                                                    # over all 1's
      if( locked[j,i] == 1 ){                                           # item is locked
        next
      }
      chck <- 0
      lck <- list()
      for( k in one ){
        if( j==k ){
          next
        }
        else{                                                           # pair j-k
          rows <- tmp[c(j,k),]
          dbl <- rows[1,] & rows[2,]
          hit <- which(dbl)
          hit <- hit[-which(hit==i)]
          dbl <- rows[1,-i,drop=FALSE] & rows[2,-i,drop=FALSE]
          if( any(dbl) ){
            chck <- chck + 1
            lck[[chck]] <- list(c(j,hit[length(hit)]), c(k,hit[length(hit)]))      # record items which have to be locked, use last column if multiple hits
          }
        }
      }
      if( (chck == (length(one)-1)) && chck != 0 ){                     # item is redundant
        for( k in 1:length(lck) ){                                      # lock items
          locked[ lck[[k]][[1]][1], lck[[k]][[1]][2] ] <- 1
          locked[ lck[[k]][[2]][1], lck[[k]][[2]][2] ] <- 1
        }
        mat[j,i] <- FALSE                                               # delete redundant entry
      }
    }
    if(all(mat[,i]==FALSE)){                                           # delete column where each entry is FALSE and restart
      mat <- mat[,-i,drop=FALSE]
      colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
      return(sweepLetters(mat, Letters=Letters, separator=separator))
    }
  }
  onlyF <- apply(mat, 2, function(x) return(all(!x)))
  if( any(onlyF) ){                                                     # There are columns with just FALSE entries
    mat <- mat[,-which(onlyF),drop=FALSE]
    colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
  }
  return( mat )
}

# Create a set of letters for a letter display. If "n" exceeds the number of letters
# specified in "Letters", they are recycled with one or more separating character(s)
# preceding each recycled letter.
# e.g. get_letters(10, Letters=letters[1:4]) produces:  "a"   "b"   "c"   "d"   ".a"  ".b"  ".c"  ".d"  "..a" "..b"
#
# n             ... number of letters
# Letters       ... the set of characters to be used
# separator     ... a character to be used as separator e.g.
#                   n=5, Letters=c("a","b") => "a", "b", ".a", ".b", "..a"

get_letters <- function( n, Letters=c(letters, LETTERS), separator="." ){

  n.complete <- floor(n / length(Letters))        # number of complete sets of Letters
  n.partial <- n %% length(Letters)               # number of additional Letters
  lett <- character()
  separ=""
  if( n.complete > 0 ){
    for( i in 1:n.complete ){
      lett <- c(lett, paste(separ, Letters, sep="") )
      separ <- paste( separ, separator, sep="" )
    }
  }
  if(n.partial > 0 )
    lett <- c(lett, paste(separ, Letters[1:n.partial], sep="") )
  return(lett)
}


##################################################################
### THE END
##################################################################
