power.ladesign <- function(gsize, odds.ratio, sig.level = 0.05, statistic = c("Kruskal-Wallis", "Jonckheere-Terpstra"), alternative = c("two.sided", "one.sided"), nrep=1e+6) {
  ng <- length(gsize)
  nOR <- length(odds.ratio)
  if (nOR != ng-1) {
    errmesg <- paste(" There are", ng, "groups and", nOR,"odds ratios, there should be", ng-1)
    stop(message=errmesg)
  }
  statistic <- match.arg(statistic)
  kw <- statistic=="Kruskal-Wallis"
  alternative <- match.arg(alternative)
  if (kw & alternative == "one.sided") stop("Only two-sided alternative is allowed for Kruskal-Wallis test")
  ng <- length(gsize)
  zz <- .Fortran("lehman",
                 as.integer(ng),
                 as.integer(gsize),
                 double(ng),
                 as.double(rep(1, ng)),
                 as.double(sum(gsize)),
                 double(ng),
                 as.logical(kw),
                 as.integer(nrep),
                 tstat=double(nrep))
  if (kw) {
    ldq <- quantile(zz$tstat, 1-sig.level)
  } else {
    if (alternative=="one.sided") {
      ldq <- quantile(zz$tstat, c(sig.level, 1-sig.level))
    } else {
      ldq <- quantile(zz$tstat, c(sig.level/2, 1-sig.level/2))
    }
  }
  zz <- .Fortran("lehman",
                 as.integer(ng),
                 as.integer(gsize),
                 double(ng),
                 as.double(c(1,odds.ratio)),
                 as.double(sum(gsize*c(1,odds.ratio))),
                 double(ng),
                 as.logical(kw),                 
                 as.integer(nrep),
                 tstat=double(nrep))
  if (kw) {
    ldpow <- sum(zz$tstat >= ldq)/nrep
  } else {
    if (alternative=="one.sided") {
      ldpow <- max(sum(zz$tstat <= ldq[1]), sum(zz$tstat >= ldq[2]))/nrep
    } else {
      ldpow <- (sum(zz$tstat <= ldq[1]) + sum(zz$tstat >= ldq[2]))/nrep
    }
  }
  out <- list()
  out$group.size <- gsize
  out$odds.ratio <- odds.ratio
  out$statistic <- statistic
  out$alternative <- alternative
  out$sig.level <- sig.level
  out$power <- ldpow
  class(out) <- "ladesign"
  out
}

print.ladesign <- function(x, ...) {
  cat("         Number of groups =", length(x$group.size), "\n")
  cat("               group size =", x$group.size, "\n")
  cat("odds ratios w.r.t group 1 =", round(x$odds.ratio, 3), "\n")
  cat("           test statistic =", x$statistic, "\n")
  cat("              alternative =", x$alternative, "\n")
  cat("       significance level =", x$sig.level, "\n")
  cat("                    power =", round(x$power, 3), "\n")
  invisible(x)
}
