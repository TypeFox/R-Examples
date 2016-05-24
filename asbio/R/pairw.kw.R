pairw.kw <- function (y, x, conf = 0.95)
{
    comp <- !is.na(y)
    y <- y[comp]; x <- x[comp]
    ranks <- rank(y)
    fitted <- tapply(y, x, median)
    ni <- tapply(ranks, x, length)
    N <- length(ranks)
    r <- length(ni)

    z <- sort(match(ranks,ranks))
      t <- 0
      for (i in 1:(N - 1)) {
          tie.ranks <- sum(match(z, i),na.rm=T)
          if (tie.ranks > 1) t <- t + tie.ranks^3 - tie.ranks
          }

    mean.ranks <- tapply(ranks, x, mean)
    dif.mat <- outer(mean.ranks, mean.ranks, "-")
    diffs <- dif.mat[upper.tri(dif.mat)]
    
    ni.mat <- outer(1/ni, 1/ni, "+")
    ni.den <- ni.mat[upper.tri(ni.mat)]
    
    SE.diff <- sqrt(((N * (N + 1))/12 - t/(12 * (N - 1))) * ni.den)
    B <- qnorm(1 - ((1 - conf)/(r^2 - r)))
    p.val <- 2 * pnorm(abs(diffs)/SE.diff, lower.tail = FALSE)
    p.adj <- ifelse(p.val * ((r^2 - r)/2) >= 1, 1, round(p.val *
        ((r^2 - r)/2), 6))
    hwidths <- B * SE.diff
    val <- round(cbind(diffs, diffs - hwidths, diffs + hwidths),
        5)
    Decision <- ifelse((val[, 2] > 0 & val[, 3] > 0) | val[,
        2] < 0 & val[, 3] < 0, "Reject H0", "FTR H0")
    val <- as.data.frame(cbind(val, Decision, p.adj))
    lvl <- outer(levels(x), levels(x), function(x1, x2) {
        paste(paste("Avg.rank", x1, sep = ""), paste("Avg.rank",
            x2, sep = ""), sep = "-")
    })
    dimnames(val) <- list(lvl[upper.tri(lvl)], c("Diff", "Lower",
        "Upper", "Decision", "Adj. P-value"))
   head <- paste(paste(as.character(conf*100),"%",sep=""),c("Confidence intervals for Kruskal-Wallis comparisons"))
   ###
   res <- list()
   res$head <- head
   res$conf <- conf
   comp <- outer(levels(x),levels(x),function(x1,x2){paste(x1, x2, sep="-")})
   res$comp <- comp[upper.tri(comp)]
   res$summary <- val
   res$band <- cbind(diffs-hwidths, diffs+hwidths)
   res$fitted <- fitted
   res$x <- x
   res$y <- y
   res$method <- "kBonferroni"
   class(res)<-"pairw"
   res
    }