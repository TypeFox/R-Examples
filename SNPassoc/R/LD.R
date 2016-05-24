LD<-function(g1, ...) UseMethod("LD")

LD.snp<-
function (g1, g2, ...)
{
# adapted from LD.genotype in package genetics by Gregory Warnes et al
    if (!is.snp(g1) || !is.snp(g2))
        stop("Please supply snp objects")
    prop.A <- summary(g1)$allele.freq[, 2]/100
    prop.B <- summary(g2)$allele.freq[, 2]/100
    major.A <- names(prop.A)[which.max(prop.A)]
    major.B <- names(prop.B)[which.max(prop.B)]
    pA <- max(prop.A, na.rm = TRUE)
    pB <- max(prop.B, na.rm = TRUE)
    if (pA<1 & pB<1){
      pa <- 1 - pA
      pb <- 1 - pB
      Dmin <- max(-pA * pB, -pa * pb)
      pmin <- pA * pB + Dmin
      Dmax <- min(pA * pb, pB * pa)
      pmax <- pA * pB + Dmax

  #    counts <- table(allele.count(g1, major.A), allele.count(g2, major.B))

      counts <- table( as.numeric(reorder(g1, "common")), as.numeric(reorder(g2,"common")) )

      n3x3 <- matrix(0, nrow = 3, ncol = 3)
      colnames(n3x3) <- rownames(n3x3) <- 0:2
      for (i in rownames(counts)) for (j in colnames(counts)) n3x3[3 -
          as.numeric(i), 3 - as.numeric(j)] <- counts[i, j]
      loglik <- function(pAB, ...) {
          (2 * n3x3[1, 1] + n3x3[1, 2] + n3x3[2, 1]) * log(pAB) +
              (2 * n3x3[1, 3] + n3x3[1, 2] + n3x3[2, 3]) * log(pA -
                  pAB) + (2 * n3x3[3, 1] + n3x3[2, 1] + n3x3[3,
              2]) * log(pB - pAB) + (2 * n3x3[3, 3] + n3x3[3, 2] +
              n3x3[2, 3]) * log(1 - pA - pB + pAB) + n3x3[2, 2] *
              log(pAB * (1 - pA - pB + pAB) + (pA - pAB) * (pB -
                  pAB))
      }
      solution <- optimize(loglik, lower = pmin + .Machine$double.eps,
          upper = pmax - .Machine$double.eps, maximum = TRUE)
      pAB <- solution$maximum
      estD <- pAB - pA * pB
      if (estD > 0)
          estDp <- estD/Dmax
      else estDp <- estD/Dmin
      n <- sum(n3x3)
      corr <- estD/sqrt(pA * pB * pa * pb)
      dchi <- (2 * n * estD^2)/(pA * pa * pB * pb)
      dpval <- 1 - pchisq(dchi, 1)
      retval <- list(call = match.call(), D = estD, "D'" = estDp,
          r = corr, "R^2" = corr^2, n = n, "X^2" = dchi, "P-value" = dpval)
    } else
      retval <- list(call = match.call(), D = NA, "D'" = NA,
          r = NA, "R^2" = NA, n = sum(table(g1,g2)) , "X^2" = NA, "P-value" = NA)

    class(retval) <- "LD"
    retval
}

LD.setupSNP<-
function (g1, SNPs, ...)
{
# adapted from LD.data.frame in package genetics by Gregory Warnes et al
    if(missing(SNPs) & inherits(g1,"setupSNP")) SNPs<-labels(g1)
    g1 <- g1[, SNPs]
    
    P <- matrix(nrow = ncol(g1), ncol = ncol(g1))
    rownames(P) <- colnames(g1)
    colnames(P) <- colnames(g1)
    P <- D <- Dprime <- nobs <- chisq <- p.value <- corr <- R.2 <- P
    for (i in 1:(ncol(g1) - 1)) for (j in (i + 1):ncol(g1)) {
        ld <- LD(g1[, i], g1[, j])
        D[i, j] <- ld$D
        Dprime[i, j] <- ld$"D'"
        corr[i, j] <- ld$r
        R.2[i, j] <- ld$"R^2"
        nobs[i, j] <- ld$n
        chisq[i, j] <- ld$"X^2"
        p.value[i, j] <- ld$"P-value"
    }
    retval <- list(call = match.call(), D = D, "D'" = Dprime,
        r = corr, "R^2" = R.2, n = nobs, "X^2" = chisq, "P-value" = p.value)
    class(retval) <- "LD.data.frame"
    retval
}

LDtable<-
function (x, colorcut = c(0, 0.01, 0.025, 0.05, 0.1, 1), colors = heat.colors(length(colorcut)),
    textcol = "black", digits = 3, show.all = FALSE, which = c("D",
        "D'", "r", "X^2", "P-value", "n"), colorize = "P-value",
    cex, ...)
{
# adapted from LDtable in package genetics by Gregory Warnes et al
    if(!inherits(x,"LD.data.frame"))
       stop("Object should be LD.data.frame, output of LD")
    if (!colorize %in% names(x))
        stop(colorize, " not an element of ", deparse(substitute(x)))
    datatab <- summary(x)
    missmatch <- which[!(which %in% names(x))]
    if (length(missmatch) > 0)
        stop(missmatch, " not an element of ", deparse(substitute(x)))
    matform <- function(value, template) {
        dim(value) <- dim(template)
        dimnames(value) <- dimnames(template)
        value
    }
    tmp <- cut(x[[colorize]], colorcut, include.lowest = TRUE)
    colormat <- matform(as.numeric(tmp), x[[colorize]])
    n <- matform(paste("(", x$n, ")", sep = ""), x$n)
    if (!show.all) {
        colormat <- colormat[-nrow(colormat), -1, drop = FALSE]
        n <- n[-nrow(n), -1, drop = FALSE]
    }
    image(x = 1:ncol(colormat), y = 1:ncol(colormat), z = t(colormat[nrow(colormat):1,
        ]), col = colors, xlab = "Marker 2\n\n", ylab = "Marker 1",
        xaxt = "n", yaxt = "n", ...)
    abline(v = -0.5 + 1:(ncol(colormat) + 1))
    abline(h = -0.5 + 1:(nrow(colormat) + 1))
    axis(3, 1:ncol(colormat), colnames(colormat))
    axis(2, 1:nrow(colormat), rev(rownames(colormat)))
    cex.old <- par("cex")
    if (missing(cex))
        cex <- min(c(1/10, 1/(length(which) + 1))/c(strwidth("W"),
            strheight("W") * 1.5))
    par(cex = cex)
    lineheight <- strheight("W") * 1.5
    center <- lineheight * length(which)/2
    for (i in 1:length(which)) {
        displaymat <- x[[which[i]]]
        if (!show.all)
            displaymat <- displaymat[-nrow(displaymat), -1, drop = FALSE]
        if (which[i] == "P-value")
            displaymat <- format.pval(displaymat, digits = digits)
        else if (which[i] != "n")
            displaymat <- format(displaymat, digits = digits)
        displaymat[] <- gsub("NA.*", "", as.character(displaymat))
        text(x = col(colormat), y = nrow(colormat) - row(colormat) +
            1 + center - lineheight * (i - 1), displaymat, col = textcol,
            adj = c(0.5, 1))
    }
    text(x = 1, y = 1, paste(which, collapse = "\n"), adj = c(0.5,
        0.5))
    par(cex = cex.old)
    title(main = "Linkage Disequilibrium\n")
    invisible(colormat)
}


LDplot<-
function (x, digits = 3, marker, distance, which = c("D", "D'",
    "r", "X^2", "P-value", "n", " "), ...)
{
# adapted from LDplot in package genetics by Gregory Warnes et al
    if(!inherits(x,"LD.data.frame"))
       stop("Object should be LD.data.frame, output of LD")
    which = match.arg(which)
    if (missing(marker))
        marker <- colnames(x[[which]])
    else if (is.numeric(marker))
        marker <- colnames(x[[which]])[marker]
    datamat <- ifelse(is.na(x[[which]]), t(x[[which]]), x[[which]])
    if (which %in% c("D'", "r"))
        diag(datamat) <- 1
    else if (which == "P-value")
        diag(datamat) <- 0
    dimnames(datamat) <- dimnames(x[[which]])
    if (missing(distance))
        distance <- 1:ncol(datamat)
    distance <- matrix(distance, ncol = ncol(datamat), nrow = nrow(datamat),
        byrow = TRUE)
    dimnames(distance) <- dimnames(datamat)
    matplot(x = t(distance[marker, , drop = FALSE]), t(datamat[marker,
        , drop = FALSE]), type = "b", xlab = "Marker", ylab = paste("Linkage Disequilibrium: ",
        which, sep = ""), xaxt = "n", ...)
    axis(1, distance[1, ], paste(1:ncol(datamat), colnames(datamat),
        sep = ": "))
    title("Pairwise Disequilibrium Plot")
    invisible()
}

