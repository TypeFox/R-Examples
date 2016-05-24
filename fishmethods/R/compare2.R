compare2 <- function (readings, usecols = c(1, 2), twovsone = TRUE, plot.summary = TRUE, 
                         barplot = TRUE, chi = TRUE, pool.criterion = 1, cont.cor = TRUE, 
                         correct = "Yates", first.name = "Reader A", second.name = "Reader B") 
{
  if (length(usecols) < 2) 
    stop("Two column values required in usecols")
  if (any(is.na(usecols))) 
    stop("Two column values required in usecols")
  if (any(is.null(usecols))) 
    stop("Two column values required in usecols")
  if (ncol(readings) < 2) 
    stop("Only 1 column present in readings")
  if (!is.numeric(readings[, usecols[1]])) 
    stop("usecols[1] column is not numeric")
  if (!is.numeric(readings[, usecols[2]])) 
    stop("usecols[2] column is not numeric")
  thecall = match.call()
  if (plot.summary == TRUE | chi == TRUE) {
    above = sum(readings[, usecols[2]] > readings[, usecols[1]])
    below = sum(readings[, usecols[2]] < readings[, usecols[1]])
    on = sum(readings[, usecols[2]] == readings[, usecols[1]])
    a_summary = as.data.frame(c(above, on, below))
    name1 = paste(second.name, ">", first.name)
    name2 = paste(second.name, "=", first.name)
    name3 = paste(second.name, "<", first.name)
    rownames(a_summary) = c(name1, name2, name3)
  }
  if (twovsone == TRUE) {
    xx = table(readings[, c(usecols[1], usecols[2])])
    xx[xx == 0] = NA
    r1_ages = as.numeric(rownames(xx))
    r2_ages = as.numeric(colnames(xx))
    xmin = min(c(r1_ages, r2_ages))
    xmax = max(c(r1_ages, r2_ages))
    xlimits = c(xmin, xmax)
    opar = par(mar = c(5.1, 5.1, 4.1, 2.1), cex.lab = 1.3, 
               cex.axis = 1.2, yaxs = "r", las = 1)
    plot(r1_ages[1], r2_ages[1], xlim = xlimits, ylim = xlimits, 
         pch = as.character(xx[1, 1]), xlab = first.name, 
         ylab = second.name, axes = F)
    axis(1, at = seq(xmin, xmax, 2), labels = seq(xmin, xmax, 
                                                  2))
    axis(2, at = seq(xmin, xmax, 2), labels = seq(xmin, xmax, 
                                                  2), las = 1)
    box()
    for (i in 1:dim(xx)[1]) {
      for (j in 1:dim(xx)[2]) {
        text(r1_ages[i], r2_ages[j], labels = as.character(xx[i, 
                                                              j]), xlim = xlimits, ylim = xlimits, xlab = "", 
             ylab = "")
      }
    }
    points(xlimits, xlimits, type = "l", lty = 3)
    if (plot.summary == TRUE) {
      maxstr <- max(nchar(as.character(above)), nchar(as.character(on)), 
                    nchar(as.character(below)))
      a1 <- ifelse(nchar(as.character(above)) != maxstr, 
                   paste(rep(" ", maxstr - nchar(as.character(above))), 
                         above), above)
      o1 <- ifelse(nchar(as.character(on)) != maxstr, paste(rep(" ", 
                                                                maxstr - nchar(as.character(on))), on), on)
      b1 <- ifelse(nchar(as.character(below)) != maxstr, 
                   paste(rep(" ", maxstr - nchar(as.character(below))), 
                         below), below)
      text(xlimits[1], (xlimits[2] - 2), paste("\n", second.name, 
                                               ">", first.name, ": ", a1, "\n", second.name, 
                                               "=", first.name, ": ", o1, "\n", second.name, 
                                               "<", first.name, ": ", b1, sep = " "), pos = 4)
    }
    par(opar)
  }
  if (barplot == TRUE | chi == TRUE) {
    diffs = readings[, usecols[2]] - readings[, usecols[1]]
    maxdif = max(diffs)
    mindif = min(diffs)
    rangeofdifs = maxdif - mindif + 1
    compare = matrix(c(mindif:maxdif, rep(NA, rangeofdifs)), 
                     ncol = 2)
    for (i in 1:rangeofdifs) {
      compare[i, 2] = sum((readings[, usecols[2]] - readings[, 
                                                             usecols[1]]) == compare[i, 1])
    }
    if (barplot == TRUE) {
      opar = par(cex.lab = 1.3, cex.axis = 1.2)
      barplot(compare[, 2], names.arg = as.character(compare[, 
                                                             1]), xlab = paste(second.name, "-", first.name), 
              ylab = "Frequency", main = "Difference in readings", 
              las = 1)
      par(opar)
    }
    if (chi == TRUE) {
      if (correct == "Yates") {
        coeff = 0.5
      }
      else {
        if (correct == "Edwards") {
          coeff = 1
        }
        else return("continuity factor invalid (must by Yates or Edwards (in quotes))")
      }
      chisq = (above - below)^2/(above + below)
      chisq.cor = (abs(above - below) - coeff)^2/(above + below)
      signif = 1 - pchisq(chisq, 1)
      signif.cor = 1 - pchisq(chisq.cor, 1)
      needrows = 2 * max(abs(maxdif), abs(mindif)) + 1 ### this line has been changed
      if (maxdif == (-mindif)) 
        augmat = compare
      if (maxdif > abs(mindif)) {
        deficit = needrows - dim(compare)[1]
        deficit
        aug1 = (mindif - deficit):(mindif - 1)
        augmat = matrix(c(aug1, rep(0, deficit)), ncol = 2)
        augmat = rbind(augmat, compare)
      }
      if (maxdif < abs(mindif)) {
        deficit = abs(mindif) - maxdif
        aug1 = (maxdif + 1):abs(mindif)
        augmat = matrix(c(aug1, rep(0, deficit)), ncol = 2)
        augmat = rbind(compare, augmat)
      }
      numcomparisons = floor(dim(augmat)[1]/2)
      expected = NULL
      for (i in 1:numcomparisons) {
        expected[i] = (augmat[i, 2] + augmat[(needrows + 
                                                1 - i), 2])/2
      }
      cumexpected = cumsum(expected)
      cumexpected
      combine = which(cumexpected >= pool.criterion)[1]
      if (is.na(combine)) {
        message = paste("no cells left if you impose minimum expected cell count of", 
                        pool.criterion)
        return(message)
      }
      new_numcomparisons = numcomparisons - combine + 1
      pooled_augmat = augmat
      if (combine > 0) {
        a = c(NA, sum(augmat[1:combine, 2]))
        b = c(NA, sum(augmat[(dim(augmat)[1] - combine + 
                                1):dim(augmat)[1], 2]))
        pooled_augmat = rbind(a, augmat[(combine + 1):(dim(augmat)[1] - 
                                                         combine), ], b)
      }
      E.numcomparisons = floor(dim(pooled_augmat)[1]/2)
      E.nrows = dim(pooled_augmat)[1]
      E.chisq = NULL
      for (i in 1:E.numcomparisons) {
        E.expected = (pooled_augmat[i, 2] + pooled_augmat[(E.nrows + 
                                                             1 - i), 2])
        E.chisq[i] = (pooled_augmat[i, 2] - pooled_augmat[(E.nrows + 
                                                             1 - i), 2])^2/E.expected
      }
      E.chisquare = sum(E.chisq)
      E.df = length(E.chisq)
      E.pvalue = 1 - pchisq(E.chisquare, E.df)
    }
    sample_size = sum(compare[, 2])
    percentages = round(100 * compare[, 2]/sample_size, 1)
    compare = cbind(compare, percentages)
    colnames(compare) = c("difference", "frequency", "percentage")
    if (chi == TRUE) {
      if (cont.cor == TRUE) {
        MN = data.frame(Chisq = chisq, pvalue = signif)
        MNc = data.frame(Chisq = chisq.cor, pvalue = signif.cor, 
                         Correction = correct)
        EH = data.frame(Chisq = E.chisquare, df = round(E.df, 
                                                        0), pvalue = E.pvalue)
        return(list(thecall = thecall, McNemar = MN, 
                    McNemar_continuity_correction = MNc, Evans_Hoenig = EH, 
                    difference_frequency = compare, sample_size = sample_size))
      }
      else {
        MN = data.frame(Chisq = chisq, pvalue = signif)
        EH = data.frame(Chisq = E.chisquare, df = round(E.df, 
                                                        0), pvalue = E.pvalue)
        return(list(thecall = thecall, McNemar = MN, 
                    Evans_Hoenig = EH, difference_frequency = compare, 
                    sample_size = sample_size))
      }
    }
    else return(list(thecall = thecall, difference_frequency = compare, 
                     sample_size = sample_size))
  }
}