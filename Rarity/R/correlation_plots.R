######################################
######    Correlation plots    #######
######################################


corPlot <- function(df, method = "spearman", digits = 2, na.action = "keep", ties.method = "average",
                    title = "", xlab = "variable.name", ylab = "variable.name", ...)
{
  rankDf <- apply(df, 2, .rankings, ties.method = ties.method, na.action = na.action)
  
  if(length(xlab) > 1) if(length(xlab) != (ncol(df))) stop("Provide as many labels for x axes as there are variables, in the same order as the column order.")
  if(length(ylab) > 1) if(length(ylab) != (ncol(df))) stop("Provide as many labels for y axes as there are variables, in the same order as the column order.")
  
  crossCongR <- NULL
  crossCongP <- NULL
  for (i in seq_len(ncol(rankDf)))
  {
    r <- NULL
    p <- NULL
    for (j in 1:ncol(df))
    {
      r <- c(r, cor.test(df[, i], df[, j], method = method)$estimate)
      p <- c(p, cor.test(df[, i], df[, j], method = method)$p.value)
    }
    crossCongR <- rbind(crossCongR, r)
    crossCongP <- rbind(crossCongP, p)
  }
  colnames(crossCongP) <- colnames(rankDf)[1:ncol(rankDf)]
  rownames(crossCongP) <- colnames(rankDf)
  colnames(crossCongR) <- colnames(rankDf)[1:ncol(rankDf)]
  rownames(crossCongR) <- colnames(rankDf)
  
  crossCongR <- round(crossCongR, 2)
  
  
  par(mfrow = c(ncol(df), ncol(df)), oma = c(4.1, 4.1, 0.6, 0.6), mar = c(.5, .5, .5, .5))
  
  for (i in 1:ncol(df))
  {
    for (j in 1:ncol(df))
    {
      if (j == i)
      {
        plot(1, 1, pch ="", bty = "o", xaxt = "n", yaxt = "n")
        text(x = 1, y = 1, paste(ifelse(is.null(colnames(df)),
                                        paste("Var", i, sep = "."),
                                        colnames(df)[i]),
                                 "\nn = ", length(which(!is.na(df[, i])))), 
             font = 2)
      } else
        if (j > i)
        {
          plot(1, 1, pch ="", xaxt = "n", yaxt = "n", bty = "n")
          text(x = 1, y = 1, paste(formatC(crossCongR[i, j], format = "f", digits = digits), 
                                   ifelse(crossCongP[i, j] > 0.1, "",
                                          ifelse(crossCongP[i, j] > 0.05, "-",
                                                 ifelse(crossCongP[i, j] > 0.01, "*",
                                                        ifelse(crossCongP[i, j] > 0.001, "**", "***"))))), font = 2)
        } else
        {
          if (length(ylab) > 1)
          {
            if (j == 1) mtext(paste(ylab[i]), outer = T, side = 2, line = 2, at = (ncol(df)-i + .5)/ncol(df), cex = .8)
          } else if (ylab == "variable.name")
          {
            if (j == 1) mtext(paste(ifelse(is.null(colnames(rankDf)),
                                           ifelse(method == "pearson",
                                                  paste("Var. ", i, sep = ""),
                                                  paste("Var. ", i, " rank", sep = "")),
                                           ifelse(method == "pearson",
                                                  colnames(rankDf)[i],
                                                  paste(colnames(rankDf)[i], " rank", sep = "")))), outer = T, side = 2, line = 2, at = (ncol(df)-i + .5)/ncol(df), cex = .8)
          } else
          {
            if (j == 1) mtext(paste(ylab), outer = T, side = 2, line = 2, at = (ncol(df)-i + .5)/ncol(df), cex = .8)
          }
          if (length(xlab) > 1)
          {
            if (i == ncol(df)) mtext(paste(xlab[j]), outer = T, side = 1, line = 2, at = (j - .5)/ncol(df), cex = .8)
          } else if(xlab == "variable.name")
          {
            if (i == ncol(df)) mtext(paste(ifelse(is.null(colnames(rankDf)),
                                                  ifelse(method == "pearson",
                                                         paste("Var. ", j, sep = ""),
                                                         paste("Var. ", j, " rank", sep = "")),
                                                  ifelse(method == "pearson",
                                                         colnames(rankDf)[j],
                                                         paste(colnames(rankDf)[j], " rank", sep = "")))), outer = T, side = 1, line = 2, at = (j - .5)/ncol(df), cex = .8)
          } else
          {
            if (i == ncol(df)) mtext(paste(xlab), outer = T, side = 1, line = 2, at = (j - .5)/ncol(df), cex = .8)
          }
          if(method == "spearman" | method == "kendall")
          {
            plot(rankDf[, i] ~ rankDf[, j], xaxt = "n", yaxt = "n", xlim = c(0, nrow(df)), ylim = c(0, nrow(df)), xlab = "", ylab = "", ...)
            if (i == ncol(df)) axis(1, at = c(0, round(nrow(df)/2, 0), nrow(df)), las = 1)
            if (j == 1) axis(2, at = c(0, round(nrow(df)/2, 0), nrow(df)), las = 1)
          } else if (method == "pearson")
          {
            dx <- diff(c(min(df[, j], na.rm = T), max(df[, j], na.rm = T)))
            dy <- diff(c(min(df[, i], na.rm = T), max(df[, i], na.rm = T)))
            plot(df[, i] ~ df[, j], xaxt = "n", yaxt = "n", xlim = c(min(df[, j], na.rm = T) - 0.05 * dx, 
                                                                     max(df[, j], na.rm = T) + 0.05 * dx), 
                 ylim = c(min(df[, i], na.rm = T) - 0.05 * dy, 
                          max(df[, i], na.rm = T) + 0.05 * dy), 
                 xlab = "", ylab = "", ...)
            if (i == ncol(df)) axis(1, las = 1)
            if (j == 1) axis(2, las = 1)
          }
        }
    }
  }
  mtext(paste(title), outer = T, side = 3, line = -1, at = .55, cex = 1.2)
}

.rankings <- function(x, ties.method = "average", na.action = "keep") length(x) - length(which(is.na(x))) + 1 - rank(x, ties.method = ties.method, na.last = na.action)
