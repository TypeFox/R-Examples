dandelion <- function (fact_load, bound = 0.5, mcex = c(1, 1), palet) 
{
  if (class(fact_load) != "loadings") {
    cat(" Example : dandelion(loading object,bound=0)\n")
    stop("please use a loadings object")
  }
  if ((bound < 0) || (bound > 1)) 
    stop("bound must be between 0 and 1")
  load_grid <- function(fact_load1, fact_load2, x3, y3, x2, 
                        y2, tempx, tempy, col1) {
    coln <- length(col1)
    col2 <- rev(col1)
    x4 <- seq(tempx, x3, by = (x3 - tempx) * (1/coln))
    y4 <- seq(tempy, y3, by = (y3 - tempy) * (1/coln))
    if (fact_load1 > 0 && fact_load2 < 0) {
      for (k in 2:coln) polygon(c(x4[k], x2, x4[k - 1]), 
                                c(y4[k], y2, y4[k - 1]), col = col2[k], border = col2[k])
    }
    else if (fact_load1 < 0 && fact_load2 > 0) {
      for (k in 2:coln) polygon(c(x4[k], x2, x4[k - 1]), 
                                c(y4[k], y2, y4[k - 1]), col = col1[k], border = col1[k])
    }
    else if (fact_load1 < 0 && fact_load2 < 0) 
      polygon(c(x3, x2, tempx), c(y3, y2, tempy), col = col2[1], 
              border = col2[1])
    else if (fact_load1 > 0 && fact_load2 > 0) 
      polygon(c(x3, x2, tempx), c(y3, y2, tempy), col = col1[1], 
              border = col1[1])
  }
  old_par <- par(no.readonly = TRUE)
  factor <- ncol(fact_load)
  commun_fact <- apply(fact_load, 1, function(x) sum(x^2))
  lambda <- apply(fact_load, 2, function(x) sum(x^2))
  lambda2 <- sort(lambda, decreasing = TRUE)
  count <- NULL
  for (i in lambda2) count <- c(count, which(i == lambda))
  count <- unique(count)
  fact_load <- fact_load[, count]
  lambda <- lambda2
  unique_fact <- 1 - commun_fact
  aci <- 360 * (lambda/nrow(fact_load))
  if ((max(aci)/2) > (360 - sum(aci))) 
    aci <- (aci/360) * (360 - (max(aci)/2))
  degreef <- (270 + c(0, cumsum(aci)))%%360
  if (aci[1] > 180) 
    aci[1] <- 360 - aci[1]
  limitt <- sqrt(2 * (1 - cos(aci * pi/180)))/2
  limit <- limitt * 0.9
  maxcex <- (lambda/max(lambda)) * mcex[1]
  count <- NULL
  count2 <- NULL
  count3 <- NULL
  for (i in 1:factor) {
    for (j in 1:nrow(fact_load)) {
      if (abs(fact_load[j, i]) == max(abs(fact_load[j, 
                                                    ]))) {
        count <- c(count, j)
        if (length(count2) == (i - 1)) 
          count2 <- c(count2, j)
      }
    }
  }
  fact_load <- fact_load[count, ]
  commun_fact <- commun_fact[count]
  unique_fact <- 1-commun_fact
  for (i in 1:factor) {
    temp <- which(count == count2[i])
    if (length(temp) != 0) 
      count3 <- c(count3, which(count == count2[i]))
  }
  xlimit <- cos(degreef * pi/180)
  temp <- c(limitt, 0)
  limord <- which(xlimit < 0)
  temp[limord] <- -1 * temp[limord]
  xlimit <- cos(degreef * pi/180) + temp
  ylimit <- sin(-degreef * pi/180)
  temp <- c(limitt, 0)
  limord <- which(ylimit < 0)
  temp[limord] <- -1 * temp[limord]
  ylimit <- sin(-degreef * pi/180) + temp
  xmax <- max(xlimit)
  xmin <- min(xlimit)
  ymax <- max(ylimit)
  ymin <- min(ylimit)
  if (max(limitt) > xmax) 
    xmax <- max(limitt)
  if (max(limitt) > abs(xmin)) 
    xmin <- -1 * max(limitt)
  if (max(limitt) > ymax) 
    ymax <- max(limitt)
  if (max(limitt) > abs(ymin)) 
    ymin <- -1 * max(limitt)
  par(mar = c(0, 0, 0, 0))
  layout(matrix(c(1, 1, 1, 1, 2, 3, 4, 5), nrow = 4), widths = c(7, 
                                                                 3), heights = c(0.5, 1, 1, 1))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  degreev <- (270 + 0:(nrow(fact_load) - 1) * (360/nrow(fact_load)))%%360
  x2 = cos(degreef * pi/180)
  y2 = sin(-degreef * pi/180)
  srt2 <- 360 - degreev
  bloc <- which((degreev - 270)%%360 > 180)
  srt2[bloc] <- srt2[bloc] - 180
  for (i in 1:factor) {
    lines(c(0, x2[i]), c(0, y2[i]), type = "l", col = "black")
    datanew <- abs(fact_load[, i])
    bloc <- which(datanew <= bound)
    if (length(bloc) > 0) {
      datanew[bloc] <- rep(0, length(datanew[bloc]))
      datanew[-bloc] <- (datanew[-bloc] - bound)/(1 - 
                                                    bound)
    }
    else datanew <- (datanew - bound)/(1 - bound)
    xm = x2[i] + limit[i] * cos(degreev * pi/180)
    ym = y2[i] + limit[i] * sin(-degreev * pi/180)
    x3 = x2[i] + (limit[i] * datanew) * (cos(degreev * pi/180))
    y3 = y2[i] + (limit[i] * datanew) * (sin(-degreev * 
                                               pi/180))
    if (fact_load[1, i] > 0) 
      col2 = palet[1]
    if (fact_load[1, i] < 0) 
      col2 = palet[length(palet)]
    polygon(c(x3[1], x2[i]), c(y3[1], y2[i]), col = col2, 
            border = col2)
    lines(c(x3[1], xm[1]), c(y3[1], ym[1]), type = "l", 
          col = "grey")
    for (j in 2:(nrow(fact_load) + 1)) {
      jnew <- ((j - 1)%%(nrow(fact_load))) + 1
      load_grid(fact_load[jnew, i], fact_load[j - 1, i], 
                x3[jnew], y3[jnew], x2[i], y2[i], x3[j - 1], 
                y3[j - 1], palet)
      lines(c(x3[jnew], xm[jnew]), c(y3[jnew], ym[jnew]), 
            type = "l", col = "grey")
    }
    if ((length(count3) + 1) > i) {
      if (length(count3) == i) 
        text_space <- count3[i]:nrow(fact_load)
      else text_space <- count3[i]:(count3[i + 1] - 1)
      for (k in text_space) {
        x4 = x2[i] + limitt[i] * cos(degreev[k] * pi/180)
        y4 = y2[i] + limitt[i] * sin(-degreev[k] * pi/180)
        text(x4, y4, paste(abbreviate(rownames(fact_load)[k])), 
             cex = maxcex[i], srt = srt2[k])
      }
    }
  }
  x2 = cos(degreef[factor + 1] * pi/180)
  y2 = sin(-degreef[factor + 1] * pi/180)
  lines(c(0, x2), c(0, y2), type = "l", col = "black", lty = 2)
  plot(1:10, type = "n", xlim = c(-1.5, 1.5), ylim = c(-1.5, 
                                                       1.5), axes = FALSE)
  legend(0, 0, c("pos. load.", "neg. load."), 
         col = c(palet[1], palet[length(palet)]), 
         text.col = "black", bg = "white", bty="n",
         xjust = 0.5, yjust = 0.5, pch = c(15, 15), cex = 1.5)
  par(mar = c(0, 0, 0.9, 0))
  plot(1:10, type = "n", xlim = c(-1.5, 1.5), ylim = c(-1.5, 
                                                       1.5), axes = FALSE)
  title(main = "uniquenesses", cex.main = mcex[2])
  x3 = cos(degreev * pi/180)
  y3 = sin(-degreev * pi/180)
  x4 = 1.25 * cos(degreev * pi/180)
  y4 = 1.25 * sin(-degreev * pi/180)
  x5 = unique_fact * cos(degreev * pi/180)
  y5 = unique_fact * sin(-degreev * pi/180)
  for (i in 1:nrow(fact_load)) {
    lines(c(0, x3[i]), c(0, y3[i]), type = "l", col = "grey")
    text(x4[i], y4[i], paste(abbreviate(rownames(fact_load)[i])), 
         cex = mcex[2], srt = srt2[i])
  }
  polygon(x5, y5, col = palet[1], border = palet[1])
  plot(1:10, type = "n", xlim = c(-1.5, 1.5), ylim = c(-1.5, 
                                                       1.5), axes = FALSE)
  title(main = "communalities", cex.main = mcex[2])
  x5 = commun_fact * cos(degreev * pi/180)
  y5 = commun_fact * sin(-degreev * pi/180)
  for (i in 1:nrow(fact_load)) {
    lines(c(0, x3[i]), c(0, y3[i]), type = "l", col = "grey")
    text(x4[i], y4[i], paste(abbreviate(rownames(fact_load)[i])), 
         cex = mcex[2], srt = srt2[i])
  }
  polygon(x5, y5, col = palet[1], border = palet[1])
  par(mar = c(0, 0, 0, 0))
  lam <- cbind(round(lambda, digits = 2), round(cumsum(lambda/nrow(fact_load)), 
                                                digits = 2))
  lam <- round(lambda, digits = 2)
  r.lam <- round(cumsum(lambda/nrow(fact_load)), digits = 2)
  names(r.lam) <- paste("F", 1:factor, sep = ".")
  par(mar=c(2, 2, 2, 2))
  bar.g <- barplot(r.lam,mgp=c(3,0.5,0),ylim=c(0,1),main="Cum. Ratio")
  text(bar.g,r.lam,lam,pos=3, offset=.2,font=2)
  par(old_par)
  invisible()
}
