#' Generate IndX and mR charts.
#'
#' @references Engineering Statistics Handbook 6.3.2, NIST/SEMATECH e-Handbook of Statistical Methods National Institute of Standards and Technology, Dec 2006
#' @references https://en.wikipedia.org/wiki/Nelson_rules
#' @references Lloyd S. Nelson, "Technical Aids," Journal of Quality Technology 16, no. 4 (October 1984), 238-239.
#' @references The 8 rules are:
#' 1 One point is more than 3 standard deviations from the mean.
#' 2 Nine (or more) points in a row are on the same side of the mean.
#' 3 Six (or more) points in a row are continually increasing (or decreasing).
#' 4 Fourteen (or more) points in a row alternate in direction, increasing then decreasing.
#' 5 Two (or three) out of three points in a row are more than 2 standard deviations from the mean in the same direction.
#' 6 Four (or five) out of five points in a row are more than 1 standard deviation from the mean in the same direction.
#' 7 Fifteen points in a row are all within 1 standard deviation of the mean on either side of the mean.
#' 8 Eight points in a row exist with none within 1 standard deviation of the mean and the points are in both directions from the mean.
#' @param x (mandatory) A data frame with the individual values in the first column and the time in the second column. It can be either a factor, a date or a string and it will be ordered automatically. \cr See \code{?spcTimeSeries}
#' @param linesColors (optional) A vector with 7 colors in order from the average + 3 standard deviations to the average - 3 standard deviations, including the average itself in the center. \cr Default value is \code{c("gray50", "gray65", "gray85", "black", "gray85", "gray65", "gray50")}
#' @param applyRules (optional) A vector with 8 boolean values indicating which rules must be applied. \cr Default value is \code{c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)}
#' @param rulesColors (optional) A vector with colors, one for each rule. The last point of each violating run will be colored indicating the violation and the corresponding rule. \cr Default value is \code{c("red", "yellow2", "green", "magenta", "blue", "orange", "brown", "cyan")}
#' @param seg (optional) A vector with the positions of the points where there should be breaks and another pair of charts should be plotted. It may be used for better visualization when the series is too long. \cr Default value is \code{c()}
#' @param keepStats (optional) A boolean indicating if each segment's plot should be considered as part of the same series or independelty, as different series. If TRUE, it will be considered as part of the same series. If FALSE, each plot will have it's limits calculated independently, as well as the application of the rules. It's useful to compare different scenarios. \cr Default value is \code{TRUE}
#' @param verbose (optional) A boolean indicating if mean, standard deviation/UCL and number of violations should be printed. \cr Default value is \code{FALSE}
#' @examples
#' data("spcTimeSeries")
#' six_sigma_ctrl_chart(spcTimeSeries)
#' six_sigma_ctrl_chart(spcTimeSeries, verbose=TRUE)
#' six_sigma_ctrl_chart(spcTimeSeries, seg=c(25, 50, 75))
#' six_sigma_ctrl_chart(spcTimeSeries, seg=c(25, 50, 75), keepStats=FALSE, verbose=TRUE)
#' @return None
#' @export
six_sigma_ctrl_chart <- function(x,
                                 linesColors = c("gray50", "gray65", "gray85", "black", "gray85", "gray65", "gray50"),
                                 applyRules = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                                 rulesColors = c("red", "yellow2", "green", "magenta", "blue", "orange", "brown", "cyan"),
                                 seg = c(), keepStats=TRUE, verbose=FALSE) {
  x = x[order(x[2]),]
  x[,2] <- as.factor(x[,2])

  if (length(seg) == 0) {
    mR = abs(x[,1][1:(length(x[,1])-1)] - x[,1][2:length(x[,1])])
    meanmR = mean(mR)
    mean = mean(x[,1])
    sd = sd(x[,1])
    .plotSixSigma(x, mR, linesColors, rulesColors, applyRules, meanmR, mean, sd, verbose)
  } else {
    if (seg[1] != 1) {
      seg = append(seg, 1, after=0)
    }
    if (seg[length(seg)] != length(x)) {
      seg = append(seg, length(x[,1]))
    }
    seg[1] <- -1
    for (i in 1:(length(seg)-1)) {
      grDevices::dev.new()
      xseg = x[c((seg[i]+1):seg[i+1]),]
      if (keepStats) {
        mR = abs(x[,1][1:(length(x[,1])-1)] - x[,1][2:length(x[,1])])
        meanmR = mean(mR)
        mean = mean(x[,1])
        sd = sd(x[,1])
        mRseg = abs(xseg[,1][1:(length(xseg[,1])-1)] - xseg[,1][2:length(xseg[,1])])
        .plotSixSigma(xseg, mRseg, linesColors, rulesColors, applyRules, meanmR, mean, sd, verbose)
      } else {
        mR = abs(xseg[,1][1:(length(xseg[,1])-1)] - xseg[,1][2:length(xseg[,1])])
        meanmR = mean(mR)
        mean = mean(xseg[,1])
        sd = sd(xseg[,1])
        .plotSixSigma(xseg, mR, linesColors, rulesColors, applyRules, meanmR, mean, sd, verbose)
      }
    }
  }
}

.plotSixSigma <- function (x, mR, linesColors, rulesColors, applyRules, meanmR, mean, sd, verbose) {
  zonesBoundaries = .findZones(x[,1], mean, sd)

  graphics::par(mfrow=c(2,1))

  graphics::plot(x[,1], xlim=c(x[,2][1],x[,2][length(x[,2])]), ylab="IndX", xlab="", cex.lab=0.8, yaxt="n", xaxt="n", ylim=c(min(zonesBoundaries), max(zonesBoundaries)), type="n", pch=16, cex=0.7)
  graphics::axis(x[,2],side=1, las=2, at=x[,2], cex.axis=0.8)
  graphics::axis(side=2, cex.axis=0.8)

  for (i in 1:7) {
    graphics::abline(zonesBoundaries[,i][1] + sd, b=0, col=linesColors[i])
    graphics::mtext(at=zonesBoundaries[,8-i][1] + sd,text=paste(" ",round(zonesBoundaries[,8-i][1] + sd,digits=2)), side=4, cex = 0.8, las=1)
  }

  colors <- .paintViolators(x[,1], rulesColors, applyRules, zonesBoundaries)

  graphics::points(y=x[,1], x=x[,2], col=colors$color, pch=16, cex=0.7)

  for (i in 1:(length(colors)-1)) {
    for (j in 1:length(colors$color)) {
      if (colors[j,i]) {
        graphics::text(x=x[,2][j], y=x[,1][j], col = rulesColors[i], labels = round(x[,1][j],2), pos=1, cex=0.8)
      }
    }
  }

  graphics::mtext(text=paste("Mean = ", round(mean,2)), side=1, padj=7, adj=0, cex = 0.8)
  graphics::mtext(text=paste("Std. dev.  = ", round(sd,2)), side=1, padj=8, adj=0, cex = 0.8)
  graphics::mtext(text=paste("Violations = ", sum(colors == TRUE)), side=1, padj=9, adj=0, cex = 0.8)

  if (verbose) {
    print(paste("IndX Mean = ", mean))
    print(paste("IndX Std. dev.  = ", sd))
    print(paste("IndX Violations = ", sum(colors == TRUE)))
  }

  ucl = 3.267*meanmR
  violations <- ifelse(mR > ucl, rulesColors[1], "black")

  graphics::plot(ylab="mR", xlab="", mR, yaxt="n", ylim=c(-1, max(ucl, max(mR)+1)), xaxt="n", type="n", cex.lab = 0.8)
  graphics::axis(side=2, cex.axis=0.8)
  graphics::axis(side=1, las=2, at=c(1:length(mR)), cex.axis=0.8)

  graphics::abline(a=ucl,b=0, col=linesColors[1])
  graphics::mtext(at=ucl,text=paste(" ",round(ucl,digits=2)), side=4, cex = 0.8, las=1)

  graphics::points(mR, col=violations, pch=16, cex=0.7)

  for (i in 1:length(violations)) {
    if (violations[i] != "black") {
      graphics::text(x=i, y=mR[i], col = rulesColors[1], labels = round(mR[i],2), pos=1, cex=0.8)
    }
  }

  graphics::mtext(text=paste("Mean = ", round(mean(mR),2)), side=1, padj=3, adj=0, cex = 0.8)
  graphics::mtext(text=paste("UCL  = ", round(ucl,2)), side=1, padj=4, adj=0, cex = 0.8)
  graphics::mtext(text=paste("Violations = ", sum(violations == rulesColors[1])), side=1, padj=5, adj=0, cex = 0.8)

  if (verbose) {
    print(paste("mR Mean = ", mean(mR)))
    print(paste("mR UCL  = ", ucl))
    print(paste("mR Violations = ", sum(violations == rulesColors[1])))
  }
}

.paintViolators <- function(x, rulesColors, applyRules, zones) {
	zones <- rowSums(x >= zones) - 4
	results <- plyr::ldply(1:length(x), function(i) {
		.runTests(x, zones, i, applyRules)
	})

	results$color <- ifelse(results$rule1!=0, rulesColors[1],
					ifelse(results$rule2!=0, rulesColors[2],
						ifelse(results$rule3!=0, rulesColors[3],
							ifelse(results$rule4!=0, rulesColors[4],
								ifelse(results$rule5!=0, rulesColors[5],
									ifelse(results$rule6!=0, rulesColors[6],
									  ifelse(results$rule7!=0, rulesColors[7],
									    ifelse(results$rule8!=0, rulesColors[8],
										    "black"))))))))
	results
}

.findZones <- function(x, mean, sd) {
  boundaries <- seq(-4, 4)
  zones <- sapply(boundaries, function(i) {
    i * rep(sd, length(x)) + mean
  })
  zones
}

.runTests <- function(x, zones, i, applyRules) {
  if (applyRules[1]) {
    values <- zones[i]
    rule1 <- any(values > 3) || any(values < -2)
  } else {
    rule1 = FALSE
  }
  if (applyRules[2]) {
    values <- zones[max(i-8, 1):i] # Rule 2
    rule2 <- length(values) == 9 && (all(values > 0) || all(values < 1))
  } else {
    rule2 = FALSE
  }
  if (applyRules[3]) {
    values <- x[max((i-5),1):i] # Rule 3
    rule3 <- length(values) == 6 && (all(values == cummax(values)) || all(values == cummin(values)))
  } else {
    rule3 = FALSE
  }
  if (applyRules[4]) {
    values <- x[max(i-13, 1):i] # Rule 4
    rule4 <- length(values) == 14 && (all((values[1:(length(values)-1)] - values[2:length(values)] < 0) == rep(c(TRUE, FALSE), length.out=13)) || all((values[1:(length(values)-1)] - values[2:length(values)] < 0) == rep(c(FALSE, TRUE), length.out=13)))
  } else {
    rule4 = FALSE
  }
  if (applyRules[5]) {
    values <- zones[max(i-2, 1):i] # Rule 5
    rule5 <- length(values) == 3 && (sum(values >= 3) >= 2 || sum(values <= -2) >= 2)
  } else {
    rule5 = FALSE
  }
  if (applyRules[6]) {
    values <- zones[max(i-4, 1):i] # Rule 6
    rule6 <- length(values) == 5 && (sum(values >= 2) >= 4 || sum(values <= -1) >= 4)
  } else {
    rule6 = FALSE
  }
  if (applyRules[7]) {
    values <- zones[max(i-14, 1):i] # Rule 7
    rule7 <- length(values) == 15 && sum(values == 0) + sum( values == 1) == 15
  } else {
    rule7 = FALSE
  }
  if (applyRules[8]) {
    values <- zones[max(i-7, 1):i] # Rule 8
    rule8 <- sum(values < 0) + sum(values > 1) == 8
  } else {
    rule8 = FALSE
  }

  c("rule1"=rule1, "rule2"=rule2, "rule3"=rule3, "rule4"=rule4, "rule5"=rule5, "rule6"=rule6, "rule7"=rule7, "rule8"=rule8)
}
