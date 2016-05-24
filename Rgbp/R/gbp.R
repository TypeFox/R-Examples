######
gbp <- function(y, se.or.n, covariates, mean.PriorDist, model, intercept, confidence.lvl, n.AR, n.AR.factor,
                trial.scale, save.result, normal.CI, t, u) UseMethod("gbp")
######
gbp <- function(y, se.or.n, covariates, mean.PriorDist, model = "gaussian", 
                        intercept = TRUE, confidence.lvl = 0.95, 
                        n.AR = 0, n.AR.factor = 4, trial.scale = NA, save.result = TRUE, 
                        normal.CI = FALSE, t = 0, u = 1) {

  ##input checks
  if(model == "poisson" & missing(mean.PriorDist)) {
    warning("Model is Poisson and the prior mean is unknown. This program can not yet give reliable results for this data. If the Poisson is being used as an approximation to the Binomial and the exposures are known then please assume a Binomial model")
  }  
  
  if(intercept == "FALSE" & missing(covariates)) {
    warning("When there are no covariates and no intercept term, model cannot fit the model. Please set the argument mean.PriorDist or use the intercept term.")
    stop()
  }  
  
  ######
  res <- switch(model, 
       gaussian = gr(y, se.or.n, X = covariates, mu = mean.PriorDist, confidence.lvl = confidence.lvl, intercept = intercept, 
                     normal.CI = normal.CI), 
       binomial = br(y, se.or.n, X = covariates, prior.mean = mean.PriorDist, intercept = intercept, confidence.lvl = confidence.lvl,
                     n.AR = n.AR, n.AR.factor = n.AR.factor, trial.scale = trial.scale, save.result = TRUE,
                     t = t, u = u), 
       poisson = pr(y, se.or.n, X = covariates, prior.mean = mean.PriorDist, intercept = intercept, confidence.lvl = confidence.lvl))
  
  class(res) <- "gbp"	
  res
}

print.gbp <- function(x, sort = TRUE, ...) {
  
  if (any(is.na(x$prior.mean)) & !identical(x$X, NA)) {

	cova <- as.matrix(x$X)
	colnames(cova) <- paste(1 : dim(cova)[2])
    
    if (x$model == "gr") {
      temp <- data.frame(obs.mean = x$sample.mean, se = x$se, x = cova, 
                         prior.mean = x$prior.mean.hat, shrinkage = x$shrinkage, 
                         low.intv = x$post.intv.low, post.mean = x$post.mean, 
                         upp.intv = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(obs.mean = x$sample.mean, n = x$se, x = cova, 
                         prior.mean = x$prior.mean.hat, shrinkage = x$shrinkage, 
                         low.intv = x$post.intv.low, post.mean = x$post.mean, 
                         upp.intv = x$post.intv.upp, post.sd = x$post.sd)
    }

  } else if (any(is.na(x$prior.mean)) & identical(x$X, NA)) { 
  # if there are neither prior.mean and X assigned

    if (x$model == "gr") {
      temp <- data.frame(obs.mean = x$sample.mean, se = x$se, 
                         prior.mean = x$prior.mean.hat, shrinkage = x$shrinkage, 
                         low.intv = x$post.intv.low, post.mean = x$post.mean, 
                         upp.intv = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(obs.mean = x$sample.mean, n = x$se, 
                         prior.mean = x$prior.mean.hat, shrinkage = x$shrinkage, 
                         low.intv = x$post.intv.low, post.mean = x$post.mean, 
                         upp.intv = x$post.intv.upp, post.sd = x$post.sd)
    }


  } else if (!any(is.na(x$prior.mean))) {
    if (x$model == "gr") {
      temp <- data.frame(obs.mean = x$sample.mean, se = x$se, 
                         prior.mean = x$prior.mean, shrinkage = x$shrinkage, 
                         low.intv = x$post.intv.low, post.mean = x$post.mean, 
                         upp.intv = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(obs.mean = x$sample.mean, n = x$se, 
                         prior.mean = x$prior.mean, shrinkage = x$shrinkage, 
                         low.intv = x$post.intv.low, post.mean = x$post.mean, 
                         upp.intv = x$post.intv.upp, post.sd = x$post.sd)
    }
  }
  
  if (sort == TRUE) {
    if (x$model == "gr") {
      temp <- temp[order(temp[, 2], decreasing = TRUE), ]
    } else {
      temp <- temp[order(temp[, 2]), ]
    }
  }

  temp.mean <- colMeans(temp)
  temp <- data.frame(rbind(temp, temp.mean), row.names = c(rownames(temp), "Mean"))
  temp[, 1] <- format.default(temp[, 1], digits = 3)
  temp[dim(temp)[1], 1] <- ""

  if (x$model == "gr") {
    temp[, 2] <- round(temp[, 2], 1)
  } else {
    temp[, 2] <- round(temp[, 2])
  }

  if (any(is.na(x$prior.mean)) & !identical(x$X, NA)) {
    temp[, 3] <- round(temp[, 3], 2)
  }

  if (sort == TRUE) {
    if (x$model == "gr") {
      cat("Summary for each group (sorted by the descending order of se): \n")
    } else {
      cat("Summary for each group (sorted by  the ascending order of n): \n")
    }
  } else {
    cat("Summary for each group: \n")
  }

  cat("\n")
  options(digits = 3)
  print(temp)
  options(digits = 7)
}

summary.gbp <- function(object, ...) {

  if (any(is.na(object$prior.mean)) & !identical(object$X, NA)) {

	cova <- as.matrix(object$X)
	colnames(cova) <- paste(1 : dim(cova)[2])
    
    if (object$model == "gr") {
      temp <- data.frame(obs.mean = object$sample.mean, se = object$se, x = cova, 
                         prior.mean = object$prior.mean.hat, shrinkage = object$shrinkage, 
                         low.intv = object$post.intv.low, post.mean = object$post.mean, 
                         upp.intv = object$post.intv.upp, post.sd = object$post.sd)
    } else {
      temp <- data.frame(obs.mean = object$sample.mean, n = object$se, x = cova, 
                         prior.mean = object$prior.mean.hat, shrinkage = object$shrinkage, 
                         low.intv = object$post.intv.low, post.mean = object$post.mean, 
                         upp.intv = object$post.intv.upp, post.sd = object$post.sd)
    }

  } else if (any(is.na(object$prior.mean)) & identical(object$X, NA)) { 
  # if there are neither prior.mean and X assigned

    if (object$model == "gr") {
      temp <- data.frame(obs.mean = object$sample.mean, se = object$se, 
                         prior.mean = object$prior.mean.hat, shrinkage = object$shrinkage, 
                         low.intv = object$post.intv.low, post.mean = object$post.mean, 
                         upp.intv = object$post.intv.upp, post.sd = object$post.sd)
    } else {
      temp <- data.frame(obs.mean = object$sample.mean, n = object$se, 
                         prior.mean = object$prior.mean.hat, shrinkage = object$shrinkage, 
                         low.intv = object$post.intv.low, post.mean = object$post.mean, 
                         upp.intv = object$post.intv.upp, post.sd = object$post.sd)
    }


  } else if (!any(is.na(object$prior.mean))) {
    if (object$model == "gr") {
      temp <- data.frame(obs.mean = object$sample.mean, se = object$se, 
                         prior.mean = object$prior.mean, shrinkage = object$shrinkage, 
                         low.intv = object$post.intv.low, post.mean = object$post.mean, 
                         upp.intv = object$post.intv.upp, post.sd = object$post.sd)
    } else {
      temp <- data.frame(obs.mean = object$sample.mean, n = object$se, 
                         prior.mean = object$prior.mean, shrinkage = object$shrinkage,
                         low.intv = object$post.intv.low, post.mean = object$post.mean, 
                         upp.intv = object$post.intv.upp, post.sd = object$post.sd)
    }
  }



  if (min(object$se) == max(object$se)) {  # if se or n are all the same

    temp2 <- temp[order(temp$obs.mean), ]

    if (length(object$se) %% 2) {  # if number of groups is odd
      summary.table <- temp2[c(1, (dim(temp2)[1] + 1) / 2, dim(temp2)[1]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      if (number.of.medians == 1) {
        row.names(summary.table) <- c("Group with min(obs.mean)", 
                                      "Group with median(obs.mean)", 
                                      "Group with max(obs.mean)")
      } else {  # if there are more than one median
        row.names(summary.table) <- c("Group with min(obs.mean)", 
                                     paste("Group with median(obs.mean)", 1 : number.of.medians, sep = ""),
                                     "Group with max(obs.mean)")
      }
    } else {  # if number of groups is even
      summary.table <- temp2[c(1, dim(temp2)[1] / 2, dim(temp2)[1] / 2 + 1, dim(temp2)[1]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      if (number.of.medians == 1) {
        row.names(summary.table) <- c("Group with min(obs.mean)", 
                                      "Group with median(obs.mean)", 
                                      "Group with max(obs.mean)")
      } else {  # if there are more than one median
        row.names(summary.table) <- c("Group with min(obs.mean)", 
                                     paste("Group with median(obs.mean)", 1 : number.of.medians, sep = ""),
                                     "Group with max(obs.mean)")
      }
    }
  } else { # if n or se are different from each group

    temp2 <- temp[order(temp[, 2], temp$obs.mean), ]

    if (length(object$se) %% 2) {  # if number of groups is odd
      summary.table <- temp2[c(1, (dim(temp2)[1] + 1) / 2, dim(temp2)[1]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      if (object$model == "gr") {
        if (number.of.medians == 1) {
          row.names(summary.table) <- c("Group with min(se)", "Group with median(se)", "Group with max(se)")
        } else {
          row.names(summary.table) <- c("Group with min(se)", 
                                       paste("Group with median(se)", 1 : number.of.medians, sep = ""),
                                       "Group with max(se)")
        }
      } else {  # if model is not "gr"
        if (number.of.medians == 1) {
          row.names(summary.table) <- c("Group with min(n)", "Group with median(n)", "Group with max(n)")
        } else {
          row.names(summary.table) <- c("Group with min(n)", 
                                       paste("Group with median(n)", 1 : number.of.medians, sep = ""),
                                       "Group with max(n)")
        }
      }

    } else {   # if number of groups is even
      summary.table <- temp2[c(1, dim(temp2)[1] / 2, dim(temp2)[1] / 2 + 1, dim(temp2)[1]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      if (object$model == "gr") {
        if (number.of.medians == 1) {
          row.names(summary.table) <- c("Group with min(se)", "Group with median(se)", "Group with max(se)")
        } else {
          row.names(summary.table) <- c("Group with min(se)", 
                                       paste("Group with median(se)", 1 : number.of.medians, sep = ""),
                                       "Group with max(se)")
        }
      } else {  # if model is not "gr"
        if (number.of.medians == 1) {
          row.names(summary.table) <- c("Group with min(n)", "Group with median(n)", "Group with max(n)")
        } else {
          row.names(summary.table) <- c("Group with min(n)", 
                                       paste("Group with median(n)", 1 : number.of.medians, sep = ""),
                                       "Group with max(n)")
        }
      }
    }
  }

  temp.mean <- colMeans(temp)
  
  if (object$model == "gr") {
    temp.mean[2] <- round(temp.mean[2], 1)
    summary.table[, 2] <- round(summary.table[, 2], 1)
  } else {
    temp.mean[2] <- round(temp.mean[2])
    summary.table[, 2] <- round(summary.table[, 2])
  }

  summary.table <- data.frame(rbind(summary.table, temp.mean), 
                              row.names = c(rownames(summary.table), "Overall Mean"))
  summary.table[, 1] <- format.default(summary.table[, 1], digits = 3)
  summary.table[dim(summary.table)[1], 1] <- ""

  post.mode.alpha <- object$a.new
  post.sd.alpha <- sqrt(object$a.var)
  if (object$model == "gr") {
    post.mode.A <- exp(object$a.new)
    result2 <- data.frame(post.mode.alpha, post.sd.alpha, post.mode.A, row.names = "")
  } else {
    post.mode.r <- exp(-object$a.new)
    result2 <- data.frame(post.mode.alpha, post.sd.alpha, post.mode.r, row.names = "")
  }

  if (object$model == "br" & length(object$weight) != 1) {
    post.median.r <- median(exp(-object$alpha))
    post.median.alpha <- median(object$alpha)
    post.sd.alpha <- sd(object$alpha)
    result2 <- data.frame(post.median.alpha, post.sd.alpha, post.median.r, row.names = "")
  }
  
  if (any(is.na(object$prior.mean))) {
    estimate <- as.vector(object$beta.new)
    names(estimate) <- paste("beta", 1 : length(estimate), sep = "")

    if ((sum(is.na(object$weight)) == 1)) {
      se <- as.vector(sqrt(diag(object$beta.var)))
    } else {
      se <- as.vector(sqrt(object$beta.var))
    }
    z.val <- estimate / se
    p.val <- 2 * pnorm(-abs(z.val))
    beta.result <- data.frame(estimate, se, z.val, p.val)
    res <- list(main = summary.table, sec.var = result2, reg = beta.result)
  } else {
    res <- list(main = summary.table, sec.var = result2, reg = NA)
  }
  
  class(res) <- "summary.gbp"
  res
}


print.summary.gbp <- function(x, ...) {

  if (identical(x$reg, NA)) {
    options(digits = 3)
    cat("Main summary:\n")
    cat("\n")
    print(x$main)
    cat("\n")
    cat("\n")
    cat("Estimation summary for the second-level variance component:\n")
    cat("alpha = log(A) for Gaussian or alpha = log(1/r) for Binomial and Poisson data:\n")
    cat("\n")
    print(x$sec.var)
    options(digits = 7)
  } else {
    cat("Main summary:\n")
    cat("\n")
    options(digits = 3)
    print(x$main)
    cat("\n")
    cat("\n")
    cat("Estimation summary for the second-level variance component:\n")
    cat("alpha = log(A) for Gaussian or alpha =  log(1/r) for Binomial and Poisson data:\n")
    cat("\n")
    print(x$sec.var)
    options(digits = 7)
    cat("\n")
    cat("\n")
    cat("Estimation summary for the regression coefficient :\n")
    cat("\n")
    print(round(x$reg, 3))
  }

}

plot.gbp <- function(x, sort = TRUE, ...) {
  y <- x$sample.mean
  se <- x$se

  if (any(is.na(x$prior.mean))) {
    pr.m <- x$prior.mean.hat
  } else {
    pr.m <- x$prior.mean
  }
  po.m <- x$post.mean
  po.sd <- x$post.sd
  po.low <- x$post.intv.low
  po.upp <- x$post.intv.upp

  if (sort == TRUE) {
    temp.data <- as.data.frame(cbind(y, se, pr.m, po.m, po.sd, po.low, po.upp))
    if (x$model == "gr") {
      temp.data <- temp.data[order(temp.data$se, decreasing = TRUE), ]
    } else {
      temp.data <- temp.data[order(temp.data$se), ]
    }
    y <- temp.data$y
    se <- temp.data$se
    pr.m <- temp.data[, 3]
    po.m <- temp.data[, 4]
    po.sd <- temp.data[, 5]
    po.low <- temp.data[, 6]
    po.upp <- temp.data[, 7]
  }

  index <- 1 : length(se)
  ylim.low <- ifelse(min(po.low, y) >= 0, 0.8 * min(po.low, y), 1.2 * min(po.low, y))
  ylim.upp <- ifelse(max(po.upp, y) >= 0, 1.2 * max(po.upp, y), 0.8 * max(po.upp, y))
  
  par(fig = c(0.25, 1, 0.5, 1), xaxs = "r", yaxs = "r", mai = c(0.5, 0.5, 0.5, 0.3), las = 1, ps = 13)

  if (x$model != "gr") {
    se <- sqrt(y * (1 - y) / se)
  }

  sqrtV <- se
  sdlens <- sqrtV / max(sqrtV)
  postlens <- po.sd / max(sqrtV)
  xmin <- min(c(y, po.m, pr.m))
  xmax <- max(c(y, po.m, pr.m))
  
  sunflowerplot(rep(4, length(y)) ~ y, ylim = c(-1, 5), xlim = c(xmin - abs(xmin) * 0.1, 
                xmax + abs(xmax) * 0.1), yaxt = "n", col.lab = "white", main = "Shrinkage plot", pch=1,cex=1)
  
  if (length(unique(pr.m)) == 1) {
    abline(v = pr.m, col = 4)
  } else {
    points(pr.m, rep(0,length(pr.m)), col = 4, pch = "-", cex = 2)
  }
  
  sunflowerplot(rep(0, length(y)) ~ po.m, add = TRUE,col="red",cex=1,pch=16)
  abline(h = 4)
  abline(h = 0)
 # axis(2, c(0, 4), c(expression(hat(theta)), expression(bar(y))), cex.axis = 1.1)
  sapply(1 : length(y), function(i) {
    lines(c(y[i], po.m[i]), c(4, 0))
    lines(c(y[i], y[i] + sdlens[i] * sd(y) * 0.4), c(4, 4 + sdlens[i]), col = "darkviolet")
    ##posterior variance lines
    lines(c(po.m[i] - postlens[i] * sd(y) * 0.4, po.m[i]), c(0 - postlens[i], 0), col = "darkgreen")
    xcord <- ((4 * po.m[i] / (y[i] - po.m[i]) - 4 * po.m / (y - po.m)) / 
              (4 / (y[i] - po.m[i]) - 4 / (y - po.m)))
    ycord <- 4 / (y - po.m) * xcord - 4 / (y - po.m) * po.m
    coords <- subset(cbind(xcord, ycord), ycord > 0 & ycord < 4)
    points(coords, pch=0)
  })

  par(fig = c(0, 1, 0, 0.5), xaxs = "r", yaxs = "r", mai = c(0.8, 0.7, 0.5, 0.3), las = 1, 
      ps = 13, new = TRUE)  


  if (sort == TRUE) {
    if (x$model == "gr") {
      xl <- c("Groups sorted by the descending order of se")
    } else {
      xl <- c("Groups sorted by the ascending order of n")
    }
  } else {
    xl <- c("Groups in the order of data input")
  }

  plot(index, po.m, ylim = c(ylim.low, ylim.upp), xlab = xl, ylab = "",
       main = paste(100 * x$confidence.lvl, "% Interval plot"), 
       col = 2, pch = 19)
  sapply(1 : length(y), function(j) {
    lines(rep(index[j], 2), c(po.low[j], po.upp[j]), lwd = 0.5)
  })
  points(index, po.low, cex = 1.5, pch = "-")
  points(index, po.upp, cex = 1.5, pch = "-")
  points(index, y)
  if (length(unique(pr.m)) == 1) {
    abline(h = pr.m, col = 4)
  } else {
    points(index, pr.m, col = 4, pch = "-", cex = 2)
  }

  ## legend
  if (x$model == "gr") {
    se.or.n <- "Standard error"
  } else {
    se.or.n <- "Group size (n)"
  }

  if (x$model == "br" & length(x$weight) != 1) { 
    par(fig = c(0, 0.35, 0.5, 1), xaxs = "r", yaxs = "r", mai = c(0.4, 0.1, 0.5, 0), las = 1, ps = 9,
        oma = c(0, 0, 0, 0), new = TRUE)  
    plot(1, type="n", axes=F, xlab="", ylab="")
    legend("topleft", pch = c(19, 1, NA, NA, NA,0), col = c(2, 1, 4,"darkviolet", "darkgreen",1), 
           lwd = c(NA, NA, 2, 2, 2), 
           c("Posterior mean \nof random effect", "Sample mean", 
             "Posterior mean of\nexpected random effect", se.or.n, 
             "Posterior sd \nof random effect", "Crossover"),
           seg.len = 0.5, bty = "n",xpd = TRUE)
  } else {
    par(fig = c(0, 0.35, 0.5, 1), xaxs = "r", yaxs = "r", mai = c(0.4, 0.1, 0.5, 0), las = 1, ps = 13,
        oma = c(0, 0, 0, 0), new = TRUE)  
    plot(1, type="n", axes=F, xlab="", ylab="")
    legend("topleft", pch = c(19, 1, NA, NA, NA,0), col = c(2, 1, 4,"darkviolet", "darkgreen",1), 
           lwd = c(NA, NA, 2, 2, 2), 
           c("Posterior mean", "Sample mean", "Prior mean", se.or.n, "Posterior sd", "Crossover"),
           seg.len = 0.5, bty = "n",xpd = TRUE)
  }


}
