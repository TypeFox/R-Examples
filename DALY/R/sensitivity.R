## sensitivity analysis for DALY estimate

sensitivity <-
function(x, method = c("src", "pcc"), rank = FALSE, mapped = TRUE){
  ## evaluate inputs
  method <- match.arg(method)
  if (!is.logical(rank) && !is.na(rank))
    stop("'rank' must be a logical value", call. = FALSE)
  if (!is.logical(rank) && !is.na(rank))
    stop("'rank' must be a logical value", call. = FALSE)

  ## aggregate DALYs
  y <- aggregate(x, by = "total")
  daly <- y$DALY

  ## max number of parameters
  ## n_outcomes * 8
  n_outcomes <- length(x) - 2

  ## merge unique columns
  listNames <- c("inc", "trt", "ons", "dur", "DWt", "DWn", "mrt", "dth")
  data <- numeric(length(daly))
  for (i in seq(n_outcomes)){
    for (j in seq(8)){
      data <-
        cbind(data,
              as_column(x[[i]]$input[[j]], i, listNames[j]))
    }
  }

  ## remove fixed values
  fixed <- apply(data, 2, var) == 0
  data <- data[, !fixed]

  ## convert values to ranks
  if (rank){
    daly <- rank(daly)
    data <- apply(data, 2, rank)
  }

  if (method == "src"){
    ## scale variables
    data <- as.data.frame(apply(data, 2, scale))
    if (!mapped) daly <- scale(daly)

    ## linear regression
    out <- summary(lm(daly ~ ., data = data))

  } else {
    ## partial correlation coefficients
    data <- as.matrix(data)
    out <- matrix(ncol = 2, nrow = ncol(data))
    colnames(out) <- c("rho", "p")
    rownames(out) <- colnames(data)

    for (i in seq(ncol(data))){
      lm_y <- lm(daly ~ data[, -i])      # regress y to other x's
      lm_x <- lm(data[, i] ~ data[, -i]) # regress x to other x's
      out[i, ] <-
        unlist(cor.test(lm_y$residuals, lm_x$residuals)[4:3],
               use.names = FALSE)
    }
  }

  ## define S3 class
  out_list <-
    list(method = list(method = method, rank = rank, mapped = mapped),
         out = out)
  class(out_list) <- "DALY_sensitivity"

  ## return output
  return(out_list)
}

## PLOT method
plot.DALY_sensitivity <-
function(x, alpha = 0.05, main = "Sensitivity analysis",
         show_values = FALSE, value_digits = 3, value_cex = 0.6, ...){
  ## standardized (rank) regression coefficients
  if (x$method$method == "src"){
    ## sort estimates
    cf <- coef(x$out)[-1, ]
    signif <- cf[, 4] < alpha
    order <- order(abs(cf[signif, 1]))
    est <- cf[signif, 1][order]

    ## define 'xlab'
    xlab <-
      ifelse(x$method$mapped,
             ifelse(x$method$rank,
                    "change in rank per sd change",
                    "change in DALY per sd change"),
             ifelse(x$method$rank,
                    "standardized rank regression coefficient",
                    "standardized regression coefficient"))
  }

  ## partial (rank) correlation coefficientts
  if (x$method$method == "pcc"){
    ## sort estimates
    signif <- x$out[, 2] < alpha
    order <- order(abs(x$out[signif, 1]))
    est <- x$out[signif, 1][order]

    ## define 'xlab'
    xlab <-
      ifelse(x$method$rank,
             "partial rank correlation coefficient",
             "partial correlation coefficient")
  }

  ## plot estimates
  par(mar = c(4, 6, 2, 2) + .5)
  ifelse(
    show_values,
    xlim <- c(0, max(est) * 1.1),
    xlim <- c(0, max(est)))
  bp <-
    barplot.default(est, horiz = T, las = 1,
                    main = main, xlab = xlab, xlim = xlim, ...)
  if (show_values)
    text(est, bp, round(est, value_digits),
         pos = 4, cex = value_cex, xpd = NA)
}

## PRINT method
print.DALY_sensitivity <-
function(x, digits = 3, signif_stars = getOption("show.signif.stars"), ...){
  ## standardized (rank) regression coefficients
  if (x$method$method == "src"){
    ## sort estimates
    cf <- coef(x$out)[-1, ]
    order <- order(abs(cf[, 1]), decreasing = TRUE)
    est <- cf[order, 1:4]
    if (signif_stars){
      signif <-
        symnum(est[, 4], corr = FALSE, na = FALSE, 
               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
               symbols = c("***", "**", "*", ".", " "))
      est <-
        data.frame(formatC(est[, 1:2], digits = digits,
                   format = ifelse(x$method$mapped, "f", "g")),
                   formatC(est[, 3], digits = digits, format = "f"),
                   formatC(est[, 4], digits = digits, format = "g"),
                   format(signif))
      colnames(est) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
    }

    ## define title
    title <-
      ifelse(x$method$mapped,
             ifelse(x$method$rank,
                    "Mapped standardized rank regression coefficients:",
                    "Mapped standardized regression coefficients:"),
             ifelse(x$method$rank,
                    "Standardized rank regression coefficients:",
                    "Standardized regression coefficients:"))
  }

  ## partial (rank) correlation coefficientts
  if (x$method$method == "pcc"){
    ## sort estimates
    order <- order(abs(x$out[, 1]), decreasing = TRUE)
    est <- x$out[order, ]
    if (signif_stars){
      signif <-
        symnum(est[, 2], corr = FALSE, na = FALSE, 
               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
               symbols = c("***", "**", "*", ".", " "))
      est <-
        data.frame(formatC(est, digits = digits, format = "g"),
                   format(signif))
      colnames(est) <- c("Rho", "Pr(>|0|)", "")
    }

    ## define title
    title <-
      ifelse(x$method$rank,
             "Partial rank correlation coefficients:",
             "Partial correlation coefficients:")
  }

  ## print results
  cat(title, "\n")
  print(est, ...)
  if (signif_stars)
    cat("---\nSignif. codes:  ", attr(signif, "legend"), "\n", sep = "")
  if (x$method$method == "src")
    cat("\nAdjusted R-squared:",
        formatC(x$out$adj.r.squared, digits = digits),
        "\n")
}