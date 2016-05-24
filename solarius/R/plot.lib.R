#----------------------------------
# Pedigree plots
#----------------------------------

#' Plot the pedigree
#'
#' Plot the pedigree tree based on the ID-fields given in phenotypes.
#' Plotting is based on \code{kinship2} package.
#' 
#' The ID fields are extracted from the data via \code{matchIdNames} function.
#' The required fields are ID, FA, MO, SEX and FAMID.
#'
#' @param data A data.frame of phenotypes.
#' @param ped An integer or a character value, that indicates the pedigree.
#'
#' @seealso \code{\link{matchIdNames}}
#' @examples
#' data(dat30)
#' plotPed(dat30, 1)
#'
#' @export
plotPed <- function(data, ped)
{
  # checks
  stopifnot(!missing(data))
  stopifnot(class(data) == "data.frame")
  stopifnot(!missing(ped))

  # load required R package kinship2
  stopifnot(requireNamespace("kinship2", quietly = TRUE))
  
  ret <- NULL # R CMD check: no visible binding
  cmd <- "ret <- require(kinship2)"
  eval(parse(text = cmd))
  if(!ret) {
    stop("`kinship2` package is required for plotting pedigrees.")
  }
      
  renames <- matchIdNames(names(data))
  data <- rename(data, renames)
  stopifnot(all(c("ID", "FA", "MO", "SEX", "FAMID") %in% names(data)))
  
  if(sum(with(data, is.na(FAMID)))) {
    warning("plotPed: rows with FAMID == NA are filtered out.")
    data <- subset(data, !is.na(FAMID))
  }
  
  peds <- unique(data$FAMID)
  
  # filter by `ped`
  ped <- switch(class(ped),
    "integer" = peds[ped],
    "numeric" = peds[ped],
    "character" = ped,
    stop("switch error"))
  stopifnot(ped %in% peds)  
  
  ### create `df`
  FAMID <- NULL # R CMD check: no visible binding

  df <- subset(data, FAMID == ped)
  if(nrow(df) == 0) {
    warning("plotPed: no families found with the given value of `ped`.")
    return(invisible()) 
  }
  
  ### filter `df`
  ids <- df$ID
  df <- within(df, {
    ind <- !(FA %in% ids)
    FA[ind] <- ""
    MO[ind] <- ""
    
    ind <- !(MO %in% ids)
    FA[ind] <- ""
    MO[ind] <- ""
  })

  # fix for `SEX
  df <- within(df, {
    if(class(SEX) == "character" & all(c("1", "2") %in% SEX)) { 
      SEX <- as.integer(SEX)
      SEX[is.na(SEX)] <- 3
    }
  })

  # make `pedigree` object
  ped <- with(df, 
    kinship2::pedigree(id = ID, dadid = FA, momid = MO, sex = SEX, famid = FAMID, missid = ""))  

  plot(ped[1])
}

#----------------------------------
# Association plots
#----------------------------------

#' Plot the association results
#'
#' Two Manhattan and quantile-quantile (QQ) plots 
#' are standard to explore the results from association studies.
#'
#' \code{plotManh} function produces the Manhattan plot based on \code{qqman} package.
#' The two red and blue lines, default in the original \code{manhattan} function 
#' of \code{qqman} package, are preserved.
#' An additional black dashed line is added, that depicts the significance level
#' according to Bonferroni multiple-test correction with \code{alpha} argument.
#'
#' \code{plotQQ} function produces the QQ plot based on the same \code{qqman} package.
#' 
#' @name plotManh
#' @rdname plotQQManh
#'
#' @param A An object of class \code{solarAssoc}.
#' @param alpha A numeric value from 0 to 1, the p-value cut-off 
#'    for Bonferroni multiple-test correction.
#'    The default value is 0.05.
#' @param main A character string, \code{main} argument (title) to a plot function.
#'    If the argument is missing, the title contains information about the model formula and the number of SNPs used for plotting.
#' @param ... additional argument to a plot function.
#'
#' @seealso \code{\link{solarAssocClass}}
#'
#' @examples
#' \dontrun{
#' data(dat50)
#'
#' assoc <- solarAssoc(trait ~ 1, phenodata, 
#'  snpdata = genodata, snpmap = snpdata, kinship = kin)
#'
#' plotManh(assoc) # equivalent to plot(assoc)
#'
#' plotQQ(assoc) # equivalent to plot(assoc, "qq")
#'
#' }
#'
#' @export
plotManh <- function(A, alpha = 0.05, main, ...)
{
  stopifnot(requireNamespace("qqman", quietly = TRUE))
  
  stopifnot(all(c("chr", "pSNP", "pos") %in% names(A$snpf)))
  
  df <- as.data.frame(A$snpf)

  ind <- with(df, !is.na(chr) & !is.na(pSNP) & !is.na(pos))
  num.na <- sum(!ind)
  num.snps <- sum(ind)
  
  stopifnot(num.snps > 0)
  
  df <- df[ind, ]
  chr <- NULL # R CMD check: no visible binding
  df <- mutate(df,
    chr = as.integer(chr))
  ### plot parameters
  if(missing(main)) {
    main1 <- paste(paste(A$traits, collapse = "+"), "~", paste(A$covlist, collapse = "+"))
    main2 <- paste("#SNPs:", num.snps) 
    if(num.na > 0) {
      main2 <- paste(main2, " (#SNPs missing: ", num.na, ")", sep = "")
    }
    
    main <- paste(main1, main2, sep = "\n")
  }
  
  ### plot    
  qqman::manhattan(df, chr = "chr", bp = "pos", p = "pSNP", snp = "SNP", 
    main = main, ...)
  abline(h = -log10(alpha / num.snps), col = "black", lty = 2)
}

#' @name plotManh
#' @rdname plotQQManh
#'
#' @param df An integer value, the degree of freedom.
#'  The default value is 1.
#'
#' @export
plotQQ <- function(A, df = 1, main, ...)
{
  stopifnot(requireNamespace("qqman", quietly = TRUE))
  
  # get p-values
  stopifnot("pSNP" %in% names(A$snpf))
  
  snpf <- as.data.frame(A$snpf)

  ind <- with(snpf, !is.na(chr) & !is.na(pSNP) & !is.na(pos))
  num.na <- sum(!ind)
  num.snps <- sum(ind)
  
  stopifnot(num.snps > 0)
  
  pvals <- snpf[ind, "pSNP"]
  
  # estimate lambda inflation-statistics
  lambda <- median(pvals, na.rm = TRUE) / qchisq(0.5, df)
  sub <- paste("lambda", round(lambda, 2))
  
  # plot par
  if(missing(main)) {
    main1 <- paste(paste(A$traits, collapse = "+"), "~", paste(A$covlist, collapse = "+"))
    main2 <- paste("#SNPs:", num.snps) 
    if(num.na > 0) {
      main2 <- paste(main2, " (#SNPs missing: ", num.na, ")", sep = "")
    }
    
    main <- paste(main1, main2, sep = "\n")
  }
  
  # plot
  qqman::qq(pvals, main = main, sub = sub, ...)
}

#----------------------------------
# Plot Kinship2
#----------------------------------

#' Plot the double kinship matrix
#'
#' The main function that calls \code{imageKinship2} or \code{histKinship2}
#' depending on value of \code{y} argument.
#'
#' \code{imageKinship2} function calls \code{image} function from \code{Matrix} package.
#'
#' \code{histKinship2} function plots a histogram based on \code{ggplot2} package.
#' 
#' @name plotKinship2
#' @rdname plotKinship2
#'
#' @param x A square matrix of double kinship coefficients.
#' @param y A character, the type of the plot.
#'    Possible values are \code{"image"} and \code{"hist"}. 
#'    The default value is \code{"image"}.
#'
#' @examples
#' # load `kin` kinship matrix from `dat50` data set
#' data(dat50)
#' kin2 <- 2* kin # double kinship matrix
#'
#' plotKinship2(kin2) # equivalent to `imageKinship2(kin2)`
#'
#' plotKinship2(kin2, "hist") # equivalent to `histKinship2(kin2)`
#'
#' @export
plotKinship2 <- function(x, y = c("image", "hist"))
{
 y <- match.arg(y)
 
  switch(y,
    "image" = imageKinship2(x),
    "hist" = histKinship2(x),
    stop())
}

#' @name imageKinship2
#' @rdname plotKinship2
#'
#' @export imageKinship2
imageKinship2 <- function(x)
{ 
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  
  Matrix::image(Matrix::Matrix(x))
}

#' @name histKinship2
#' @rdname plotKinship2
#'
#' @export histKinship2
histKinship2 <- function(x)
{
  #stopifnot(require(ggplot2))
  
  df <- data.frame(kin = as.vector(x))
  
  kin <- NULL # R CMD check: no visible binding
  ggplot(df, aes(kin)) + geom_histogram() + 
    labs(x = "") +
    theme_bw()
}

#----------------------
# Plot polygenic model
#----------------------

#' Plot the residuals of a polygenic model
#'
#' Plot the residuals on scatter or quantile-quantile plots.
#' 
#' \code{plotRes} function makes a scatter plot of fitted values vs. residuals.
#' Note that the residuals returned by SOLAR include both random effects,
#' i.e. house-hold, genetic and residuals itself.
#'
#' \code{plotResQQ} function plots  quantile-quantile (QQ) plot of the residuals.
#' 
#' @name plotRes
#' @rdname plotRes
#'
#' @param x An object of class \code{solarPolygenic}.
#' @param labels A logical value for \code{plotRes} function, indicating if the labels of IDs 
#'    (which residuals are outside the 3 * sd interval) are to be plotted.
#'    A logical value for \code{plotResQQ} function, 
#'    indicating if the samples (their IDs) outside the confidence intervals are to be plotted.
#' @param text.size An integer, the text size of labels.
#' @param ... additional arguments.
#'
#' @seealso \code{\link{solarPolygenicClass}}
#' @examples
#'
#' \dontrun{
#' ### basic (univariate) polygenic model
#' mod <- solarPolygenic(trait1 ~ age + sex, dat30) 
#'
#' plotRes(mod)
#'
#' plotResQQ(mod)
#'
#' }
#'
#' @export
plotRes <- function(x, 
  labels = FALSE, text.size = 4, ...)
{
  #stopifnot(require(ggplot2))
  
  # var
  r <- residuals(x)
  yh <- residuals(x, trait = TRUE)
  labs <- names(r)
  
  # sd  
  r.sd <- sd(r, na.rm = TRUE)
  ind <- which(!(r > 3*r.sd | r < -3*r.sd)) # residuals inside [-3 sd; 3 sd]
  labs[ind] <- ""
  
  # data for plotting
  ord <- order(yh)
  df <- data.frame(ord = ord, yh = yh[ord], r = r[ord], label = labs[ord])
  
  # plot
  p <- ggplot(df, aes(x = ord, y = r)) + geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -3*r.sd, linetype = "dashed") + 
    geom_hline(yintercept = 3*r.sd, linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE) + #, se = TRUE, level = conf)
    labs(title = "Residuals",  
      x = "Trait order", y = "Residuals")
  
  label <- NULL # R CMD check: no visible binding
  if(labels) {
    p <- p + geom_text( aes(label = label), size = text.size) # hjust = 0, vjust = 0
  }
  
  # print
  if(labels) {
    labs <- with(df, label[label != ""])
    if(length(labs) > 0) {
      cat(" * Sample(s) outside the 3*sd interval: ", 
        paste(labs, collapse = ", "), "\n", sep = "")
    } else {
      cat(" * All sampes are within the 3*sd interval\n", sep = "")
    }
  }
  
  # return
  return(p)
}
 
 
#' @name plotRes
#' @rdname plotRes
#'
#' @param conf 
#'    A numeric value between 0 and 1, that represents the confidence boundary.
#' @param distribution 
#'    A character, name of distribution of the residuals.
#'    The default value is \code{"norm"}.
#' @param line.estimate
#'    Function for estimation of QQ-line.
#'
#' @export
plotResQQ <- function(x, distribution = "norm", ..., line.estimate = NULL, 
  conf = 0.90, 
  labels = FALSE, text.size = 4)
{
  # source: http://stackoverflow.com/a/27191036/551589
  # alternatives: 
  #  -- qqnorm(res); qqline(res)
  #  -- cars::qqPlot(res, id.method = "identify")

  #stopifnot(require(ggplot2))

  z <- ord.r <- lower <- upper <- label <- NULL # R CMD check: no visible binding
  
  ### var
  r <- residuals(x)
  labs <- names(r)
  
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))

  r <- na.omit(r)
  ord <- order(r)
  n <- length(r)
  P <- ppoints(length(r))
  df <- data.frame(ord.r = r[ord], z = q.function(P, ...))

  if(is.null(line.estimate)) {
    Q.r <- quantile(df$ord.r, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.r)/diff(Q.z)
    coef <- c(Q.r[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.r ~ z))
  }

  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE

  df$label <- ifelse(df$ord.r > df$upper | df$ord.r < df$lower, labs[ord], "")

  p <- ggplot(df, aes(x=z, y=ord.r)) + geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2]) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    labs(title = "Q-Q plot", 
      x = paste0("Theoretical quantiles (", distribution, " distribution)"),
      y = "Sample quantiles")

  if(labels) {
    p <- p + geom_text( aes(label = label), size = text.size) # hjust = 0, vjust = 0
  }
  
  ### print
  if(labels) {
    labs <- with(df, label[label != ""])
    if(length(labs) > 0) {
      cat(" * Sample(s) outside the confidence interval (", conf, "): ", 
        paste(labs, collapse = ", "), "\n", sep = "")
    } else {
      cat(" * All sampes are within the confidence interval (", conf, ")\n", sep = "")
    }
  }
  
  return(p)
}
