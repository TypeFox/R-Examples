##==============================================================================
## Author Y. Deville
##
## This function is for internal use only and should not appear in the NAMESPACE
## file.
##
## It groups expected counts in a one-dimensional chi-square contingency table
## in order to reach the fatidic number of 5. The counts are grouped from left
## waiting to reach a cumulated number >= 5 in each group. The last
## group can still have count less than 5 in which case it has to be grouped
## with the group on its left. 
##
##==============================================================================

groupchi2 <- function(f.theo) {

  nf <- length(f.theo)

  group <- rep(NA, nf)
  g <- 1
  Cum <- 0

  for (i in 1:nf) {
    Cum <- Cum + f.theo[i]
    if ( Cum >= 5 ) {
      group[i] <- g
      Cum <- 0
      g <- g + 1
    } else group[i] <- g
  }

  G <- (1:nf)[group == group[nf]]
  
  if (sum(f.theo[G]) < 5) group[G] <- group[G] - 1

  names(group) <- 1:nf
  invisible(group)

}


##==============================================================================
## Author Y. Deville
##
## Perfoms simple analyses for counts in blocks. These counts are hoped to
## follow an homogenous Poisson Process (HPP) and form what is sometimes called
## a partially observed HPP.
##
##==============================================================================

barplotRenouv <- function(data,
                          blockname = colnames(data)[1],
                          varname = colnames(data)[2],
                          threshold = quantile(data[ , varname], 0.2),
                          na.block = NULL,
                          plot = TRUE,
                          main = NULL,
                          xlab = NULL,
                          ylab = NULL,
                          mono  = FALSE,
                          prob.theo = 0.999,
                          ...) {
  
  if ( is(data[ , blockname], "factor") ) {
    warning("'blockname' of class factor converted to 'integer'")
    b <- as.integer(format(data[ , blockname]))
  } else b <- data[ , blockname]

  ## 
  block.start <- min(b)
  block.end   <- max(b)
  block <- factor(b, levels = block.start:block.end)

  tab1  <- tapply(data[ , varname], block, function(x) sum(x>threshold))
  tab1[is.na(tab1)] <- 0
  
  ## if some blocks are declared NA, they are eliminated.
  if (length(na.block)) {
    na.ind <- names(tab1) %in% as.character(na.block)
    if (all(na.ind)) stop("all observations lie in NA blocks")
    cat("number of obs. in NA blocks:", sum(na.ind),"\n")
    tab1 <- tab1[!na.ind]
  }

  ## Poisson parameter and values 
  mu   <- mean(tab1)
  xmax <- max(c(max(tab1, na.rm = TRUE),
                qpois(p = prob.theo, lambda = mu)))
              
  ## Caution: non-observed levels give NA values
  x <- 0:xmax
  nx <- length(x)
  tabf  <- factor(tab1, levels = x)
  f.obs <- table(tabf)

  ## Compute the theoretical frequencies
  ## Caution: the last class is "xmax and over"
  tot <- sum(f.obs)
  p <- dpois(x, lambda = mu)
  p[nx] <- ppois(xmax-1, lambda = mu, lower.tail = FALSE)
  f.theo <- p*tot
  
  group <- groupchi2(f.theo)

  f.all <- cbind(f.obs, f.theo, group)
  
  m <- group[length(group)]
  names(m) <- NULL
  colnames(f.all) <- c("obs.", "theo.", "group")
  
  ## chi-square test the df is m - 2 since one parameter is
  ## estimated from the data
  
  f.theo.group <- tapply(f.theo, as.factor(group), sum)
  f.obs.group  <- tapply(f.obs, as.factor(group), sum)
  
  chisq.stat <- sum( (f.theo.group-f.obs.group)^2/ f.theo.group ) +
    (1- sum(p))*tot
  chisq.df <- m - 2
  chisq.pvalue <- 1 - pchisq(chisq.stat, df = chisq.df)
  
  chisq.test <-
    list(statistic = chisq.stat, df = chisq.df, p.value = chisq.pvalue)
  
  cat("\nGoodness-of-fit test (Poisson). Stat =", chisq.test$stat,
      "df =", chisq.test$df, "p-value =", chisq.test$p.value,"\n")
  
  ## (over) dispersion test
  ID <- var(tab1, na.rm = TRUE) / mean(tab1, na.rm = TRUE)

  disp.stat <- (length(tab1)-1)*ID
  disp.df <- length(tab1) - 1
  disp.pvalue <- 1-pchisq(disp.stat, length(tab1)-1)

  disp.test <-
    list(statistic = disp.stat, df = disp.df, p.value = disp.pvalue)

  ## group the two tests
  tests <- matrix(c(disp.stat, chisq.stat, disp.df, chisq.df, disp.pvalue, chisq.pvalue),
                  nrow = 2, ncol = 3, byrow = FALSE)
  
  rownames(tests) <- c("disp", "chi2")
  colnames(tests) <- c("statistic", "df", "p.value")
  
  cat("\nDispersion index", ID, "p-value", disp.pvalue,"\n")

  if (plot) {
    
    if (mono) cols <- c("white", "darkgray")
    else cols <- c("orange", "SteelBlue3")

    if (is.null(main)) main <- ""
    if (is.null(xlab)) xlab <-  "number of evts in block"
    if (is.null(ylab)) ylab <-  "number of blocks"
    
    barplot(t(f.all[ , 1:2]),
            beside = TRUE,
            legend = TRUE,
            col = cols,
            main = main,
            xlab = xlab,
            ylab = ylab,
            ...)
  }

  
  invisible(list(freq = f.all,
                 mu = mu,
                 overdispersion = ID,
                 disp.test = disp.test,
                 chisq.test = chisq.test,
                 tests = tests))
            
}


