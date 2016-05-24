# CPChisq.R -- comparative chi-square test for heterogeneous associations
#
# MS
# Created: July 3, 2015 generalized from two contingency tables in
#                       "hetero-chisq-steel.R"
# Modified: Feb 5, 2016. MS
#   Modified the method name
#   Added log.p argument to specify log p-value

cp.chisq.test <- function(x, method="chisq", log.p=FALSE)
{
  DNAME <- deparse(substitute(x))
  if(mode(x)!="list" || length(x)<2 )
  {
    stop("only accept list of 2 or more matrices as input!")
  }

  pooled <- x[[1]]
  for(i in 2:length(x)) {
    pooled <- pooled + x[[i]]
  }

  # Row sum of matrix pooled
  pooled.row.sum <- apply(pooled, 1, sum)

  # Column sum of matrix pooled
  pooled.col.sum <- apply(pooled, 2, sum)

  #Sum of matrix pooled
  pooled.sum <- sum(pooled.col.sum)

  pooled.chisq <- 0

  if(pooled.sum > 0) {
    #Expected value in pooled
    pooled.expected <- pooled.row.sum %*% t(pooled.col.sum) / pooled.sum

    #Chi-square value of pooled
    pooled.chisq <- sum((pooled - pooled.expected)^2 / pooled.expected, na.rm = TRUE)

    total.chisq <- 0

    for(i in seq_along(x)) {
      #Expected values in matrix i
      expected <- pooled.expected / pooled.sum * sum(x[[i]])

      #Chi-square value of matrix i
      chisq <- sum((x[[i]] - expected)^2 / expected, na.rm = TRUE)

      total.chisq <- total.chisq + chisq
    }

    hetero.chisq  <- total.chisq - pooled.chisq
    df <- (length(x) - 1) * (sum(pooled.row.sum != 0) - 1) *
      (sum(pooled.col.sum != 0) - 1)

  } else {

    hetero.chisq <- 0
    df <- 0

  }

  if(method=="default" || method=="chisq") {

    finalStat <- hetero.chisq
    names(finalStat) <- "statistic"

    finalDf <- df
    names(finalDf) <- "parameter"

    p.val <- pchisq(finalStat, df = finalDf, lower.tail=FALSE, log.p=log.p)
    names(p.val) <- "p.value"

    return(structure(list( statistic=finalStat, parameter=finalDf, p.value=p.val,
                           method = "Comparative chi-square heterogeneity test",
                           data.name= DNAME),
                     class = "htest"))

  } else if(method=="normalized" || method=="nchisq" ) {

    finalStat <- ifelse(df > 0, ( hetero.chisq - df) / sqrt( 2 * df ), 0)
    names(finalStat) <- "statistic"

    finalDf <- df
    names(finalDf) <- "parameter"

    p.val <- ifelse(df > 0, pnorm( finalStat, lower.tail=FALSE, log.p=log.p), 1)
    names(p.val) <- "p.value"

    return(structure(list(statistic = finalStat, parameter = finalDf, p.value = p.val,
                          method = "Nomalized comparative chi-square heterogeneity test",
                          data.name= DNAME),
                     class = "htest"))

  } else {
    stop("method can only be \"default\", \"normalized\", or \"exact\".\n")
  }
}
