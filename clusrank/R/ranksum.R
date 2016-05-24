################################################################################
##
##   R package clusrank by Mei-Ling Ting Lee, Jun Yan, and Yujing Jiang
##   Copyright (C) 2015
##
##   This file is part of the R package clusrank.
##
##   The R package clusrank is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package clusrank is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##    with the R package clusrank. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


#'Wilcoxon Rank Sum Test for Clustered Data
#'
#'This is the rank sum test for clustered data.
#' The cluster size can be either
#'identitical or variable. Effect of stratification on the test
#'is also adjusted for if in presence.
#'@param x a numeric vector of data values.
#'@param cluster a numeric vector indicating the cluster ids of the scores.
#'If not specified, each score has its own id, i.e., there is no
#'cluster in the data.
#'@param group a numeric vector indicating group id.
#'@param strats a numeric vector indicatg the strats ids of
#'the scores. The default is that the data has no stratum.
#'@param alternative a character string specifying the
#' alternative hypothesis, must be one of "two.sided" (default),
#'  "greater" or "less". You can specify just the initial letter.
#'@param mu null value of the hypothesis
#'@param DNAME a character string, inheritated from
#' \code{cluswilcox.test.model}, for result output.
#' @param METHOD a character string, inheritated from
#' \code{cluswilcox.test.model}, for result output.

#'@return  a list with the following components
#'\item{rstatistic}{Clustered Wilcoxon ranksum statistic.}
#'\item{erstatistics}{Expected value clustered Wilcoxon ranksum statistic.}
#'\item{vrstatistics}{Variance of clustered Wilcoxon ranksum statistic.}
#'\item{statistics}{the value of the test statistic.}
#'\item{p.value}{the p-value for the test}
#'\item{data.name}{a character string giving the names of the data.}
#'\item{method}{the name of the method}
#'\item{balance}{a logical, indicating if the data is balanced.}
#'@importFrom  stats pnorm
#'@examples
#'data(crd)
#'cluswilcox.test(z ~ group(group) + cluster(id), data = crd)
#'data(crdStr)
#'cluswilcox.test(z ~ group(group) + cluster(id) + stratum(stratum), data = crdStr)
#'
#'@author Yujing Jiang
#'@references
#'Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2003)
#' \emph{Incorporation of Clustering Effects for the Wilcoxon Rank
#' Sum Test: A Large-Sample Approach.} Biometrics, \bold{59}, 1089-1098.



cluswilcox.test.ranksum <-
  function(x, cluster, group, strats,
           alternative, mu, DNAME = NULL, METHOD = NULL) {
    # Incoporating clustering effects for the WilcoxonRank Sum Test
    # for stratified balanced or unbalanced designs. In addition,
    # one can control for confounding by forming strata which are
    # defined based on cross-classfication of one or more categorical
    # covariates which are defined in terms of a single categorical
    # variable denoted by strata.
    #
    #
    # Requirement:
    #   An ASCII data file with 4 variables per record. The data file
    #   does not have to be sorted in ID order.
    #   The 4 variables need to be space selimited and arranged in the
    #   following order:
    #   1. id
    #   2. score
    #   3. group (X, Y) indicators : need to be 1 and 2.
    #   4. strats
    #
    # Args:
    #   1. name of data set.
    #   2. number of strats.
    #
    # Returns:
    #   1. Clustered Wilcoxon RankSum Statistic
    #   2. Expected Clustered Wilcoxon RankSum Statistic
    #   3. Variance of clustered Wilcoxon RankSum Statistic
    #   4. Z statistic for Clustered Wilcoxon RankSum Statistic
    #   5. P-value for CLustered Wilcoxon RankSum Z Statistic




  METHOD <- "Wilcoxon rank sum test for clutered data"

  data <- as.data.frame(cbind(x, cluster, group, strats))
  ## group is assumed to take value 1 and 2; could be made
  ## more general
  group.uniq <- unique(data$group) ## Record possible groups

  data <- data[with(data, order(cluster)), ]
  ## Reorder the observations by the order of id.
  xrank <- rank(data$x, na.last = NA)
  ## Compute the rank of x score

  data <- cbind(data, xrank)
  ## Add z score to the dataframe crd1

### calculate ranksum within each cluster within each strats

  cluster.size <- table(data$cluster)
  ## cluster.size is the cluster size
  sumrank <- c(by(data$xrank, data$cluster, sum))
  ## Compute rank sum within each cluster

  strats <- c(by(data$strats, data$cluster, mean))
  ## Compute strats  of each cluster

  group <- c(by(data$group, data$cluster, mean))
  ## Compute group  of each cluster


  #count number of subunits within cluster size group within each strats

  n.csize.strats <- as.data.frame(table(cluster.size, strats))
  ## Count the objects for each group and each strats
  n.csize.xy <- as.data.frame(table(cluster.size, strats, group))
  n.csize.x <- subset(n.csize.xy, group == 1)
  ## Take out the count for group 1
  n.csize.y <- subset(n.csize.xy, group == 2)

  ## Ngv is the number of clusters which have the same
  ## subunit, controlled by stratum.
  colnames(n.csize.strats)[3] <- "Ngv"

  ## mgv is number of clusters which have the same
  ## subunit in group x, ontrolled by stratum.
  colnames(n.csize.x)[4] <- "mgv"

  psumrnk <- cbind(cluster.size, strats, group, sumrank)

  str.uniq <- sort(unique(strats))
  n.str.uniq <- length(str.uniq)
  cluser.size.uniq <- unique(cluster.size)
  n.cluser.size.uniq <- length(cluser.size.uniq)

  if(n.cluser.size.uniq > 1L) {
    balance <- FALSE
  } else {
    balance <- TRUE
  }

  psumrnk_int <- matrix(0,  n.str.uniq * cluser.size.uniq , 3)

  k <- 1
  for(i in 1:n.str.uniq ){
    for(j in 1:n.cluser.size.uniq){
      ind <- strats == str.uniq[i] & cluster.size == cluser.size.uniq[j]
      foo <- sum(psumrnk[ind, "sumrank"])
      if (foo) {
        psumrnk_int[k,] <- c(str.uniq[i], cluser.size.uniq[j], foo)
        k <- k + 1
      }
    }
  }
  colnames(psumrnk_int) <- c("strats", "cluster.size", "psumrank")
  psumrnk_int <- as.data.frame(psumrnk_int[psumrnk_int[, "strats"] != 0, , drop = FALSE])

  WC <- sum(psumrnk[psumrnk[,"group"] == 1,"sumrank"])

  expwc <- merge(merge(n.csize.x, psumrnk_int), n.csize.strats)
  ExpWc <- sum(expwc[ , "mgv"] * expwc[ , "psumrank"] / expwc[, "Ngv"])

  varwc <- merge(merge(merge(psumrnk, n.csize.xy), n.csize.strats), psumrnk_int)
  varwc_int <- (varwc[ ,"sumrank"] - varwc[ ,"psumrank"] / varwc[, "Ngv"])^2
  VarWc <- cbind(varwc, varwc_int)
  varwc_pre_final <- VarWc[,"Freq"] * (VarWc[, "Ngv"]- VarWc[, "Freq"]) / (VarWc[,"Ngv"] * (VarWc[,"Ngv"] -1)) * VarWc[,"varwc_int"]
  varwc_final <- sum(varwc_pre_final[VarWc[,"Ngv"] > 1])
  ## Drop the value where there is only 1 cluster in that strats,
  ## and that cluster contains only 1

  zc <- (WC - ExpWc)/sqrt(varwc_final)
  pval <- switch(alternative,
                  less = pnorm(abs(zc)),
                  greater = pnorm(abs(zc), lower.tail = FALSE),
                  two.sided = 2 * min(pnorm(abs(zc)),
                                      pnorm(abs(zc), lower.tail = FALSE)))

  names(WC) <- "Rank sum statistic"
  names(ExpWc) <- "Expected value of rank sum statistic"
  names(varwc_final) <- "Variance of rank sum statistic"
  names(zc) <- "Test statistic"
  names(mu) <- "difference in locations"
  result <- list(rstatistic = WC, erstatistic = ExpWc,
                 vrstatistic = varwc_final,
                 statistic = zc, p.value = pval,
                 alternative = alternative,
                 null.value = mu,
                 data.name = DNAME,
                 method = METHOD,
                 balance = balance)
  class(result) <- "ctest"
  result
}
