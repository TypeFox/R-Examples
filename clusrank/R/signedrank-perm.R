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
##   along with the R package reda. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

#' The Wilcoxon Signed Rank Test for Clustered Data
#'
#' Performs one-sample Wilcoxon test on vectors of data using
#' large sample.
#'
#' @param x  numeric vector of data values. Non-finite (e.g.,
#' infinite or missing) values will be omitted.
#' @param cluster numeric or charater vector, the id of clusters.
#'  If not specified, each observation will
#' be assigned a distinct cluster, i.e., no cluster in the data.
#' @param alternative a character string specifying the
#' alternative hypothesis, must be one of "two.sided" (default),
#'  "greater" or "less". You can specify just the initial letter.
#'@param mu null value of the hypothesis
#'@param n.rep number of samples generated for permutation test.
#' @param DNAME a character string, inheritated from
#' \code{cluswilcox.test.formula}, for result output.
#' @param METHOD a character string, inheritated from
#' \code{cluswilcox.test.formula}, for result output.
#' @return  a list with class "\code{ctest}" containing
#' the following components:
#' \item{rstatistic}{the value of the signed rank statistic
#'  with a name describing it.}
#' \item{vrstatistic}{Variance of \code{rstatistic}.}
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{n}{Total number of observations.}
#' \item{cn}{Total number of clusters.}
#' \item{data.name}{a character string giving the names of the data.}
#' \item{method}{the type of test applied.}
#' \item{adjusted}{indicator of whether adjusted signed rank statistic is used.}
#' @note This function is able to deal with data with
#' clusterentitical or variable cluster size. When the data
#' is unbalanced, adjusted signed rank statistic is used.
#' Ties are dropped in the test.
#' @importFrom stats ecdf lm
#' @examples
#' data(crsd)
#' cluswilcox.test(z, cluster = id, data = crsd, permutation = TRUE)
#' data(crsdUnb)
#' cluswilcox.test(z, cluster = id, data = crsdUnb, permutation = TRUE)
#' @author Yujing Jiang
#' @references
#' Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006)
#' \emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#'  Clustered Data.} Biometrics, \bold{62}, 185-192.
cluswilcox.test.signedrank.permutation <-
  function(x, cluster, alternative, mu, n.rep = 500,
           DNAME = NULL , METHOD = NULL){

    #Calculate number of observations per cluster

    names(mu) <- "location shift"

    data <- data.frame(x, cluster)
    cluster.size <- table(data$cluster)
    m <- length(cluster.size)
    n <- nrow(data)
    if (length(table(cluster.size)) != 1) {
      balance = FALSE
    } else {
      balance = TRUE
    }
    xrank <- rank(abs(data$x))
    data <- cbind(data, xrank)
    signrank <- ifelse(data$x > 0, 1, -1) * data$xrank
    data <- cbind(data, signrank)
    colnames(data)[4] <- "signrank"

    if(balance == TRUE){
      T_c <- sum(data$signrank)
      sumrank <- c(by(data$signrank, data$cluster, sum))
      delta <- replicate(n.rep, sample(c(-1, 1), m, TRUE))
      T_c.sim <- colSums(delta * sumrank)
      ecdf.tc <- ecdf(T_c.sim)


      P_val <- switch(alternative,
                      less = ecdf.tc(abs(T_c)),
                      greater = 1 - ecdf.tc(abs(T_c)),
                      two.sided = 2 * min(ecdf.tc(abs(T_c)),
                                          1 - ecdf.tc(abs(T_c)),
                                          0.5))


      names(T_c) <- "rank statistic"
      names(n) <- "total number of observations"
      names(m) <- "total number of clusters"
      result <- list(rstatistic = T_c,
                     p.value = P_val,
                     n = n,  cn = m, permutation = TRUE,
                     alternative = alternative,
                     null.value = mu,
                     method = METHOD, data.name = DNAME)
      class(result) <- "ctest"
      return(result)
    } else {
      if (balance == FALSE) {
        sumclusterrank <- c(by(data$signrank, data$cluster, sum))
        sumsq <- sum(sumclusterrank ^ 2)
        meansumrank <- sumclusterrank / cluster.size
        sumsqi <- sum(cluster.size ^ 2)

        # calculate intraclass correlation between signed ranks within the same cluster
        data$cluster.f <- as.factor(data$cluster)
        mod <- lm(signrank ~ cluster.f, data, y = TRUE)
        errordf <- mod$df.residual
        errorss <- sum(mod$residuals ^ 2)
        modeldf <- n - errordf - 1
        modelss <- sum((mod$y - mean(mod$y)) ^ 2) - errorss
        sumi <- n

        m0 <- (sumi - (sumsqi / sumi)) / (m - 1)
        totalss <- errorss + modelss
        totaldf <- errordf + modeldf
        wthms <- errorss / errordf
        betms <- modelss / modeldf
        vars <- totalss / totaldf
        s2b <- (betms - wthms) / m0
        s2w <- wthms
        rosglm <- s2b / (s2b + s2w)
        if (rosglm < 0) {
          rosglm = 0
        }
        ros <- rosglm
        roscor <- ros * (1 + (1 - ros ^ 2) / (m - 2.5))
        if (roscor > 1) {
          roscor = 1
        }
        wi <- cluster.size / (vars * (1 + (cluster.size - 1) * roscor))
        T_c <- sum(meansumrank * wi)

        delta <- replicate(n.rep, sample(c(-1, 1), m, TRUE))
        T_c.sim <- colSums(delta * c(meansumrank) * c(wi))

        ecdf.tc <- ecdf(T_c.sim)


        P_val <- switch(alternative,
                        less = ecdf.tc(abs(T_c)),
                        greater = 1 - ecdf.tc(abs(T_c)),
                        two.sided = 2 * min(ecdf.tc(abs(T_c)),
                                            1 - ecdf.tc(abs(T_c)),
                                            0.5))


        names(T_c) <- "adjusted rank statistic"
        names(n) <- "total number of observations"
        names(m) <- "total number of clusters"
        result <- list(rstatistic = T_c,
                       p.value = P_val,
                       n = n,  cn = m, permutation = TRUE,
                       alternative = alternative,
                       method = METHOD, data.name = DNAME,
                       null.value = mu,
                       balance = balance)
        class(result) <- "ctest"
        return(result)
      }
    }



  }
