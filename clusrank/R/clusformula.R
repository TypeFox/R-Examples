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
##   along with the R package clusrank. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


#'Wilcoxon Rank Sum Test for Clustered Data
#'
#'This is the sum rank test to compare the means of scores from
#'two groups for clustered data. The cluster size can be either
#'identitical or variable. Effect of stratification on the test
#'is also adjusted for if in presence.
#'
#'@author Yujing Jiang
#'@references
#'Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2003)
#' \emph{Incorporation of Clustering Effects for the Wilcoxon Rank
#' Sum Test: A Large-Sample Approach.} Biometrics, \bold{59}, 1089-1098.
#' @describeIn cluswilcox.test formula interface for sum rank test.
#' @importFrom  stats na.omit terms complete.cases model.extract
#' @export

cluswilcox.test.formula <- function(formula, data = NULL,
                                    subset = NULL, na.action = na.omit,
                                    alternative = c("two.sided", "less", "greater"),
                                    mu = 0,
                                    group.x = NULL, permutation = FALSE,
                                    n.rep = 500, ...) {
  ## This function is only used for rank sum test.
  ## Mainly to process data.
  ## Inputs:
  ##  formula: a formula of the form "lhs ~ rhs",
  ##    where the lhs is a numeric variable giving the
  ##    data values, the observed score and the rhs
  ##    is of the form with special terms:
  ##    cluster(x1) + group(x2) + stratum(x3), where
  ##    x1, x2, x3 are the corresponding variables.
  ##  data: an optional matrix or dataframe of data
  ##    used in the formula.
  ##  subset: an optional vector specifying
  ##   a subset of observations to be used.
  ##  na.action: a function which indicates what should happen when
  ##    the data contain NAs. Defaults to getOption("na.action").
  ##
  ##
  ##

  METHOD <- "Wilcoxon rank sum test for clutered data"

  alternative <- match.arg(alternative)
  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")
  Call <- match.call()

  if(!missing(data)) {
    DNAME <- paste("from", Call$data)
  }


  indx <- match(c("formula", "data", "subset", "na.action"),
                names(Call), nomatch = 0)
  if(indx[1] == 0)
    stop(" A formula argument is required for clustered rank sum test")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  special <- c("stratum", "cluster", "group")
  temp$formula <- if(missing(data))
    terms(formula, special)
  else terms(formula, special, data = data)
  x.name <- rownames(attr(temp$formula, "factors"))[1]
  DNAME <- paste0(paste(x.name ,DNAME), ",")




  mf <- eval(temp, parent.frame())

  if(nrow(mf) == 0)
    stop("No (non-missing) observations")

  Terms <- terms(mf)
  extraArgs <- list(...)

  x <- model.extract(mf, "response")
  if(is.vector(x)) {
    data.n <- length(x)
  } else {
    data.n <- nrow(x)
  }


  strats <- attr(Terms, "specials")$stratum

  if(length(strats)) {
    stemp <- untangle.specials(Terms, "stratum", 1)
    strats.name <- gsub("[\\(\\)]", "",
                         regmatches(stemp$vars,
                                    gregexpr("\\(.*?\\)", stemp$vars))[[1]])
        DNAME <- paste0(  DNAME, " stratum: ", strats.name, ",")

    if(length(stemp$vars) == 1) {
      strats.keep <- mf[[stemp$vars]]
    } else {
      stop("more than one variable are set as the stratum id")
    }
    strats.uniq <- unique(strats.keep)
    strats.uniq.l <- length(strats.uniq)

    if(is.character(strats.uniq)) {
      strats <- recoderFunc(strats.keep, strats.uniq, c(1 : strats.uniq.l))
    } else {
      if(!is.numeric(strats.uniq)) {
        stop("stratum id should be numeric or character")
      }
      strats <- strats.keep
    }
  } else {
    strats <- rep(1, data.n)
  }

  cluster <- attr(Terms, "specials")$cluster
  if(length(cluster)) {
    ctemp <- untangle.specials(Terms, "cluster", 1)
    cluster.name <- gsub("[\\(\\)]", "",
                         regmatches(ctemp$vars,
                                    gregexpr("\\(.*?\\)", ctemp$vars))[[1]])
    DNAME <- paste0(DNAME, " cluster: ", cluster.name, ",")


    if(length(ctemp$vars) == 1) {
      cluster.keep <- mf[[ctemp$vars]]
    } else {
      stop("more than one variable are set as the cluster id")
    }
    cluster.uniq <- unique(cluster.keep)
    cluster.uniq.l <- length(cluster.uniq)

    if(is.character(cluster.uniq)) {
      cluster <- recoderFunc(cluster.keep, cluster.uniq, c(1 : cluster.uniq.l))
    } else {
      if(!is.numeric(cluster.uniq)) {
        stop("cluster id should be numeric or character")
      }
      cluster <- cluster.keep
    }
  } else {
    cluster <- c(1 : data.n)
  }

  group <- attr(Terms, "specials")$group
  if(length(group)) {
    gtemp <- untangle.specials(Terms, "group", 1)
    group.name <- gsub("[\\(\\)]", "",
                         regmatches(gtemp$vars,
                                    gregexpr("\\(.*?\\)", gtemp$vars))[[1]])
        DNAME <- paste0(DNAME, " group: ", group.name, ",")


    if(length(gtemp$vars) == 1) {
      group.keep <- mf[[gtemp$vars]]
    } else {
      stop("more than one variable are set as the group id")
    }
    group.uniq <- unique(group.keep)
    group.uniq.l <- length(group.uniq)

    if(!is.character(group.uniq) && !is.numeric(group.uniq)) {
      stop("group id has to be numeric or character")
    }

    if(group.uniq.l > 2L) {
      stop("can only handle 2 groups in rank sum test")
    }
    if(!is.null(group.x)) {
      if(! group.x %in% group.uniq) {
        stop("label of group x is not valid")
      }
    } else {
      group.x <- min(group.uniq)
    }
    group.y <- group.uniq[group.uniq != group.x]
    group <- recoderFunc(group.keep, c(group.x, group.y), c(1, 2))
  } else {
    stop("no group id for rank sum test")
  }

  OK <- complete.cases(x)
  finite.x <- is.finite(x)
  x <- x[OK && finite.x]
  group <- group[OK && finite.x]
  stratum <- group[OK && finite.x]
  cluster <- cluster[OK && finite.x]
  mu.vec <- (group == 1)* mu
  x <- x - mu.vec

  if(permutation == FALSE) {
    return(cluswilcox.test.ranksum(x,  cluster,
                                   group, strats,
                                   alternative,
                                   mu,
                                   DNAME, METHOD))
  } else {
    return(cluswilcox.test.ranksum.permutation(x,
                                               cluster,
                                               group,
                                               strats,
                                               alternative,
                                               mu,
                                               n.rep,
                                               DNAME,
                                               METHOD))
  }



}
