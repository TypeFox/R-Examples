# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

#' Compares output observations from two or more groups
#'
#' Compares output observations from two or more groups.
#'
#' @param name Comparison name (useful when calling this function to perform
#' multiple comparisons).
#' @param ve_npcs Percentage (\code{0 < ve_npcs < 1}) of variance explained by
#' the \emph{q} principal components (i.e. number of dimensions) used in MANOVA,
#' or the number of principal components (\code{ve_npcs > 1}, must be integer).
#' Can be a vector, in which case the MANOVA test will be applied multiple
#' times, one per specified variance to explain / number of principal
#' components.
#' @param data A \emph{n} x \emph{m} matrix, where \emph{n} is the total number
#' of output observations (runs) and \emph{m} is the number of variables (i.e.
#' output length).
#' @param obs_lvls Levels or groups associated with each observation.
#' @param lim_npcs Limit number of principal components used for MANOVA to
#' minimum number of observations per group?
#' @param mnv_test The name of the test statistic to be used in MANOVA, as
#' described in \code{\link[stats]{summary.manova}}.
#'
#' @return  Object of class \code{cmpoutput} containing the following data:
#' \describe{
#'  \item{scores}{\emph{n} x \emph{n} matrix containing projections of
#'        output data in the principal components space. Rows correspond to
#'        observations, columns to principal components. }
#'  \item{obs_lvls}{Levels or groups associated with each observation.}
#'  \item{varexp}{Percentage of variance explained by each principal component.}
#'  \item{npcs}{Number of principal components specified in \code{ve_npcs} OR
#'        which explain the variance percentages given in \code{ve_npcs}.}
#'  \item{ve}{Percentage (between 0 and 1) of variance explained by the \emph{q}
#'        principal components (i.e. number of dimensions) used in MANOVA.}
#'  \item{name}{Comparison name (useful when calling this function to perform
#'        multiple comparisons).}
#'  \item{p.values}{\emph{P}-values for the performed statistical tests, namely:
#'    \describe{
#'      \item{manova}{List of \emph{p}-values for the MANOVA test for each
#'            number of principal component in \code{npcs}.}
#'      \item{parametric}{Vector of \emph{p}-values for the parametric test
#'            applied to groups along each principal component (\emph{t}-test
#'            for 2 groups, ANOVA for more than 2 groups).}
#'      \item{nonparametric}{Vector of \emph{p}-values for the non-parametric
#'            test applied to groups along each principal component
#'            (Mann-Whitney U test for 2 groups, Kruskal-Wallis test for more
#'            than 2 groups).}
#'      \item{parametric_adjusted}{Same as field \code{parametric}, but
#'            \emph{p}-values are adjusted using weighted Bonferroni procedure.
#'            Percentages of explained variance are used as weights.}
#'      \item{nonparametric_adjusted}{Same as field \code{nonparametric}, but
#'            \emph{p}-values are adjusted using weighted Bonferroni procedure.
#'            Percentages of explained variance are used as weights.}
#'    }
#'  }
#'  \item{tests}{
#'    \describe{
#'      \item{manova}{Objects returned by the \code{\link[stats]{manova}}
#'            function for each value specified in \code{ve_npcs}.}
#'      \item{parametric}{List of objects returned by applying
#'            \code{\link[stats]{t.test}} (two groups) or
#'            \code{\link[stats]{aov}} (more than two groups) to each principal
#'            component.}
#'      \item{nonparametric}{List of objects returned by applying
#'            \code{\link[stats]{wilcox.test}} (two groups) or
#'            \code{\link[stats]{kruskal.test}} (more than two groups) to each
#'            principal component.}
#'    }
#'  }
#' }
#'
#' @export
#'
#' @examples
#'
#' # Comparing the first output ("Pop.Sheep") of one the provided datasets.
#' cmp <-
#'  cmpoutput("SheepPop", 0.8, pphpc_ok$data[["Pop.Sheep"]], pphpc_ok$obs_lvls)
#'
#' # Compare bogus outputs created from 2 random sources, 5 observations per
#' # source, 20 variables each, yielding a 10 x 20 data matrix.
#' data <- matrix(c(rnorm(100), rnorm(100, mean = 1)), nrow = 10, byrow = TRUE)
#' olvls <- factor(c(rep("A", 5), rep("B", 5)))
#' cmp <- cmpoutput("Bogus", 0.7, data, olvls)
#'
cmpoutput <- function(name, ve_npcs, data, obs_lvls, lim_npcs = TRUE,
                      mnv_test = "Pillai") {

  # Check parameters
  if (any(ve_npcs <= 0))
    stop("'ve_npcs' parameter must only have positive values.")
  if (length(obs_lvls) != dim(data)[1])
    stop("Number of observations in 'data' and 'obs_lvls' does not match.")
  if (nlevels(obs_lvls) < 2)
    stop("At least two levels are required to perform model comparison.")

  # Perform PCA
  pca <- stats::prcomp(data)

  # Explained variances
  eig <- (pca$sdev) ^ 2
  varexp <- eig / sum(eig)
  cumvar <- cumsum(varexp)
  nve <- length(ve_npcs)

  # Minimum number of observations per group
  min_obs <- min(table(obs_lvls));

  # Pre-allocate vectors for Manova test
  npcs <- vector(mode = "integer", length = nve)
  mnvtest <- vector(mode = "list", length = nve)
  mnvpval <- vector(mode = "numeric", length = nve)

  # Perform a Manova test for each specified ve_npcs (variance to explain or
  # number of principal components)
  for (i in 1:nve) {

    # Percentage of variance to explain or number of PCs?
    if (ve_npcs[i] < 1) {
      # Percentage of variance to explain

      # Determine number of PCs required to explain the specified percentage of
      # variance
      npcs[i] <- which(cumvar > ve_npcs[i])[1]

    } else {
      # Number of PCs

      # Number of PCs given directly, use only integer part
      npcs[i] <- floor(ve_npcs[i])

      # Keep the variance explained by the specified number of PCs
      ve_npcs[i] <- cumvar[npcs[i]]

    }

    # Check number of dimensions (PCs), which should lower than the number of
    # observations
    if (min_obs < npcs[i]) {
      # Not enough observations for the specified number of PCs....

      # What action to take?
      if (lim_npcs) {

        # Limit number of PCs
        warning(paste0("Number of principal components for MANOVA test (",
                      npcs[i],") is higher than the size of the smallest ",
                      "group (", min_obs, "). Using only the first ", min_obs,
                      " principal components."))
        npcs[i] <- min_obs
        ve_npcs[i] <- cumvar[min_obs]

      } else {

        # Don't limit number of PCs and let MANOVA be performed in this less
        # than ideal situation
        warning(paste0("Number of principal components for MANOVA test (",
                       npcs[i],") is higher than the size of the smallest",
                       "group (", min_obs, ")."))

      }
    }

    # Manova
    if (npcs[i] > 1) {

      # Can only use Manova with more than one dimension
      mnvtest[[i]] <- stats::manova(pca$x[, 1:npcs[i]] ~ obs_lvls)
      mnvpval[i] <- summary(mnvtest[[i]], test = mnv_test)$stats[1, 6]

    } else {

      # Only one dimension, can't use Manova
      mnvtest[[i]] <- NULL
      mnvpval[i] <- NA

    }
  }

  # Total number of PCs returned by PCA operation
  tpcs <- length(eig);

  # Univariate tests
  parpvals <- vector(mode = "numeric", length = tpcs)
  parpvals_adjusted <- vector(mode = "numeric", length = tpcs)
  partests <- list()
  nonparpvals <- vector(mode = "numeric", length = tpcs)
  nonparpvals_adjusted <- vector(mode = "numeric", length = tpcs)
  nonpartests <- list()

  if (nlevels(obs_lvls) == 2) {
    # Use two-group tests

    # Cycle through each PC
    for (i in 1:tpcs) {

      # Parametric test (t-test) for each PC
      partests[[i]] <- stats::t.test(pca$x[, i] ~ obs_lvls, var.equal = T)
      parpvals[i] <- partests[[i]]$p.value

      # Non-parametric test (Mann-Whitney) for each PC
      nonpartests[[i]] <- stats::wilcox.test(pca$x[, i] ~ obs_lvls)
      nonparpvals[i] <- nonpartests[[i]]$p.value

    }

  } else {
    # Use multi-group tests (npcs > 2)

    # Cycle through each PC
    for (i in 1:tpcs) {

      # Parametric test (ANOVA) for each PC
      partests[[i]] <- stats::aov(pca$x[, i] ~ obs_lvls)
      parpvals[i] <- summary(partests[[i]])[[1]]$"Pr(>F)"[1]

      # Non-parametric test (Kruskal-Wallis) for each PC
      nonpartests[[i]] <- stats::kruskal.test(pca$x[, i] ~ obs_lvls)
      nonparpvals[i] <- nonpartests[[i]]$p.value

    }

  }

  # Determine adjusted univariate p-values using the weighted Bonferroni
  # procedure, explained variances used as weights
  parpvals_adjusted <- pmin(parpvals / varexp, 1)
  nonparpvals_adjusted <- pmin(nonparpvals / varexp, 1)

  # Return cmpoutput object
  cmpout <- list(scores = pca$x,
                 obs_lvls = obs_lvls,
                 varexp = varexp,
                 npcs = npcs,
                 ve = ve_npcs,
                 name = name,
                 p.values = list(manova = mnvpval,
                                 parametric = parpvals,
                                 nonparametric = nonparpvals,
                                 parametric_adjusted = parpvals_adjusted,
                                 nonparametric_adjusted = nonparpvals_adjusted),
                 tests = list(manova = mnvtest,
                              parametric = partests,
                              nonparametric = nonpartests))
  class(cmpout) <- "cmpoutput"
  cmpout

}

#' Print information about comparison of an output
#'
#' Print information about objects of class \code{cmpoutput}.
#'
#' @param x Object of class \code{cmpoutput}.
#' @param ... Currently ignored.
#'
#' @return The argument \code{x}, invisibly, as for all \code{\link{print}}
#' methods.
#'
#' @export
#'
#' @examples
#'
#' # Comparing the fifth output of the pphpc_diff dataset, which contains
#' # simulation output data from two implementations of the PPHPC model executed
#' # with a different parameter.
#'
#' cmpoutput("WolfPop", 0.7, pphpc_diff$data[[5]], pphpc_diff$obs_lvls)
#'
print.cmpoutput <- function(x, ...) {

  if (length(unique(x$obs_lvls)) == 2) {
    test_names <- c("t-test", "Mann-Whitney U test")
  } else {
    test_names <- c("ANOVA test", "Kruskal-Wallis test")
  }

  if (length(x$npcs) == 1) {
    chopen <- ""
    chclose <- ""
  } else {
    chopen <- "["
    chclose <- "]"
  }

  cat("Output name:", x$name, "\n")
  cat("Number of PCs which explain ", chopen,
      paste(sprintf("%2.1f", x$ve * 100), collapse = ", "), chclose,
      "% of variance: ", chopen, paste(x$npcs, collapse = ", "), chclose,
      "\n", sep = "")
  cat("P-Value for MANOVA along ", chopen, paste(x$npcs, collapse = ", "),
      chclose, " dimensions: ", chopen,
      paste(sprintf("%g", x$p.values$manova), collapse = ", "),
      chclose, "\n", sep = "")
  cat("P-Value for", test_names[1], "(1st PC):",
      x$p.values$parametric[1], "\n")
  cat("P-Value for", test_names[2], "(1st PC):",
      x$p.values$nonparametric[1], "\n")
  cat("Adjusted p-Value for", test_names[1], "(1st PC):",
      x$p.values$parametric_adjusted[1], "\n")
  cat("Adjusted p-Value for", test_names[2], "(1st PC):",
      x$p.values$nonparametric_adjusted[1], "\n")

  invisible(x)

}

#' Summary method for comparison of an output
#'
#' Summary method for objects of class \code{cmpoutput}.
#'
#' @param object Object of class \code{cmpoutput}.
#' @param ... Currently ignored.
#'
#' @return A list with the following components:
#' \describe{
#'  \item{output.name}{Output name.}
#'  \item{num.pcs}{Number of principal components which explain \code{var.exp}
#'        percentage of variance.}
#'  \item{var.exp}{Minimum percentage of variance which must be explained by the
#'        number of principal components used for the MANOVA test.}
#'  \item{manova.pvals}{\emph{P}-value of the MANOVA test.}
#'  \item{parametric.test}{Name of the used parametric test.}
#'  \item{parametric.pvals}{Vector of $p$-values returned by applying the
#'        parametric test to each principal component.}
#'  \item{parametric.pvals.adjusted}{Vector of $p$-values returned by applying
#'        the parametric test to each principal component, adjusted with the
#'        weighted Bonferroni procedure, percentage of explained variance used
#'        as weight.}
#'  \item{nonparametric.test}{Name of the used non-parametric test.}
#'  \item{nonparametric.pvals}{Vector of $p$-values returned by applying the
#'        non-parametric test to each principal component.}
#'  \item{nonparametric.pvals.adjusted}{Vector of $p$-values returned by
#'        applying the non-parametric test to each principal component, adjusted
#'        with the weighted Bonferroni procedure, percentage of explained
#'        variance used as weight.}
#' }
#'
#' @export
#'
#' @examples
#'
#' # Comparing the concatenated output of the pphpc_noshuff dataset, which
#' # contains simulation output data from two implementations of the PPHPC model
#' # executed with a minor implementation difference.
#'
#' summary(
#'   cmpoutput("All", 0.6, pphpc_noshuff$data[["All"]], pphpc_noshuff$obs_lvls)
#' )
#'
summary.cmpoutput <- function(object, ...) {

  if (length(unique(object$obs_lvls)) == 2) {
    test_names <- c("t-test", "Mann-Whitney")
  } else {
    test_names <- c("ANOVA", "Kruskal-Wallis")
  }

  list(output.name = object$name,
       num.pcs = object$npcs,
       var.exp = object$ve,
       manova.pvals = object$p.values$manova,
       parametric.test = test_names[1],
       parametric.pvals = object$p.values$parametric,
       parametric.pvals.adjusted = object$p.values$parametric_adjusted,
       nonparametric.test = test_names[2],
       nonparametric.pvals = object$p.values$nonparametric,
       nonparametric.pvals.adjusted = object$p.values$nonparametric_adjusted)
}

#' Plot comparison of an output
#'
#' Plot objects of class \code{cmpoutput}.
#'
#' This method produces four sub-plots, namely:
#' \itemize{
#'   \item Scatter plot containing the projection of output observations on the
#'         first two dimensions of the principal components space.
#'   \item Bar plot of the percentage of variance explain per principal
#'         component.
#'   \item Bar plot of \emph{p}-values for the parametric test for each
#'         principal component.
#'   \item Bar plot of \emph{p}-values for the non-parametric test for each
#'         principal component.
#' }
#'
#' @param x Object of class \code{cmpoutput}.
#' @param ... Extra options passed to \code{\link[graphics]{plot.default}}. The
#' \code{col} option determines the colors to use on observations of different
#' groups (scatter plot only).
#'
#' @return None.
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#'
#' # Comparing the concatenated output of the pphpc_ok dataset, which
#' # contains simulation output data from two similar implementations of the
#' # PPHPC model.
#'
#' plot(cmpoutput("All", 0.95, pphpc_ok$data[["All"]], pphpc_ok$obs_lvls))
#'
plot.cmpoutput <- function(x, ...) {

  # Was a color specified?
  params <- list(...)
  if (exists("col", where = params)) {
    col <- params$col
    params$col <- NULL # We don't want any color in the barplots
  } else {
    col <- plotcols()
  }

  # Set mfrow graphical parameter to setup subplots
  graphics::par(mfrow = c(3,2))

  # Total number of PCs
  tpcs <- length(x$varexp)

  # Score plot (first two PCs)
  params_sp <- params
  params_sp$x <- x$scores[, 1]
  params_sp$y <- x$scores[, 2]
  params_sp$col <- col[as.numeric(x$obs_lvls)]
  params_sp$xlab <- paste("PC1 (", round(x$varexp[1] * 100, 2), "%)", sep = "")
  params_sp$ylab <- paste("PC2 (", round(x$varexp[2] * 100, 2), "%)", sep = "")
  params_sp$main <- "Score plot"
  do.call("plot.default", params_sp)

  # Explained variance bar plot
  params_ve <- params
  params_ve$height <- x$varexp
  params_ve$names.arg <- as.character(1:tpcs)
  params_ve$main <- "Explained variance by PC"
  params_ve$xlab <- "PC"
  params_ve$ylab <- "Var. exp. (%)"
  do.call("barplot", params_ve)

  # Parametric p-values bar plot
  params_ppv <- params
  params_ppv$height <- x$p.values$parametric
  params_ppv$names.arg <- as.character(1:tpcs)
  params_ppv$main <- "Parametric p-values by PC"
  params_ppv$xlab <- "PC"
  params_ppv$ylab <- "Prob."
  do.call("barplot", params_ppv)

  # Non-parametric p-values bar plot
  params_nppv <- params
  params_nppv$height <- x$p.values$nonparametric
  params_nppv$names.arg <- as.character(1:tpcs)
  params_nppv$main <- "Non-parametric p-values by PC"
  params_nppv$xlab <- "PC"
  params_nppv$ylab <- "Prob."
  do.call("barplot", params_nppv)

  # Parametric adjusted p-values bar plot
  params_papv <- params
  params_papv$height <- x$p.values$parametric_adjusted
  params_papv$names.arg <- as.character(1:tpcs)
  params_papv$main <- "Parametric p-values by PC (adjusted)"
  params_papv$xlab <- "PC"
  params_papv$ylab <- "Prob."
  do.call("barplot", params_papv)

  # Non-parametric adjusted p-values bar plot
  params_npapv <- params
  params_npapv$height <- x$p.values$nonparametric_adjusted
  params_npapv$names.arg <- as.character(1:tpcs)
  params_npapv$main <- "Non-parametric p-values by PC (adjusted)"
  params_npapv$xlab <- "PC"
  params_npapv$ylab <- "Prob."
  do.call("barplot", params_npapv)

  invisible(NULL)
}

#' Get assumptions for parametric tests performed on output comparisons
#'
#' Get assumptions for parametric tests performed on output comparisons (i.e.
#' from objects of class \code{\link{cmpoutput}}).
#'
#' @param obj Object of class \code{cmpoutput}.
#'
#' @return Object of class \code{assumptions_cmpoutput} containing the
#' assumptions for parametric tests performed on an output comparison.
#' Basically a list containing the assumptions for the MANOVA (list of objects
#' of class \code{\link{assumptions_manova}}, one per explained variance) and
#' univariate parametric tests for each principal component (object of class
#' \code{\link{assumptions_paruv}}).
#'
#' @export
#'
#' @examples
#'
#' # Create a cmpoutput object from the provided datasets
#' cmp <- cmpoutput("All", 0.9, pphpc_ok$data[["All"]], pphpc_ok$obs_lvls)
#'
#' # Get the assumptions for the parametric tests performed in cmp
#' acmp <- assumptions(cmp)
#'
assumptions.cmpoutput <- function(obj) {

  # Create assumptions object
  assumptions <- list()
  class(assumptions) <- "assumptions_cmpoutput"

  # Get stuff from cmpoutput object
  npcs <- obj$npcs
  tpcs <- length(obj$varexp)
  obs_lvls <- obj$obs_lvls
  scores <- obj$scores
  nve <- length(npcs)

  # Allocate vector for Manova assumptions
  assumptions$manova <- vector(mode = "list", length = nve)

  # Determine Manova assumptions
  for (i in 1:nve) {
    if (npcs[i] > 1) {

      # Can only use manova if more than one variable
      assumptions$manova[[i]] <-
        assumptions_manova(scores[, 1:npcs[i]], obs_lvls)

      # Keep number of PCs
      attr(assumptions$manova[[i]], "npcs") <- npcs[i]

    } else {

      # Only one variable, can't use manova
      assumptions$manova[[i]] <- NULL

    }
  }

  # Parametric test (t-test) for each PC
  assumptions$ttest <- assumptions_paruv(scores[, 1:tpcs], obs_lvls)

  # Return assumptions
  assumptions

}

#' Plot \emph{p}-values for testing the assumptions of the parametric tests used
#' in output comparison
#'
#' Plot method for objects of class \code{assumptions_cmpoutput}
#' containing \emph{p}-values produced by testing the assumptions of the
#' parametric tests used for comparing an output.
#'
#' Several bar plots are presented, showing the \emph{p}-values yielded by the
#' Shapiro-Wilk (\code{\link[stats]{shapiro.test}}) and Royston tests
#' (\code{\link[MVN]{roystonTest}}) for univariate and multivariate normality,
#' respectively, and for the Bartlett (\code{\link[stats]{bartlett.test}}) and
#' Box's M (\code{\link[biotools]{boxM}}) for testing homogeneity of variances
#' and of covariance matrices, respectively. The following bar plots are shown:
#'
#' \itemize{
#'  \item One bar plot for the \emph{p}-values of the Bartlett test, one bar
#'        (\emph{p}-value) per individual principal component.
#'  \item \emph{s} bar plots for \emph{p}-values of the Shapiro-Wilk test, where
#'        \emph{s} is the number of groups being compared. Individual bars in
#'        each plot are associated with a principal component.
#'  \item \emph{t} bar plot for the \emph{p}-values of the Royston test with
#'        \emph{s} bars each, where \emph{t} is the number of unique MANOVA
#'        tests performed (one per requested explained variances) and \emph{s}
#'        is the number of groups being compared. These plots will not show if
#'        there is only one principal component being considered.
#'  \item One plot for the \emph{p}-values of the Box's M test, one bar
#'        (\emph{p}-value) per unique MANOVA tests performed  (one per requested
#'        explained variances).
#' }
#'
#' @param x Objects of class \code{assumptions_cmpoutput}.
#' @param ... Extra options passed to \code{\link[graphics]{plot.default}}.
#'
#' @return None.
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#'
#' # Create a cmpoutput object from the provided datasets
#' cmp <- cmpoutput("All", 0.9, pphpc_ok$data[["All"]], pphpc_ok$obs_lvls)
#'
#' # Display a bar plot with the p-values of the assumptions for the parametric
#' # tests performed in cmp
#' plot(assumptions(cmp))
#'
plot.assumptions_cmpoutput <- function(x, ...) {

  # Multivariate assumptions

  # How many PCs did each MANOVA test?
  npcs <- sapply(x$manova,
                 function(y) {
                   if (!is.null(y) && methods::is(y$mvntest[[1]], "royston")) {
                     dim(y$mvntest[[1]]@dataframe)[2]
                   } else {
                     1
                   }
                 })

  # We don't need to re-plot for the same number of PCs
  npcs[duplicated(npcs)] <- 1

  # How many are greater than 1?
  nmnvmvplt <- sum(npcs > 1)

  # Filter number of MANOVA multivariate assumptions to plot
  mnvmv <- x$manova[npcs > 1]

  # More than one? Then plot also Box's M p-values for different number of PCs
  if (length(mnvmv) > 1) {

    nmnvboxplt <- 1
    mnvboxp <- sapply(mnvmv, function(x) x$vartest$p.value)

  } else {

    # No Box's M p-values plot
    nmnvboxplt <- 0

  }

  # How many plots?
  nplots <- length(x$ttest) + 1 + nmnvmvplt + nmnvboxplt

  # Determine layout matrix side dimension
  side_dim <- ceiling(sqrt(nplots))

  # Setup subplot layout
  graphics::par(mfrow = c(side_dim, side_dim))

  # Plot univariate assumptions
  plot(x$ttest, ...)

  # Plot multivariate assumptions
  for (amnv in mnvmv) {
    plot(amnv, ...)
  }

  # Plot Box's M in case there are several p-values to plot
  if (nmnvboxplt) {

    # Was a color specified?
    params <- list(...)
    if (exists("col", where = params)) {
      col <- params$col
      params$col <- NULL
    } else {
      col <- c("darkgreen", "yellow", "red")
    }

    # Plot the p-values in a bar plot
    params$height <- mnvboxp
    params$main <- "Box's M test p-values"
    params$sub <- "Homogeneity of Covariance Matrices"
    params$xlab <- "Number of PCs"
    params$names.arg <- npcs[npcs > 1]
    params$ylab <- "Probability"
    params$col <- pvalcol(mnvboxp, col)
    do.call(get("barplot", asNamespace("graphics")), params)

  }

  invisible(NULL)

}

#' Print method for the assumptions of parametric tests used in a comparison
#' of an output
#'
#' Print method for objects of class \code{assumptions_cmpoutput}, which
#' contain the assumptions for the parametric tests used in a comparison of an
#' output.
#'
#' @param x Object of class \code{assumptions_cmpoutput}.
#' @param ... Currently ignored.
#'
#' @return None.
#'
#' @export
#'
#' @examples
#'
#' # Create a cmpoutput object from the provided datasets
#' cmp <- cmpoutput("All", c(0.7, 0.8, 0.9),
#'                  pphpc_diff$data[["All"]], pphpc_diff$obs_lvls)
#'
print.assumptions_cmpoutput <- function(x, ...) {

  # Obtain object summary
  sa <- summary(x)

  # Nice print of object summary
  cat("=== MANOVA assumptions ===\n")
  if (is.null(sa$manova)) {
    cat("No MANOVA tests were performed.\n")
  } else {
    print(sa$manova)
  }
  cat("\n=== T-test assumptions ===\n")
  print(sa$ttest[, 1, drop = F])

  invisible(NULL)

}

#' Summary method for the assumptions of parametric tests used in a comparison
#' of an output
#'
#' Summary method for objects of class \code{assumptions_cmpoutput}, which
#' contain the assumptions for the parametric tests used in a comparison of an
#' output.
#'
#' @param object Object of class \code{assumptions_cmpoutput}.
#' @param ... Currently ignored.
#'
#' @return A list with the following items:
#' \describe{
#'  \item{manova}{A matrix of \emph{p}-values for the MANOVA assumptions. All
#'        rows, expect the last one, correspond to the Royston test for
#'        multivariate normality for each group; the last row corresponds to
#'        Box's M test for homogeneity of covariance matrices. Columns
#'        correspond to number of principal components required to explain the
#'        percentage of user-specified variance.}
#'  \item{ttest}{A matrix of \emph{p}-values for the \emph{t}-test assumptions.
#'        All rows, expect the last one, correspond to the Shapiro-Wilk
#'        normality test for each group; the last row corresponds to Bartlett's
#'        for equality of variances. Columns correspond to the principal
#'        components on which the \emph{t}-test was applied.}
#' }
#'
#' @export
#'
#' @examples
#'
#' # Create a cmpoutput object from the provided datasets
#' cmp <- cmpoutput("All", c(0.5, 0.6, 0.7),
#'                  pphpc_ok$data[["All"]], pphpc_ok$obs_lvls)
#'
#' # Obtain the summary of the assumptions of the cmpoutput object
#' summary(assumptions(cmp))
#'
summary.assumptions_cmpoutput <- function(object, ...) {

  # Multivariate assumptions

  # How many PCs did each MANOVA test?
  npcs <- sapply(object$manova,
                 function(y) {
                   if (!is.null(y) && methods::is(y$mvntest[[1]], "royston")) {
                     dim(y$mvntest[[1]]@dataframe)[2]
                   } else {
                     1
                   }
                 })

  # We don't need to repeat the p-values for the same number of PCs
  npcs[duplicated(npcs)] <- 1

  # Filter number of MANOVA multivariate assumptions to consider
  mnvmv <- object$manova[npcs > 1]

  # Get p-values for multivariate assumptions
  if (length(mnvmv) > 0) {
    mvpvals <- sapply(mnvmv,
                      function(mnv) {
                        npv <- sapply(mnv$mvntest,
                                      function(roy) roy@p.value)
                        pv <- c(npv, mnv$vartest$p.value)
                        names(pv) <-
                          c(paste0("Royston (", names(mnv$mvntest), ")"),
                            "Box's M")
                        pv
                      })
    colnames(mvpvals) <- paste0("NPCs=", npcs[npcs > 1])
  } else {
    mvpvals <- NULL
  }

  # Univariate assumptions

  # How many PCs?
  tnpcs <- length(object$ttest$uvntest[[1]])

  # How many p-values per PC? Number of groups (normality) + Bartlett (variance)
  npvals <- length(object$ttest$uvntest) + 1

  # Allocate p-value matrix and set row and col names
  uvpvals <- matrix(nrow = npvals, ncol = tnpcs)
  rownames(uvpvals) <-
    c(paste0("Shapiro-Wilk (", names(object$ttest$uvntest), ")"), "Bartlett")
  colnames(uvpvals) <- paste0("PC", 1:tnpcs)

  # Get the p-values for for t-test assumptions and put them into matrix
  for (i in 1:tnpcs) {

    # Normality assumption p-values (Shapiro-Wilk)
    uvpvals[1:(npvals - 1), i] <-
      sapply(object$ttest$uvntest, function(grp) grp[[i]]$p.value)

    # Variance assumption p-value (Bartlett)
    uvpvals[npvals, i] <- object$ttest$vartest[[i]]$p.value

  }

  # Return summary
  list(manova = mvpvals, ttest = uvpvals)
}
