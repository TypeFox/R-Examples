# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

#' Multiple independent comparisons of observations
#'
#' Performs multiple independent comparisons of output observations.
#'
#' @param outputs A vector with the labels of each output, or an integer with
#' the number of outputs (in which case output labels will be assigned
#' automatically).
#' @param ve_npcs Percentage (\code{0 < ve_npcs < 1}) of variance explained by
#' the \emph{q} principal components (i.e. number of dimensions) used in MANOVA,
#' or the number of principal components (\code{ve_npcs > 1}, must be integer).
#' Can be a vector, in which case the MANOVA test will be applied multiple
#' times, one per specified variance to explain / number of principal
#' components.
#' @param comps A list of lists, where each list contains information regarding
#' an individual comparison. Each list can have one of two configurations:
#' \enumerate{
#'   \item Lists with the first configuration are used to load data from files,
#'         and require the following fields:
#'     \describe{
#'       \item{name}{A string specifying the comparison name.}
#'       \item{folders}{Vector of folder names where to read files from. These
#'             are recycled if \code{length(folders) < length(files)}.}
#'       \item{files}{Vector of filenames (with wildcards) to load in each
#'             folder.}
#'       \item{lvls}{Vector of level or group names, must be the same length as
#'             \code{files}, i.e. each file set will be associated with a
#'             different group. If not given, default group names will be set.}
#'     }
#'   \item Lists with the second configuration are used to load data from
#'         environment variables, and require the following fields:
#'     \describe{
#'       \item{name}{A string specifying the comparison name.}
#'       \item{grpout}{Either an object of class \code{\link{grpoutputs}} or a
#'             list with the following two fields:
#'         \describe{
#'           \item{data}{List of all outputs, where tags correspond to output
#'                 names and values correspond to the output data. Output data
#'                 is a \emph{n} x \emph{m} matrix, where \emph{n} is the total
#'                 number of output observations and \emph{m} is the number of
#'                 variables (i.e. output length). }
#'           \item{obs_lvls}{Levels or groups associated with each observation.}
#'         }
#'       }
#'     }
#' }
#' @param concat Create an additional, concatenated output? Ignored for sublists
#' passed in the \code{comps} which follow the second configuration.
#' @param centscal Method for centering and scaling outputs if \code{concat} is
#' TRUE. It can be one of "center", "auto", "range" (default), "iqrange",
#' "vast", "pareto" or "level". Centering and scaling is performed by the
#' \code{\link{centerscale}} function.
#' @param lim_npcs Limit number of principal components used for MANOVA to
#' minimum number of observations per group?
#' @param mnv_test The name of the test statistic to be used in MANOVA, as
#' described in \code{\link[stats]{summary.manova}}.
#' @param ... Options passed to \code{\link[utils]{read.table}}, which is used
#' to read the files specified in lists using the first configuration in the
#' \code{comp} parameter.
#'
#' @return An object of class \code{micomp}, which is a two-dimensional
#' list of \code{\link{cmpoutput}} objects. Rows are associated with individual
#' outputs, while columns are associated with separate comparisons.
#'
#' @export
#'
#' @examples
#'
#' # Create a micomp object from existing files and folders
#'
#' dir_nl_ok <-
#'   system.file("extdata", "nl_ok", package = "micompr")
#' dir_jex_ok <-
#'   system.file("extdata", "j_ex_ok", package = "micompr")
#' dir_jex_noshuff <-
#'   system.file("extdata", "j_ex_noshuff", package = "micompr")
#' dir_jex_diff <-
#'   system.file("extdata", "j_ex_diff", package = "micompr")
#' files <- glob2rx("stats400v1*.tsv")
#'
#' mic <- micomp(7, 0.8,
#'               list(list(name = "NLOKvsJEXOK",
#'                         folders = c(dir_nl_ok, dir_jex_ok),
#'                         files = c(files, files),
#'                         lvls = c("NLOK", "JEXOK")),
#'                    list(name = "NLOKvsJEXNOSHUFF",
#'                         folders = c(dir_nl_ok, dir_jex_noshuff),
#'                         files = c(files, files),
#'                         lvls = c("NLOK", "JEXNOSHUFF")),
#'                    list(name = "NLOKvsJEXDIFF",
#'                         folders = c(dir_nl_ok, dir_jex_diff),
#'                         files = c(files, files),
#'                         lvls = c("NLOK", "JEXDIFF"))),
#'               concat = TRUE)
#'
#' \donttest{
#' # Create a micomp object from package datasets (i.e. grpoutputs objects)
#' # directly
#'
#' mic <- micomp(c("o1", "o2", "o3", "o4"), 0.9,
#'               list(list(name = "NLOKvsJEXOK", grpout = pphpc_ok),
#'                    list(name = "NLOKvsJEXNOSHUFF", grpout = pphpc_noshuff),
#'                    list(name = "NLOKvsJEXDIFF", grpout = pphpc_diff)))
#'
#' # Create a micomp object using manually inserted data
#'
#' mic <- micomp(6, 0.5, list(
#'   list(name = "NLOKvsJEXOK",
#'        grpout = list(data = pphpc_ok$data,
#'                      obs_lvls = pphpc_ok$obs_lvls)),
#'   list(name = "NLOKvsJEXNOSHUFF",
#'        grpout = list(data = pphpc_noshuff$data,
#'                      obs_lvls = pphpc_noshuff$obs_lvls)),
#'   list(name = "NLOKvsJEXDIFF",
#'        grpout = list(data = pphpc_diff$data,
#'                      obs_lvls = pphpc_diff$obs_lvls))))
#' }
micomp <- function(
  outputs, ve_npcs, comps, concat = F, centscal = "range", lim_npcs = TRUE,
  mnv_test = "Pillai", ...) {

  # Determine number of comparisons
  ncomp <- length(comps)
  cmp_names <- vector(mode = "character", length = ncomp)

  # Did the user specify output names? Or should we provide some default names?
  if (length(outputs) == 1) {
    nout <- outputs
    outputs <- paste("out", 1:nout, sep = "")
  } else {
    nout <- length(outputs)
  }

  # List of grouped outputs for each comparison
  grpd_outputs <- list()

  # Results of each comparison for each output
  comp_res <- vector("list", length = nout * ncomp)
  dim(comp_res) <- c(nout, ncomp)

  # Group outputs for each comparison
  for (i in 1:ncomp) {

    # What kind of configuration does the current comparison have? File or
    # matrix?
    if (exists("name", where = comps[[i]]) &&
        exists("folders", where = comps[[i]]) &&
        exists("files", where = comps[[i]])) {

      # First configuration: load data from files

      grpd_outputs[[i]] <-
        grpoutputs(outputs = outputs,
                   folders = unlist(comps[[i]]$folders),
                   files = unlist(comps[[i]]$files),
                   lvls = comps[[i]]$lvls,
                   concat = concat,
                   centscal = centscal,
                   ...)

    } else if (exists("name", where = comps[[i]]) &&
               exists("grpout", where = comps[[i]])) {

      # Second configuration: load data from environment variables

      # Is the object in the grpout field of class grpoutputs?
      if (methods::is(comps[[i]], "grpoutputs")) {

        # Yes, so use it directly
        grpd_outputs[[i]] <- comps[[i]]

      } else {

        # No, so build it using the provided data
        olvls <- comps[[i]]$grpout$obs_lvls
        grpd_outputs[[i]] <- structure(
          list(data = comps[[i]]$grpout$data,
               groups = sapply(unique(olvls),
                               function(onelvl, lvls) sum(lvls == onelvl),
                               olvls),
               obs_lvls = olvls,
               lvls = levels(unique(olvls)),
               concat = F),
          class = "grpoutputs")

      }

    } else {

      # Unknown configuration, throw error
      stop("Invalid 'comps' parameter")

    }

    # Keep comparison names
    cmp_names[i] <-
      if (!is.null(comps[[i]]$name)) {
        comps[[i]]$name
      } else {
        paste("Comp", i)
      }


  }

  # Cycle through each output
  for (i in 1:nout) {

    # Cycle through comparisons
    for (j in 1:ncomp) {

      comp_res[[i, j]] <-
        cmpoutput(outputs[i], ve_npcs,
                  grpd_outputs[[j]]$data[[i]],
                  grpd_outputs[[j]]$obs_lvls,
                  lim_npcs = lim_npcs,
                  mnv_test = mnv_test)
    }

  }

  colnames(comp_res) <- cmp_names
  rownames(comp_res) <- outputs
  class(comp_res) <- "micomp"
  attr(comp_res, "ve_npcs") <- ve_npcs
  comp_res

}

#' Print information about multiple comparisons of outputs
#'
#' Print information about objects of class \code{\link{micomp}}.
#'
#' @param x Object of class \code{\link{micomp}}.
#' @param ... Currently ignored.
#'
#' @return The argument \code{x}, invisibly, as for all \code{\link{print}}.
#' methods.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # A micomp object from package datasets (i.e. grpoutputs objects) directly
#'
#' micomp(c("outA", "outB", "outC", "outD"), 0.9,
#'               list(list(name = "Comp1", grpout = pphpc_ok),
#'                    list(name = "Comp2", grpout = pphpc_noshuff),
#'                    list(name = "Comp3", grpout = pphpc_diff)))
#' }
print.micomp <- function(x, ...) {

  # Use summary to get the info to be printed
  smic <- summary(x)

  # Cycle through comparisons
  for (cmpname in names(smic)) {

    cat("====", cmpname, "====\n")
    print(smic[[cmpname]], digits = 5, print.gap = 2)
    cat("\n")

  }

  invisible(x)

}

#' Summary method for multiple comparisons of outputs
#'
#' Summary method for objects of class \code{\link{micomp}}.
#'
#' @param object Object of class \code{\link{micomp}}.
#' @param ... Currently ignored.
#'
#' @return A list in which each component is associated with a distinct
#' comparison. Each component contains a matrix, in which columns represent
#' individual outputs and rows have information about the outputs. More
#' specifically, each matrix has the following rows:
#' \describe{
#'  \item{#PCs (ve=\%)}{Number of principal components required to explain the
#'        specified percentage of variance. There is one row of this kind for
#'        each percentage of variance specified when creating the
#'        \code{\link{micomp}} object.}
#'  \item{MANOVA (ve=\%)}{\emph{P}-value for the MANOVA test applied to the #PCs
#'        required to explain the specified percentage of variance. There is one
#'        row of this kind for each percentage of variance specified when
#'        creating the \code{\link{micomp}} object.}
#'  \item{par.test}{\emph{P}-value for the parametric test (first principal
#'        component).}
#'  \item{nonpar.test}{\emph{P}-value for the non-parametric test (first
#'        principal component).}
#'  \item{par.test.adjust}{\emph{P}-value for the parametric test (first
#'        principal component), adjusted with the weighted Bonferroni procedure,
#'        percentage of explained variance used as weight.}
#'  \item{nonpar.test.adjust}{\emph{P}-value for the non-parametric test (first
#'        principal component), adjusted with the weighted Bonferroni procedure,
#'        percentage of explained variance used as weight.}
#' }
#'
#' @export
#'
#' @examples
#'
#' # A micomp object from package datasets (i.e. grpoutputs objects) directly
#' \donttest{
#' summary(micomp(5, 0.85,
#'                list(list(name = "CompEq", grpout = pphpc_ok),
#'                     list(name = "CompNoShuf", grpout = pphpc_noshuff),
#'                     list(name = "CompDiff", grpout = pphpc_diff))))
#' }
summary.micomp <- function(object, ...) {

  dims <- dim(object)
  nout <- dims[1]
  ncomp <- dims[2]
  cmpnames <- colnames(object)

  smic <- vector(mode = "list", length = ncomp)
  names(smic) <- cmpnames

  # How many variances to explain?
  ve_npcs <- attr(object, "ve_npcs")
  nve_npcs <- length(ve_npcs)

  # Cycle through comparisons
  for (i in 1:ncomp) {

    # Allocate summary matrix for current comparison
    smat <- matrix(nrow = 2 * nve_npcs + 4, ncol = nout)
    colnames(smat) <- rownames(object)
    rnames <- NULL

    # Obtain information from the MANOVA tests, one per requested variance to
    # explain / number of PCs
    for (j in 1:nve_npcs) {
      idx <- (j - 1) * 2 + 1

      # Get number of PCs and MANOVA p-value
      smat[idx, ] <- sapply(object[, i],
                            function(mc) mc$npcs[j])
      smat[idx + 1, ] <- sapply(object[, i],
                                function(mc) mc$p.values$manova[j])

      # Show variance to explain or number of PCs in info label?
      infolabel <-
        if (ve_npcs[j] < 1) {
          paste0("ve=", ve_npcs[j])
        } else {
          paste0("npcs=", ve_npcs[j])
        }

      # Add two more row names, for number of PCs and MANOVA, respectively
      rnames <- c(rnames,
                  paste0("#PCs (", infolabel, ")"),
                  paste0("MANOVA (", infolabel, ")"))
    }

    # Obtain information for the univariate tests for the 1st PC
    smat[idx + 2, ] <-
      sapply(object[, i],
             function(mc) mc$p.values$parametric[1])
    smat[idx + 3, ] <-
      sapply(object[, i],
             function(mc) mc$p.values$nonparametric[1])
    smat[idx + 4, ] <-
      sapply(object[, i],
             function(mc) mc$p.values$parametric_adjusted[1])
    smat[idx + 5, ] <-
      sapply(object[, i],
             function(mc) mc$p.values$nonparametric_adjusted[1])

    # Set summary matrix row names
    rownames(smat) <- c(rnames, "par.test", "nonpar.test",
                        "par.test.adjust", "nonpar.test.adjust")

    smic[[cmpnames[i]]] <- smat

  }

  smic

}

#' Plot projection of output observations on the first two dimensions of the
#' principal components space
#'
#' For each comparison and output combination, draw a scatter plot containing
#' the projection of output observations on the first two dimensions of the
#' principal components space.
#'
#' @param x An object of class \code{\link{micomp}}.
#' @param ... Extra options passed to \code{\link[graphics]{plot.default}}. The
#' \code{col} option determines the colors to use on observations of different
#' groups.
#'
#' @return None.
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#' \donttest{
#' plot(micomp(c("SheepPop", "WolfPop", "GrassQty"), 0.95,
#'             list(list(name = "I", grpout = pphpc_ok),
#'                  list(name = "II", grpout = pphpc_noshuff),
#'                  list(name = "III", grpout = pphpc_diff))))
#' }
plot.micomp <- function(x, ...) {

  # Was a color specified?
  params <- list(...)
  if (!exists("col", where = params)) {
    params$col <- plotcols()
  }
  col <- params$col

  # Useful variables
  dims <- dim(x)
  nout <- dims[1]
  ncomp <- dims[2]
  nplots <- nout * ncomp

  # Layout of subplots
  m <- matrix(1:(nplots + ncomp), nrow = ncomp, ncol = nout + 1, byrow = T)
  graphics::layout(mat = m)

  for (i in 1:ncomp) {

    # Get levels associated with each observation from the first output of the
    # current comparison
    obs_lvls <- x[[1, i]]$obs_lvls

    # Set title and legend for current comparison
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
    graphics::text(1, 1, pos = 1, labels = paste("Comp. ", i))
    graphics::legend("center", legend = unique(obs_lvls), fill = col, horiz = F)

    for (j in 1:nout) {

      # Get data
      scores <- x[[j, i]]$scores
      varexp <- x[[j, i]]$varexp

      # Score plot (first two PCs)
      params$x <- scores[, 1]
      params$y <- scores[, 2]
      params$col <- col[as.numeric(obs_lvls)]
      params$xlab <- paste("PC1 (", round(varexp[1] * 100, 2), "%)", sep = "")
      params$ylab <- paste("PC2 (", round(varexp[2] * 100, 2), "%)", sep = "")
      params$main <- x[[j, i]]$name
      do.call(get("plot.default", asNamespace("graphics")), params)

    }
  }

  invisible(NULL)

}

#' Get assumptions for parametric tests performed on each comparisons
#'
#' Get assumptions for parametric tests performed on multiple comparisons (i.e.
#' from objects of class \code{\link{micomp}}).
#'
#' @param obj Object of class \code{\link{micomp}}.
#'
#' @return Object of class \code{assumptions_micomp} containing the
#' assumptions for parametric tests performed for the multiple comparisons held
#' by the \code{mcmp} object. This object is a multi-dimensional list of
#' \code{assumptions_cmpoutput} objects. Rows are associated with individual
#' outputs, while columns are associated with separate comparisons.
#'
#' @export
#'
#' @examples
#'
#' # Create a micomp object, use provided dataset
#' mic <- micomp(6, 0.8,
#'               list(list(name = "NLOKvsJEXOK", grpout = pphpc_ok),
#'                    list(name = "NLOKvsJEXNOSHUFF", grpout = pphpc_noshuff),
#'                    list(name = "NLOKvsJEXDIFF", grpout = pphpc_diff)))
#'
#' # Create an object containing the statistic tests evaluating the assumptions
#' # of the comparisons performed in the mic object
#' a <- assumptions(mic)
#'
assumptions.micomp <- function(obj) {
  micas <- lapply(obj, function(x) assumptions(x))
  dim(micas) <- dim(obj)
  colnames(micas) <- colnames(obj)
  rownames(micas) <- rownames(obj)
  class(micas) <- "assumptions_micomp"
  attr(micas, "ve_npcs") <- attr(obj, "ve_npcs")
  micas
}

#' Print information about the assumptions concerning the parametric tests
#' performed on multiple comparisons of outputs
#'
#' Print information about objects of class \code{assumptions_micomp}, which
#' represent the assumptions concerning the parametric tests performed on
#' multiple comparisons of outputs.
#'
#' @param x Object of class \code{assumptions_micomp}.
#' @param ... Currently ignored.
#'
#' @return The argument \code{x}, invisibly, as for all \code{\link{print}}
#' methods.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Create a micomp object, use provided dataset
#' mic <- micomp(c("SheepPop", "WolfPop", "GrassQty"), 0.7,
#'               list(list(name = "NLOKvsJEXOK", grpout = pphpc_ok),
#'                    list(name = "NLOKvsJEXNOSHUFF", grpout = pphpc_noshuff),
#'                    list(name = "NLOKvsJEXDIFF", grpout = pphpc_diff)))
#'
#' # Print the results (p-values) of the statistic tests evaluating the
#' # assumptions of the comparisons performed in the mic object
#' assumptions(mic)
#' }
print.assumptions_micomp <- function(x, ...) {

  sm <- summary(x)

  # Cycle through comparisons
  for (cmp in names(sm)) {

    cat("==== ", cmp, "====\n")
    print(sm[[cmp]], digits = 5, print.gap = 2)
    cat("\n")

  }

  invisible(x)

}

#' Plot \emph{p}-values for testing the assumptions of the parametric tests used
#' in multiple output comparison
#'
#' Plot method for objects of class \code{assumptions_cmpoutput}
#' containing \emph{p}-values produced by testing the assumptions of the
#' parametric tests used for multiple output comparisons.
#'
#' Several bar plots are presented, one for each comparison and output
#' combination, showing the several statistical tests employed to verify
#' the assumptions of the parametric tests.
#'
#' @param x Object of class \code{assumptions_micomp}.
#' @param ... Extra options passed to \code{\link[graphics]{barplot}}.
#'
#' @return None.
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#' \donttest{
#' # Create a micomp object, use provided dataset
#' mic <- micomp(6, 0.65,
#'               list(list(name = "NLOKvsJEXOK", grpout = pphpc_ok),
#'                    list(name = "NLOKvsJEXNOSHUFF", grpout = pphpc_noshuff),
#'                    list(name = "NLOKvsJEXDIFF", grpout = pphpc_diff)))
#'
#' # Plot the p-values of the statistic tests evaluating the assumptions of the
#' # comparisons performed in the mic object
#' plot(assumptions(mic))
#' }
plot.assumptions_micomp <- function(x, ...) {

  # Get the assumptions summary
  sm <- summary(x)

  # Set layout for plots
  graphics::par(mfcol = dim(x))

  # Cycle through comparisons
  for (cmp in names(sm)) {

    # Current comparison
    sco <- sm[[cmp]]

    # Cycle through outputs
    for (out in colnames(sco)) {

      # Set margins
      graphics::par(mar = c(2, 4, 2, 4))

      # Get vector to plot, remove NAs
      pvals <- stats::na.omit(sco[, out])

      # Are any p-values left after filtering the NAs?
      if (length(pvals) > 0) {

        # Seems so, plot the p-values
        # Set bar plot parameters
        params <- list()
        params$height <- c(pvals)
        params$col <- pvalcol(pvals, c("darkgreen", "yellow", "red"))
        params$beside <- T
        params$main <- paste(cmp, out)
        params$las <- 2
        params$horiz <- TRUE
        do.call(get("barplot", asNamespace("graphics")), params)

      } else {

        # No, place an empty plot in this spot
        plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
        graphics::title(main = paste(cmp, out))

      }

    }
  }

  invisible(NULL)

}

#' Summary method for the assumptions of parametric tests used in multiple
#' comparisons of outputs
#'
#' Summary method for objects of class \code{assumptions_micomp}, which
#' contain the assumptions for the parametric tests used in multiple comparisons
#' of outputs.
#'
#' @param object Object of class \code{assumptions_micomp}.
#' @param ... Currently ignored.
#'
#' @return A list in which each component is associated with a distinct
#' comparison. Each component contains a matrix, in which columns represent
#' individual outputs and rows correspond to the statistical tests evaluating
#' the assumptions of the parametric tests used in each output. More
#' specifically, each matrix has rows with the following information:
#' \describe{
#'  \item{Royston (\emph{group}, \emph{ve=\%/npcs=})}{One row per group per
#'        variance to explain / number of PCs, with the \emph{p}-value yielded
#'        by the Royston test (\code{\link[MVN]{roystonTest}}) for the
#'        respective group and variance/npcs combination.}
#'  \item{Box's M (\emph{ve=\%/npcs=})}{One row per variance to explain with the
#'        \emph{p}-value yielded by Box's M test
#'        (\code{\link[biotools]{boxM}}).}
#'  \item{Shapiro-Wilk (\emph{group})}{One row per group, with the
#'        \emph{p}-value yielded by the Shapiro-Wilk test
#'        (\code{\link[stats]{shapiro.test}}) for the respective group.}
#'  \item{Bartlett}{One row with the \emph{p}-value yielded by Bartlett's test
#'        (\code{\link[stats]{bartlett.test}}).}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Create a micomp object, use provided dataset
#' mic <- micomp(5, c(0.7, 0.8, 0.9),
#'               list(list(name = "NLOKvsJEXOK", grpout = pphpc_ok),
#'                    list(name = "NLOKvsJEXNOSHUFF", grpout = pphpc_noshuff)),
#'               concat = TRUE)
#'
#' # Get the assumptions summary
#' sam <- summary(assumptions(mic))
#' }
summary.assumptions_micomp <- function(object, ...) {

  dims <- dim(object)
  nout <- dims[1]
  ncomp <- dims[2]
  cmpnames <- colnames(object)

  samic <- vector(mode = "list", length = ncomp)
  names(samic) <- cmpnames

  # How many variances to explain?
  ve_npcs <- attr(object, "ve_npcs")
  nve_npcs <- length(ve_npcs)

  # Cycle through comparisons
  for (i in 1:ncomp) {

    # Number of groups compared in current comparison
    grps <- names(object[[1, i]]$ttest$uvntest)
    ngrps <- length(grps)

    # Number of Royston tests performed: one per group per variance to explain
    # / number of PCs
    nroytests <- ngrps * nve_npcs

    # Allocate summary matrix for current comparison
    smat <- matrix(nrow = nroytests + nve_npcs + ngrps + 1, ncol = nout)
    colnames(smat) <- rownames(object)
    rnames <- NULL
    idx <- 0

    # Royston tests (ngrps) and Box's M tests (1) per requested variance to
    # explain
    for (j in 1:nve_npcs) {

      # Show variance to explain or number of PCs in info label?
      infolabel <-
        if (ve_npcs[j] < 1) {
          paste0("ve=", ve_npcs[j])
        } else {
          paste0("npcs=", ve_npcs[j])
        }

      # One Royston test per group
      for (k in 1:ngrps) {

        # Update summary matrix index
        idx <- idx + 1

        # Get Royston p-values
        smat[idx, ] <- sapply(object[, i], function(aco)
          if (length(aco$manova) >= j) {
            if (!is.null(aco$manova[[j]])) {
              if (methods::is(aco$manova[[j]]$mvntest[[k]], "royston")) {
                aco$manova[[j]]$mvntest[[k]]@p.value
              } else { NA }
            } else { NA }
          } else { NA })

        # Compose row name
        rnames <- c(rnames, paste0("Royston (", grps[k], ", ", infolabel, ")"))

      }

      idx <- idx + 1

      # One Box's M test per requested variance to explain
      smat[idx, ] <- sapply(object[, i], function(aco)
        if (length(aco$manova) >= j) {
          if (!is.null(aco$manova[[j]])) {
            aco$manova[[j]]$vartest$p.value
          } else { NA }
        } else { NA })

      # Compose Box's M test row name
      rnames <- c(rnames, paste0("Box's M (", infolabel, ")"))
    }

    # One Shapiro-Wilk test per group
    for (k in 1:ngrps) {

      # Update matrix row index
      idx <- idx + 1

      # Obtain p-values for the Shapiro-Wilk test for the 1st PC for the current
      # group
      smat[idx, ] <-
        sapply(object[, i], function(aco) aco$ttest$uvntest[[k]][[1]]$p.value)

      # Compose row name
      rnames <- c(rnames, paste0("Shapiro-Wilk (", grps[k], ")"))
    }

    # One Bartlett test per comparison
    idx <- idx + 1
    smat[idx, ] <-
      sapply(object[, i], function(aco) aco$ttest$vartest[[1]]$p.value)

    # Compose row name
    rnames <- c(rnames, "Bartlett")

    # Set summary matrix row names
    rownames(smat) <- c(rnames)

    samic[[cmpnames[i]]] <- smat

  }

  samic

}
