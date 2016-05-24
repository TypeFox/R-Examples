#' Plot carrying capacity and abundance trajectories
#' 
#' Plot each population's carrying capacity and abundance over time.
#' 
#' @param meta The R object holding population info returned by
#'   \code{\link{meta}}.
#' @param kch The R object holding kch data returned by \code{\link{kch}}.
#' @param pops (Optional) A character vector of population names. The
#'   metapopulation will be subset to these populations before plotting. If not
#'   provided, all populations will be plotted (see Note below).
#' @param samelims (logical) If \code{TRUE}, the y-axis limits will be constant
#'   across plots.
#' @param show_N (logical) If \code{TRUE}, mean population abundance will be
#'   plotted (solid lines) in addition to carrying capacity (dashed lines).
#' @param results (required only if \code{show_N} is \code{TRUE}) The R object
#'   holding simulation results returned by \code{\link{results}}.
#' @param ... Additional arguments to \code{lattice::\link[lattice]{xyplot}},
#'   e.g., \code{layout}.
#' @return A lattice object is returned invisibly, and plotted if not assigned.
#' @importFrom zoo zoo
#' @importFrom lattice xyplot panel.rect panel.text
#' @export
#' @note When plotting many populations, \code{layout} should be set 
#'   appropriately, and it may be useful to plot to, e.g., a \code{pdf} device. 
#' @examples
#' mp <- system.file('example.mp', package='mptools')
#' met <- meta(mp)
#' 
#' # Subset of populations
#' knt(met, pops=c('Pop 169', 'Pop 170', 'Pop 174', 'Pop 175'), 
#'     kch(met, dirname(mp)), show_N=TRUE, results=results(mp),
#'     layout=c(2, 2), samelims=TRUE)
knt <- function(meta, kch, pops, samelims=FALSE, show_N=FALSE, results, ...) {
  if(!missing(pops)) {
    if(any(!pops %in% meta$popName)) {
      warning('Some populations not in meta:\n', setdiff(pops, meta$popName),
              call.=FALSE)
    }
    meta <- meta[meta$popName %in% pops, ]
    kch <- kch[, intersect(pops, colnames(kch))]
    if(any(!pops %in% colnames(kch))) {
      warning('Some populations not in kch:\n', setdiff(pops, colnames(kch)),
              call.=FALSE)
    }
  }
  hasInitN <- meta$initN > 0
  Strip <- function(which.panel, factor.levels, ...) {
    lattice::panel.rect(
      0, 0, 1, 1, border = 1,
      col=ifelse(hasInitN, 'lightsteelblue2', 'gray90')[which.panel])
    lattice::panel.text(x = 0.5, y = 0.5, lab = factor.levels[which.panel],
                        font = ifelse(hasInitN, 2, 1)[which.panel])
  }
  scl <- list(alternating=FALSE, tck=c(1, 0), cex=0.8, rot=0)
  if (isTRUE(samelims)) {
    scl$y <- list(relation = "same")
  } else {
    scl$y <- list(limits=mapply(c, 0, apply(kch, 2, max), SIMPLIFY=FALSE))
  }
  if (isTRUE(show_N)) {
    if(!missing(pops)) {
      if(any(!pops %in% dimnames(results$results)[[3]])) {
        warning('Some populations not in results:\n', 
                setdiff(pops, dimnames(results$results)[[3]]),
                call.=FALSE)
      }
      results$results <- results$results[,, c('ALL', pops)]
    }
    kchn <- cbind(kch, results$results[, 'mean', -1]
                [, match(colnames(results$results[, 'mean', -1]), 
                         colnames(kch))])
    colnames(kchn) <- make.unique(colnames(kchn))
    p <- lattice::xyplot(zoo::zoo(kchn), screens=rep(colnames(kch), 2),
                         lwd=rep(c(3, 1), each=ncol(kch)),
                         scales=scl, ylab='Population size', 
                         strip=Strip, col=rep(c('gray85', 1), each=ncol(kch)), 
                         ...)
  } else {
    p <- lattice::xyplot(zoo::zoo(kch), scales=scl, 
                         ylab='Carrying capacity', strip=Strip, col=1, ...)
  }
  p
}
