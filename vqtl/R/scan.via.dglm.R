#'  @title Conduct a Scanonevar Using the DGLM Function
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'  @description \code{scan.via.dglm} should not typically be called by a user.
#'  This function is used by both \code{scanonevar} and \code{scanonevar.perm}.
#'  This function is not typically necessary for a typical user.
#'
#'  @param mean.alt.formula The formula for the trait mean in the alternative model.
#'    \code{mean.null.formula} and \code{test.mean.effect} are inferred from it.
#'  @param var.alt.formula The formula for the trait variance in the alternative model.
#'    \code{var.null.formula} and \code{test.var.effect} are inferred from it.
#'  @param genoprobs The probability of each genotype for each individual.
#'  @param mapping.df The tbl_df with the response, all covariates, and space for the focal genotype.
#'  @param chr.by.marker a vector of the chromosome name of each marker
#'  @param pos.by.marker a vector of the position of each marker
#'  @param marker.names a vector of the name of each marker
#'  @param cor.threshold Numeric between 0 and 1 indicating how tightly a locus must be correlated with a covariate to be skipped.
#'    e.g. if cor.threshold is 0.8 (it's default) any locus with \code{cor(locus, covariate) > 0.8} will be skipped.
#'  @param perm The permutation to apply to the genotypes.  Defaults to identity permutation.
#'
#'  @inheritParams scanonevar
#'
#'  @return Returns a scanonevar object.
#'
#'  @seealso  \code{\link{scanonevar}}, \code{\link{scanonevar.perm}}
#'
#'  @details none
#'
#'

scan.via.dglm <- function(mean.alt.formula,
                          var.alt.formula,
                          genoprobs,
                          mapping.df,
                          chr.by.marker,
                          pos.by.marker,
                          marker.names,
                          return.effects = FALSE,
                          return.effect.ses = FALSE,
                          return.effect.ps = FALSE,
                          cor.threshold = 0.8,
                          perm = 1:nrow(genoprobs),
                          family = 'gaussian')
{

  # hack to get R CMD CHECK to run without NOTEs that these globals are undefined
  starts_with <- 'fake.global'

  # todo: consider tighening up the interface...do we really need 3 vectors '.by.marker'?

  #### figure out which tests to do based on where the QTL term appears ####
  # also make null formulae programmatically
  # todo: allow user to specity null formulae optionally
  test.mean.effect <- test.var.effect <- test.meanvar.effect <- FALSE
  mean.terms <- terms(mean.alt.formula)
  mean.qtl.terms <- grep(pattern = 'mean.QTL', x = attr(mean.terms, 'term.labels'))
  if (any(mean.qtl.terms)) {
    mean.null.terms <- drop.terms(termobj = mean.terms,
                                  dropx = mean.qtl.terms)
    mean.null.formula <- reformulate(termlabels = attr(mean.null.terms, 'term.labels'),
                                     response = mean.alt.formula[[2]])
    test.mean.effect <- TRUE
  }

  var.terms <- terms(var.alt.formula)
  var.qtl.terms <- grep(pattern = 'var.QTL', x = attr(var.terms, 'term.labels'))
  if (any(var.qtl.terms)) {
    var.null.terms <- drop.terms(termobj = var.terms,
                                 dropx = var.qtl.terms)
    var.null.formula <- reformulate(termlabels = attr(var.null.terms, 'term.labels'))
    test.var.effect <- TRUE
  }

  if (!any(test.mean.effect, test.var.effect)) {
    stop('No genetic effects specified (use (mean|var).QTL.(add|dom) in mean.formula or var.formula).')
  }
  test.meanvar.effect <- test.mean.effect & test.var.effect

  #### set up containers for the scan ####
  full.lod <- mean.lod <- var.lod <- rep(NA, length(marker.names))
  log10lik.bothalt <- log10lik.bothnull <- log10lik.meanalt <- log10lik.meannull <- log10lik.varalt <- log10lik.varnull <- NA
  chrtype <- rep(NA, length(marker.names))

  if (any(return.effects, return.effect.ses, return.effect.ps)) {
    num.mean.effects <- length(all.vars(mean.alt.formula)) # minus 1 for response, plus one for intercept term
    num.var.effects <- length(all.vars(var.alt.formula)) + 1 # add one for the intercept
    num.markers <- length(marker.names)
    mean.indices <- 1:(num.mean.effects)
    var.indices <- (num.mean.effects + 1):(num.mean.effects + num.var.effects)
  }
  if (return.effects) {
    fitted.effects <- array(data = NA, dim = c(num.markers, num.mean.effects + num.var.effects),              # add one for null marker
                            dimnames = list(marker.names,
                                            c('mean.intercept', paste0('mean.', all.vars(mean.alt.formula)[-1]),
                                              'var.intercept', paste0('var.', all.vars(var.alt.formula)))))
  }
  if (return.effect.ses) {
    effect.ses <- array(data = NA, dim = c(num.markers, num.mean.effects + num.var.effects),              # add one for null marker
                        dimnames = list(marker.names,
                                        c('se.mean.intercept', paste0('se.mean.', all.vars(mean.alt.formula)[-1]),
                                          'se.var.intercept', paste0('se.var.', all.vars(var.alt.formula)))))
  }
  if (return.effect.ps) {
    effect.ps <- array(data = NA, dim = c(num.markers, num.mean.effects + num.var.effects),              # add one for null marker
                       dimnames = list(marker.names,
                                       c('p.mean.intercept', paste0('p.mean.', all.vars(mean.alt.formula)[-1]),
                                         'p.var.intercept', paste0('p.var.', all.vars(var.alt.formula)))))
  }

  # calculate null for joint-test
  if (test.meanvar.effect) {
    both.null.fit <- tryCatch(dglm(formula = mean.null.formula,
                                   dformula = var.null.formula,
                                   data = mapping.df,
                                   family = family),
                              error = function(e) NA)
    log10lik.bothnull <- -0.5*both.null.fit$m2loglik / log(10)
  }

  # loop through markers (including pseudomarkers)
  for (marker.idx in 1:length(marker.names)) {

    # select focal marker from genoprob df
    marker.name <- marker.names[marker.idx]
    if (grepl(pattern = 'chr0X', marker.name)) {
      chrtype[marker.idx] <- 'X'
    } else {
      chrtype[marker.idx] <- 'A'
    }
    marker.genoprobs <- select(genoprobs, starts_with(paste0(marker.name, '_')))

    # JOINT TESTING
    if (ncol(marker.genoprobs) == 3) {
      mapping.df$mean.QTL.add <- mapping.df$var.QTL.add <- get.additive.coef.from.3.genoprobs(marker.genoprobs)[perm]
      mapping.df$mean.QTL.dom <- mapping.df$var.QTL.dom <- get.dom.coef.from.3.genoprobs(marker.genoprobs)[perm]
    }
    if (ncol(marker.genoprobs) == 2) {
      mapping.df$mean.QTL.add <- mapping.df$var.QTL.add <- get.additive.coef.from.2.genoprobs(marker.genoprobs)[perm]
      mapping.df$mean.QTL.dom <- mapping.df$var.QTL.dom <- 0
    }

    both.alt.fit <- tryCatch(dglm(formula = mean.alt.formula,
                                  dformula = var.alt.formula,
                                  data = mapping.df,
                                  family = family),
                             error = function(e) NA)

    log10lik.bothalt <- ifelse(test = identical(both.alt.fit, NA), yes = NA, no = -0.5*both.alt.fit$m2loglik / log(10))

    if (!identical(both.alt.fit, NA)) {

      if (any(return.effects, return.effect.ses, return.effect.ps)) {
        if (ncol(marker.genoprobs == 3)) {

          mean.indices <- 1:(num.mean.effects)
          var.indices <- (num.mean.effects + 1):(num.mean.effects + num.var.effects)
        }
        if (ncol(marker.genoprobs) == 2) {

          mean.idx.to.rm <- grep(pattern = 'dom', x = colnames(effect.ps[,mean.indices]))
          mean.indices <- mean.indices[-mean.idx.to.rm]

          var.idx.to.rm <- grep(pattern = 'dom', x = colnames(effect.ps[,var.indices]))
          var.indices <- var.indices[-var.idx.to.rm]
        }
      }
      if (return.effects) {
        fitted.effects[marker.idx, mean.indices] <- summary(both.alt.fit)$coef[,'Estimate']
        fitted.effects[marker.idx, var.indices] <- summary(both.alt.fit$dispersion.fit)$coef[,'Estimate']
      }
      if (return.effect.ses) {
        effect.ses[marker.idx , mean.indices] <- summary(both.alt.fit)$coef[,'Std. Error']
        effect.ses[marker.idx, var.indices] <- summary(both.alt.fit$dispersion.fit)$coef[,'Std. Error']
      }
      if (return.effect.ps) {
        effect.ps[marker.idx, mean.indices] <- summary(both.alt.fit)$coef[,'Pr(>|t|)']
        effect.ps[marker.idx, var.indices] <- summary(both.alt.fit$dispersion.fit)$coef[,'Pr(>|t|)']
      }
    }

    # store results
    full.lod[marker.idx] <- log10lik.bothalt - log10lik.bothnull

    # MEAN TESTING
    if (test.mean.effect) {

      # if we are not doing a permutation (doing a real scan) the mean alt is the same as the both alt
      if (identical(perm, 1:nrow(genoprobs))) {
        log10lik.meanalt <- log10lik.bothalt

      # if we are doing a permutation, we calculate a different mean alt
      } else {

        if (ncol(marker.genoprobs) == 3) {
          mapping.df$mean.QTL.add <- get.additive.coef.from.3.genoprobs(marker.genoprobs)[perm]
          mapping.df$var.QTL.add <- get.additive.coef.from.3.genoprobs(marker.genoprobs)
          mapping.df$mean.QTL.dom <- get.dom.coef.from.3.genoprobs(marker.genoprobs)[perm]
          mapping.df$var.QTL.dom <- get.dom.coef.from.3.genoprobs(marker.genoprobs)
        }
        if (ncol(marker.genoprobs) == 2) {
          mapping.df$mean.QTL.add <- get.additive.coef.from.2.genoprobs(marker.genoprobs)[perm]
          mapping.df$var.QTL.add <- get.additive.coef.from.2.genoprobs(marker.genoprobs)
          mapping.df$mean.QTL.dom <- mapping.df$var.QTL.dom <- 0
        }

        mean.alt.fit <- tryCatch(expr = dglm(formula = mean.alt.formula,
                                             dformula = var.alt.formula,
                                             data = mapping.df,
                                             family = family),
                                 error = function(e) NA)

        log10lik.meanalt <- ifelse(test = identical(mean.alt.fit, NA), yes = NA, no = -0.5*mean.alt.fit$m2loglik / log(10))
      }

      mean.null.fit <- tryCatch(expr = dglm(formula = mean.null.formula,
                                            dformula = var.alt.formula,
                                            data = mapping.df,
                                            family = family),
                                error = function(e) NA)

      log10lik.meannull <- ifelse(test = identical(mean.null.fit, NA), yes = NA, no = -0.5*mean.null.fit$m2loglik / log(10))

      mean.lod[marker.idx] <- log10lik.meanalt - log10lik.meannull
    }

    # VAR TESTING
    if (test.var.effect) {

      # if we are not doing a permutation (doing a real scan) the mean alt is the same as the both alt
      if (identical(perm, 1:nrow(genoprobs))) {
        log10lik.varalt <- log10lik.bothalt

      # if we are doing a permutation, we calculate a different mean alt
      } else {

        if (ncol(marker.genoprobs) == 3) {
          mapping.df$mean.QTL.add <- get.additive.coef.from.3.genoprobs(marker.genoprobs)
          mapping.df$var.QTL.add <- get.additive.coef.from.3.genoprobs(marker.genoprobs)[perm]
          mapping.df$mean.QTL.dom <- get.dom.coef.from.3.genoprobs(marker.genoprobs)
          mapping.df$var.QTL.dom <- get.dom.coef.from.3.genoprobs(marker.genoprobs)[perm]
        }
        if (ncol(marker.genoprobs) == 2) {
          mapping.df$mean.QTL.add <- get.additive.coef.from.2.genoprobs(marker.genoprobs)
          mapping.df$var.QTL.add <- get.additive.coef.from.2.genoprobs(marker.genoprobs)[perm]
          mapping.df$mean.QTL.dom <- mapping.df$var.QTL.dom <- 0
        }

        var.alt.fit <- tryCatch(expr = dglm(formula = mean.alt.formula,
                                            dformula = var.alt.formula,
                                            data = mapping.df,
                                            family = family),
                                error = function(e) NA)

        log10lik.varalt <-  ifelse(test = identical(var.alt.fit, NA), yes = NA, no = -0.5*var.alt.fit$m2loglik / log(10))
      }

      var.null.fit <- tryCatch(expr = dglm(formula = mean.alt.formula,
                                           dformula = var.null.formula,
                                           data = mapping.df),
                               error = function(e) NA)

      log10lik.varnull <- ifelse(test = identical(var.null.fit, NA), yes = NA, no = -0.5*var.null.fit$m2loglik / log(10))

      var.lod[marker.idx] <- log10lik.varalt - log10lik.varnull
    }

  }

  # compile data into return format
  varscan <- data_frame(chr = factor(chr.by.marker,
                                     levels = mixedsort(unique(chr.by.marker))),
                        chrtype = chrtype,
                        pos = pos.by.marker,
                        marker.name = marker.names,
                        full.lod,
                        mean.lod,
                        var.lod)
  if (return.effects) {
    varscan <- bind_cols(varscan, data.frame(fitted.effects))
  }
  if (return.effect.ses) {
    varscan <- bind_cols(varscan, data.frame(effect.ses))
  }
  if (return.effect.ps) {
    varscan <- bind_cols(varscan, data.frame(effect.ps))
  }

  class(varscan) <- c('scanonevar', class(varscan))
  attr(varscan, 'pheno') <- as.character(mean.null.formula[[2]])
  attr(varscan, 'units') <- 'lods'
  if (test.mean.effect) { attr(varscan, 'mean.null.formula') <- mean.null.formula }
  if (test.var.effect) { attr(varscan, 'var.null.formula') <- var.null.formula }
  attr(varscan, 'mean.alt.formula') <- mean.alt.formula
  attr(varscan, 'var.alt.formula') <- var.alt.formula
  attr(varscan, 'null.fit') <- both.null.fit

  return(varscan)
}