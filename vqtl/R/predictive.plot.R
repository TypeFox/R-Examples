#'  @title Plot Predictive Interval for Categorical Genotype/Phenotype Groups
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'  @description \code{predictive.plot} should be used to visually investigate loci identified
#'    with plot.scanonevar or summary.scanonevar.  The user can specify the same mean and variance
#'    formulae that were used in the scan, or specify new formulae to investigate interactions.
#'
#'  @param cross The cross object to be plotted
#'  @param marker.name The name of the marker the effects of which we want to investigate and visualize.
#'  @param phen.name The categorical phenotype the effects of which we want to investigate and visualize.
#'  @param genotype.plotting.names Labels for the genotype groups.  Defaults to \code{c('AA', 'AB', 'BB')}.
#'  @param title Optionally, title for the plot.  Defaults to 'Predictive of [response phenotype] from
#'    [predictive phenotype (e.g. sex)] and [marker name]
#'  @param title.cex Optionally, character expansion for title.  Defaults to 1.
#'  @param ribbon.width Optionally, width of ribbon connecting same-phenotype (different genotype) groups.
#'    Defaults to 10.
#'  @param xlim Optionally specify x-axis limits.  Defaults to data-dependent.
#'  @param ylim Optionally specify y-axis limits.  Defaults to data.dependent.
#'
#'
#'  @return None.  Only makes plot.
#'
#'  @inheritParams scanonevar
#'
#'  @details none
#'
#'  @examples
#'
#'    set.seed(27599)
#'    my.cross <- sim.cross(map = sim.map(), type = 'f2')
#'    my.cross <- calc.genoprob(my.cross)
#'    my.cross$pheno$phenotype <- rnorm(n = 100,
#'                                      mean = my.cross$geno$`1`$data[,5],
#'                                      sd = my.cross$geno$`2`$data[,5])
#'    my.cross$pheno$sex <- rbinom(n = 100, size = 1, prob = 0.5)
#'    my.cross$pheno$cage <- sample(x = 1:5, size = 100, replace = TRUE)
#'
#'    predictive.plot(cross = my.cross,
#'                    mean.formula = 'phenotype ~ sex + mean.QTL.add + mean.QTL.dom',
#'                    var.formula = '~ sex + var.QTL.add + var.QTL.dom',
#'                    marker.name = 'D1M5',
#'                    phen.name = 'sex')
#'
#'    predictive.plot(cross = my.cross,
#'                    mean.formula = 'phenotype ~ sex + mean.QTL.add + mean.QTL.dom',
#'                    var.formula = '~ sex + var.QTL.add + var.QTL.dom',
#'                    marker.name = 'D2M5',
#'                    phen.name = 'sex')
#'
predictive.plot <- function(cross,
                            mean.formula,
                            var.formula,
                            marker.name,
                            phen.name,
                            title = paste('Predictive of', response.phen, 'from', phen.name, 'and', marker.name),
                            title.cex = 1,
                            genotype.plotting.names = c('AA', 'AB', 'BB'),
                            ribbon.width = 10,
                            xlim = NA,
                            ylim = NA) {

  # hack to get R CMD CHECK to run without NOTEs that these globals are undefined
  genotype <- plotting.genotype <- 'fake.global'
  indiv.mean.estim <- indiv.mean.lb <- indiv.mean.ub  <- 'fake.global'
  indiv.var.estim <- indiv.var.lb <- indiv.var.ub  <- 'fake.global'


  # store current graphical parameters and customize them for this plot
  # start.pars <- par(no.readonly = TRUE)
  # par(mar = c(2, 3, 6, 2))

  mean.formula <- formula(mean.formula)
  var.formula <- formula(var.formula)

  phenotypes <- cross$pheno[[phen.name]]
  genoprobs <- get.genoprobs.by.marker.name(cross = cross, marker.name = marker.name)
  genotypes <- get.genotypes.by.marker.name(cross = cross, marker.name = marker.name, as.matrix = FALSE)
  response.phen <- mean.formula[[2]]

  # set up display names for the genotype groups (rather than default AA, AB, BB)
#   if (length(genotype.plotting.names) != length(unique(genotypes))) {
#     stop('length of genotype.plotting.names must be equal to the number of genotypes at marker.name')
#   }
  plotting.genotypes <- mapvalues(x = genotypes, from = sort(unique(genotypes)), to = genotype.plotting.names)

  # TODO: consider a similar 'plotting' version of phenotype group names?

  if (dim(genoprobs)[2] == 3) {
    model.df <- data.frame(mean.QTL.add = get.additive.coef.from.3.genoprobs(genoprobs),
                           mean.QTL.dom = get.dom.coef.from.3.genoprobs(genoprobs),
                           var.QTL.add = get.additive.coef.from.3.genoprobs(genoprobs),
                           var.QTL.dom = get.dom.coef.from.3.genoprobs(genoprobs))
  } else {
    model.df <- data.frame(mean.QTL.add = get.additive.coef.from.2.genoprobs(genoprobs),
                           var.QTL.add = get.additive.coef.from.2.genoprobs(genoprobs))
  }


  model.df <- cbind(model.df, cross$pheno)

  # if no 'var.formula', use lm for modeling
  if (missing(var.formula)) {

    lm.fit <- lm(formula = mean.formula, data = model.df)

    mean.pred <- predict(object = lm.fit, se.fit = TRUE)
    mean.estim <- mean.pred$fit
    mean.se <- mean.pred$se.fit

    if (missing(var.formula)) {
      var.estim <- mean.pred$residual.scale
      var.se <- 0

      prediction.tbl <- data_frame(genotype = genotypes,
                                   plotting.genotype = plotting.genotypes,
                                   phen = phenotypes,
                                   indiv.mean.estim = mean.estim,
                                   indiv.mean.lb = mean.estim - mean.se,
                                   indiv.mean.ub = mean.estim + mean.se,
                                   indiv.var.estim = var.estim,
                                   indiv.var.lb = var.estim - var.se,
                                   indiv.var.ub = var.estim + var.se)
    }
  }

  # if there is a 'var.formula' specified, use it in the DGLM model
  if (!missing(var.formula)) {
    dglm.fit <- dglm(formula = mean.formula,
                     dformula = var.formula,
                     data = model.df)

    mean.pred <- predict(dglm.fit, se.fit = TRUE)
    mean.estim <- mean.pred$fit
    mean.se <- mean.pred$se.fit

    var.pred <- predict(dglm.fit$dispersion.fit, se.fit = TRUE)
    var.estim <- var.pred$fit/var.pred$residual.scale
    var.se <- var.pred$se.fit/var.pred$residual.scale

    prediction.tbl <- data_frame(genotype = genotypes,
                                 plotting.genotype = plotting.genotypes,
                                 phen = phenotypes,
                                 indiv.mean.estim = mean.estim,
                                 indiv.mean.lb = mean.estim - mean.se,
                                 indiv.mean.ub = mean.estim + mean.se,
                                 indiv.var.estim = exp(var.estim),
                                 indiv.var.lb = exp(var.estim - var.se),
                                 indiv.var.ub = exp(var.estim + var.se))
  }



  plotting.tbl <- prediction.tbl %>%
    group_by(phen, genotype, plotting.genotype) %>%
    summarise(group.mean.estim = mean(indiv.mean.estim),
              group.mean.lb = mean(indiv.mean.lb),
              group.mean.ub = mean(indiv.mean.ub),
              group.var.estim = mean(indiv.var.estim),
              group.var.lb = mean(indiv.var.lb),
              group.var.ub = mean(indiv.var.ub)) %>%
    arrange(phen, genotype)

  # set up colors -- special red/blue is phen is 'sex'
  unique.phens <- sort(unique(plotting.tbl$phen))
  num.unique.phens <- length(unique.phens)
  if (phen.name == 'sex') {
    phen.group.colors <- c('red', 'blue')
  } else {
    phen.group.colors <- brewer.pal(n = num.unique.phens, name = 'Set1')
  }
  plotting.tbl$phen.col <- mapvalues(x = plotting.tbl$phen, from = unique.phens, to = phen.group.colors)

  # blank plot of correct size
  if (any(is.na(xlim))) { xlim <- range(c(plotting.tbl$group.mean.lb, plotting.tbl$group.mean.ub)) }
  if (any(is.na(ylim))) { ylim <- range(c(plotting.tbl$group.var.lb, plotting.tbl$group.var.ub)) }
  with(plotting.tbl,
       plot(x = 1, y = 1,
            type = 'n',
            xlab = NA,
            ylab = NA,
            axes = FALSE,
            xlim = xlim,
            ylim = ylim))
  axis(side = 1)
  axis(side = 2)
  mtext(side = 3, line = 0, cex = title.cex, text = title)
  mtext(side = 1, text = paste('phenotype mean'), line = 2)
  mtext(side = 2, text = paste('phenotype SD'), line = 2)


  # horizontal lines
  with(plotting.tbl,
       segments(x0 = group.mean.lb,
                y0 = group.var.estim,
                x1 = group.mean.ub,
                y1 = group.var.estim,
                col = alpha(phen.col, 0.5),
                lwd = 3))

  # vertical lines
  with(plotting.tbl,
       segments(x0 = group.mean.estim,
                y0 = group.var.lb,
                x1 = group.mean.estim,
                y1 = group.var.ub,
                col = alpha(phen.col, 0.5),
                lwd = 3))

  # draw light lines connecting same-phenotype groups
  # e.g. connect male AA, male AB, and male BB together (same for female)
  for (phen.group.idx in 1:num.unique.phens) {
    phen <- unique.phens[phen.group.idx]
    col <- col2rgb(phen.group.colors[phen.group.idx])

    this.group.idxs <- plotting.tbl$phen == phen

    lines(x = plotting.tbl$group.mean.estim[this.group.idxs],
          y = plotting.tbl$group.var.estim[this.group.idxs],
          col = rgb(red = col[1], green = col[2], blue = col[3], maxColorValue = 255, alpha = 50),
          lwd = ribbon.width)
  }

  # white circles clear space to write genotype names
  with(plotting.tbl,
       points(x = group.mean.estim, y = group.var.estim,
              pch = 19, col = rgb(1, 1, 1, 0.5), cex = 5))

  # write genotype names
  with(plotting.tbl,
       text(x = group.mean.estim, y = group.var.estim,
            labels = plotting.genotype, cex = 1.5,
            col = phen.col))

  # reset graphical parameteers to how they were on start
#   par(start.pars)

  # return nothing
  invisible()
}
