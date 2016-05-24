#'  @title Plot Phenotype of interest Averaged (Marginalized) Across Specified Markers and Phenotypes
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'  @description \code{margin.plot} should be used to visually investigate the relationship
#'    between the phenotype of interest and other phenotypes.  \code{margin.plot} can also
#'    be used to visualize the relationship between the phenotype of interest and genetic
#'    loci of interest, but \code{predictive.plot} is usually preferrable.
#'
#'  @param cross The cross object to be plotted
#'  @param focal.phenotype.name the phenotype to put on the y-axis
#'  @param marginal.phen.names a list of phenotypes to average over (put on the x-axis).
#'  @param marginal.marker.names a list of marker names, whose values will be averaged over (put on the x-axis).
#'  @param genotype.plotting.names Labels for the genotype groups.  Defaults to \code{c('AA', 'AB', 'BB')}.
#'  @param subset the subset of individuals to use
#'  @param col optionally, color of dots, as in base R graphics.  Defaults to gray.
#'  @param pch optionally, plotting character, as in base R graphics.  Defaults to 19 (disc).
#'  @param xlab.override optionally, x axis label, as in base R graphics.  Defaults to the name of the marginal marker.
#'  @param ylab.override optionally, y axis label, as in base R graphics.  Defaults to focal phenotype name.
#'  @param title.override optionally, plot title, as in base R graphics.  Defaults to 'focal phenotype name by marginal phenotype name'.
#'  @param title.cex optionally, character expansion for title, as in base R graphics.  Defaults to 1.5.
#'  @param circle.alpha optionally, alpha (transparency) of discs.  Defaults to 0.2.
#'
#'  @return None.  Only makes plot.
#'
#'  @details none
#'
#'  @examples
#'    set.seed(27599)
#'    my.cross <- sim.cross(map = sim.map(), type = 'f2')
#'    my.cross$pheno$phenotype <- rnorm(n = 100,
#'                                      mean = my.cross$geno$`1`$data[,5],
#'                                      sd = my.cross$geno$`2`$data[,5])
#'    my.cross$pheno$sex <- rbinom(n = 100, size = 1, prob = 0.5)
#'    my.cross$pheno$cage <- sample(x = 1:5, size = 100, replace = TRUE)
#'
#'    margin.plot(cross = my.cross,
#'                focal.phenotype.name = 'phenotype',
#'                marginal.phen.name = list('sex', 'cage'),
#'                marginal.marker.name = list('D1M5', 'D2M5'))
#'
#'
margin.plot <- function(cross,
                        focal.phenotype.name,
                        marginal.phen.names = NULL,
                        marginal.marker.names = NULL,
                        genotype.plotting.names = c('A', 'H', 'B'),
                        subset = 1:nind(cross),
                        col = rep(rgb(0.5, 0.5, 0.5, 0.5), nind(cross)),
                        pch = 19,
                        xlab.override = NA,
                        ylab.override = NA,
                        title.override = NA,
                        title.cex = 1.5,
                        circle.alpha = 0.2) {

  if (any(missing(cross), missing(focal.phenotype.name), !(focal.phenotype.name %in% names(cross$pheno)))) {
    stop('Must provide a cross and a focal phenotype in that cross.')
  }

  num.plots <- sum(length(marginal.phen.names), length(marginal.marker.names))
  if (num.plots == 0) { stop('Must provide a marginal phenotype or marker.')}

  focal.phen <- cross$pheno[[focal.phenotype.name]][subset]
  if (!missing(col) & length(col) == nind(cross)) {
    col <- col[subset]
  }

  for (marginal.phen.name in marginal.phen.names) {

    marginal.phen <- factor(cross$pheno[[marginal.phen.name]][subset])
    plotting.phen <- as.numeric(marginal.phen)

    plot(x = jitter(plotting.phen),
         y = focal.phen,
         xaxt = 'n',
         xlab = NA,
         ylab = NA,
         axes = FALSE,
         col = alpha(col, circle.alpha),
         pch = pch)

    # x axis stuff
    mtext(side = 1, text = ifelse(test = is.na(xlab.override),
                                  yes = marginal.phen.name,
                                  no = xlab.override),
          line = 2)
    mtext(side = 1, at = 1:length(levels(marginal.phen)), text = levels(marginal.phen))

    # y axis stuff
    axis(side = 2)
    mtext(side = 2, text = ifelse(test = is.na(ylab.override),
                                  yes = focal.phenotype.name,
                                  no = ylab.override),
          line = 2)

    # plot title
    mtext(side = 3, text = ifelse(test = is.na(title.override),
                                  yes = paste(focal.phenotype.name, 'by', marginal.phen.name),
                                  no = title.override),
          line = 1, cex = title.cex)
  }

  for (marginal.marker.name in marginal.marker.names) {

    chr.of.interest <- which(sapply(X = cross$geno, FUN = function(chr) { marginal.marker.name %in% colnames(chr$data)}))

    genotypes <- cross$geno[[chr.of.interest]]$data[subset,marginal.marker.name]
    plot(x = jitter(genotypes),
         y = focal.phen,
         xaxt = 'n',
         xlab = NA,
         ylab = NA,
         axes = FALSE,
         col = alpha(col, circle.alpha),
         pch = pch)
    axis(side = 2)
    mtext(text = genotype.plotting.names, side = 1, line = 0, at = 1:3)
    mtext(side = 1, text = ifelse(test = is.na(xlab.override),
                                  yes = marginal.marker.name,
                                  no = xlab.override),
          line = 2)
    mtext(side = 2, text = ifelse(test = is.na(ylab.override),
                                  yes = focal.phenotype.name,
                                  no = ylab.override),
          line = 2)
    mtext(side = 3, text = ifelse(test = is.na(title.override),
                                  yes = paste(focal.phenotype.name, 'by', marginal.marker.name),
                                  no = title.override),
          line = 1, cex = title.cex)

    means <- aggregate(x = focal.phen, by = list(genotypes), FUN = mean)[,2]
    sds <- aggregate(x = focal.phen, by = list(genotypes), FUN = sd)[,2]

    x.start.wide <- 1:3 - 0.2
    x.end.wide <- 1:3 + 0.2
    x.start.narrow <- 1:3 - 0.1
    x.end.narrow <- 1:3 + 0.1

    # horizontal lines at means and means +/- 1SD
    segments(x0 = c(x.start.narrow, x.start.wide, x.start.wide),
             y0 = c(means, means - sds, means + sds),
             x1 = c(x.end.narrow, x.end.wide, x.end.wide),
             y1 = c(means, means - sds, means + sds),
             lwd = rep(c(4, 2, 2), each = 3))

    # vertical line down the middle of each genotype group
    segments(x0 = 1:3,
             y0 = means - sds,
             x1 = 1:3,
             y1 = means + sds,
             lwd = 2)

    # dotted horizontal lines connecting means and SD's
    segments(x0 = rep(1:2, 3),
             y0 = c(means[1:2], means[1:2] + sds[1:2], means[1:2] - sds[1:2]),
             x1 = rep(2:3, 3),
             y1 = c(means[2:3], means[2:3] + sds[2:3], means[2:3] - sds[2:3]),
             lty = 2)
  }

  # reset graphical parameteers to how they were on start
  # par(start.pars)

  # return nothing
  invisible()
}