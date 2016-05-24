#' Various DeLorean object plots
#'
#' @param x de.lorean object
#' @param type Type of plot:
#'   \itemize{
#'     \item 'expr.data': The expression data plotted by capture time.
#'       See \code{\link{expr.data.plot}}.
#'     \item 'Rhat': \eqn{hat{R}} convergence statistics
#'       See \code{\link{Rhat.plot}}.
#'     \item 'pseudotime': Pseudotimes in best posterior sample
#'       See \code{\link{pseudotime.plot}}.
#'     \item 'profiles': Gene expression profiles for best posterior sample
#'       See \code{\link{profiles.plot}}.
#'     \item 'tau.offsets': Offsets of pseudotimes to assess the prior
#'       See \code{\link{tau.offsets.plot}}.
#'     \item 'marg.like': Plot the posterior of the marginal likelihoods
#'       for individual genes.
#'       See \code{\link{marg.like.plot}}.
#'     \item 'roughnesses': Roughnesses of the pseudotime posterior
#'       See \code{\link{roughnesses.plot}}.
#'   }
#' @param ... Extra arguments to plot function
#'
#' @method plot de.lorean
#' @export
#'
plot.de.lorean <- function(x, type="profiles", ...) {
    result <- switch(type,
        profiles=profiles.plot(x, ...),
        pseudotime=pseudotime.plot(x, ...),
        Rhat=Rhat.plot(x, ...),
        expr.data=expr.data.plot(x, ...),
        roughnesses=roughnesses.plot(x, ...),
        marg.like=marg.like.plot(x, ...),
        orderings=orderings.plot(x, ...),
        tau.offsets=tau.offsets.plot(x, ...)
    )
    if (is.null(result)) {
        stop('Unknown plot type')
    }
    result
}

#' Calculate a suitable value for a rug plot given the
#' number of points
#'
#' @param n Number of points.
#' @param scale Scale the value.
#'
alpha.for.rug <- function(n, scale=100) {
    1 / (max(1, n / scale))
}

#' Plot posterior for marginal log likelihoods of individual gene's
#' expression profiles
#'
#' @param dl de.lorean object
#'
#' @export
#'
marg.like.plot <- function(dl) {
    with(dl, {
        gp <- (ggplot(samples.l$logmarglike %>% left_join(gene.map),
                      aes(x=gene,
                          y=logmarglike,
                          colour=is.held.out),
                      environment=environment())
            + geom_boxplot()
        )
    })
}

#' Plot pseudotime (tau) against observed capture time.
#'
#' @param dl de.lorean object
#' @param sample.iter Which sample to take pseudotimes from
#'
#' @export
#'
pseudotime.plot <- function(dl, sample.iter=dl$best.sample) {
    with(dl, {
        gp <- (ggplot(samples.l$tau %>% filter(iter == sample.iter),
                      aes(x=tau, y=obstime, color=capture),
                      environment=environment())
            + geom_point()
            + geom_vline(data=cell.meta %>% group_by(capture),
                         aes(xintercept=obstime, color=capture),
                         linetype=2,
                         alpha=.8)
            + scale_x_continuous(name="Pseudotime (tau)")
            + scale_y_continuous(name="Observed (capture) time")
        )
    })
}


#' Plot the tau offsets, that is how much the pseudotimes (tau) differ
#' from their prior means over the full posterior.
#'
#' @param dl de.lorean object
#' @param rug.alpha Alpha parameter for rug geom
#'
#' @export
#'
tau.offsets.plot <- function(dl, rug.alpha=.3) {
    with(dl,
         ggplot(samples.l$tau, aes(x=tau.offset, color=capture))
         + geom_density()
         + geom_rug(alpha=rug.alpha)
         + stat_function(fun=functional::Curry(dnorm, sd=hyper$sigma_tau),
                         linetype='dashed',
                         alpha=.7,
                         color='blue')
    )
}


#' Plot the Rhat convergence statistics. \code{\link{examine.convergence}}
#' must be called before this plot can be made.
#'
#' @param dl de.lorean object
#'
#' @export
#'
Rhat.plot <- function(dl) {
    with(dl, {
        rhat.df <- data.frame(
            rhat=rhat.sorted,
            param=names(rhat.sorted),
            parameter=stringr::str_match(names(rhat.sorted), "^[[:alpha:]]+"))
        gp <- (ggplot(rhat.df,
                      aes(y=rhat, x=parameter),
                      environment=environment())
            + geom_boxplot()
        )
    })
}

#' Plot a comparison of the profiles from several de.lorean objects
#'
#' @param ... Named de.lorean objects
#' @param genes Genes to plot (defaults to genes.high.psi of first de.lorean
#'   object)
#'
#' @export
#'
cmp.profiles.plot <- function(..., genes = NULL) {
    dls <- list(...)
    dl.levels <- names(dls)
    stopifnot(! is.null(dl.levels))  # Must have names for de.lorean objects
    if (is.null(genes)) {
        genes <- dls[[1]]$genes.high.psi
    }
    get.mean <- function(.name) {
        with(dls[[.name]], (
            predictions
            %>% filter(best.sample == iter)
            %>% left_join(gene.map)
            %>% filter(gene %in% genes)
            %>% mutate(name=factor(.name, levels=dl.levels))
        ))
    }
    means <- do.call(rbind, lapply(dl.levels, get.mean))
    gp <- ggplot(mutate.profile.data(means),
                 aes(x=tau),
                 environment=environment())
    line.alpha <- .8
    ribbon.alpha <- .2
    (
        gp
        + geom_line(aes(x=x, y=mean, color=name),
                    alpha=line.alpha)
        + geom_ribbon(aes(x=x,
                          ymin=mean-2*sqrt(var),
                          ymax=mean+2*sqrt(var),
                          fill=name),
                      alpha=ribbon.alpha)
        + facet_wrap(~ gene)
        + scale_x_continuous(name="Pseudotime",
                             breaks=unique(dls[[1]]$cell.meta$obstime))
        + scale_y_continuous(name="Expression")
    )
}

#' Plot best sample predicted expression.
#'
#' @param dl de.lorean object
#' @param genes Genes to plot (defaults to genes.high.psi)
#' @param profile.color Colour for the profile
#' @param add.data Add actual expression data to plot
#' @param sample.iter Which sample to plot
#' @param ... Extra arguments
#'
#' @export
#'
profiles.plot <- function(dl,
                          genes=dl$genes.high.psi,
                          profile.color='black',
                          add.data=T,
                          sample.iter=dl$best.sample,
                          ...) {
    varargs <- list(...)
    with(dl, {
        if (opts$periodic) {
            modulo.period <- function(t) ( t - floor(t / opts$period)
                                                * opts$period )
        } else {
            modulo.period <- function(t) { t }
        }
        gp <- (ggplot(predictions
                      %>% filter(sample.iter == iter)
                      %>% left_join(gene.map)
                      %>% filter(gene %in% genes),
                      environment=environment()))
        profile.data <- (
            predictions
            %>% filter(sample.iter == iter)
            %>% left_join(dl$gene.map)
            %>% filter(gene %in% genes)
        )
        # stopifnot(! any(is.na(profile.data %>% select(-cbRank, -cbPeaktime))))
        gp <- (
            plot.add.mean.and.variance(
                gp,
                .data=mutate.profile.data(profile.data),
                color=profile.color)
            + facet_wrap(~ gene)
            + scale_x_continuous(name="Pseudotime",
                                 breaks=unique(cell.meta$obstime))
            + scale_y_continuous(name="Expression")
        )
        if (add.data) {
            expr.data <- (
                gene.map
                %>% filter(gene %in% genes)
                %>% left_join(melt(unname(expr),
                                   varnames=c("g", "c"),
                                   value.name="expr"))
                %>% left_join(samples.l$tau
                              %>% filter(sample.iter == iter)
                              %>% mutate(tau=modulo.period(tau))))
            if (! is.null(varargs$cell.size.adj) && varargs$cell.size.adj) {
                expr.data <- (
                    expr.data
                    %>% left_join(samples.l$S)
                    %>% mutate(expr=expr - S))
            }
            gp <- plot.add.expr(gp, .data=expr.data)
        }
        gp
    })
}

#' Mutate the profile data into shape compatible with GP plot function
#'
#' @param .data The data
#'
mutate.profile.data <- function(.data) {
    (
        .data
        %>% mutate(x=tau, mean=predictedmean+phi, var=predictedvar)
        %>% dplyr::select(-tau, -predictedmean, -phi, -predictedvar)
    )
}


# Adjust the predicted mean with the predictions from the model.
#
adjust.predictions <- function(.data, adjust.model) {
    # print(names(.data))
    adjustments <- (
        .data
        %>% group_by(t)
        %>% dplyr::summarise(tau=tau[1]))
    # print(tail(adjustments))
    adjustments$adjustment <- predict(adjust.model,
                                        newdata=adjustments)
    # .T <- nrow(adjustments)
    # print(adjustments$adjustment[1:(.T-1)]-adjustments$adjustment[2:.T])
    (
        .data
        %>% left_join(dplyr::select(adjustments, -tau))
        %>% mutate(predictedmean=predictedmean+adjustment))
}

#' Add expression data to a plot
#'
#' @param gp Plot object
#' @param .data Expression data to add
#'
plot.add.expr <- function(gp, .data=NULL)
{
    (gp + geom_point(data=.data,
                     aes(x=tau,
                         y=expr,
                         color=capture),
                     size=4,
                     alpha=.7))
}


#' Plot the expression data by the capture points
#'
#' @param dl de.lorean object
#' @param genes Genes to plot. If NULL plots some random varying genes
#' @param num.genes Number of genes to plot
#'
#' @export
#'
expr.data.plot <- function(dl, genes=NULL, num.genes=12) {
    with(dl, {
         if (is.null(genes)) {
             num.to.sample <- min(nrow(expr), num.genes * 10)
             sample.genes <- sample(rownames(expr), num.to.sample)
             expr.l <- (
                 expr[sample.genes,]
                 %>% melt(varnames=c("gene", "cell"), value.name="x"))
             variation <- (
                 expr.l
                 %>% group_by(gene)
                 %>% dplyr::summarise(var=var(x))
                 %>% arrange(-var))
             if (nrow(variation) > num.genes) {
                 variation <- variation %>% head(num.genes)
             }
             genes <- variation$gene
         } else {
            expr.l <- (
                expr[genes,]
                %>% melt(varnames=c("gene", "cell"), value.name="x"))
         }
         stopifnot(all(genes %in% rownames(expr)))
         expr.l <- expr.l %>%
             mutate(cell=factor(cell, levels=levels(cell.meta$cell))) %>%
             left_join(cell.meta) %>%
             dplyr::filter(gene %in% genes)
         ggplot(expr.l, aes(x=capture, y=x)) +
             # geom_boxplot() +
             geom_violin() +
             stat_summary(fun.y=mean,
                          colour="red",
                          aes(group=gene),
                          geom="line") +
             facet_wrap(~ gene)
    })
}

#' Plot two sets of pseudotimes against each other.
#'
#' @param dl The DeLorean object
#' @param fits Fit indexes
#'
#' @export
#'
pseudotimes.pair.plot <- function(dl, fits=NULL) {
  stopifnot(all(dim(dl$best.orderings) == dim(dl$vb.tau)))
  stopifnot(! is.null(dl$vb.tau))  # Must have estimated tau.
  #
  # Create a data frame with the best orderings
  best.o <- sapply(dl$best.orderings, function(O) O$ser.order)
  best.o.m <- reshape2::melt(best.o, varnames=c('c', 'fit'),
                             value.name='ordering.idx') %>%
              dplyr::left_join(data.frame(
                ordering.idx=1:dl$stan.data$C,
                ordering=even.tau.spread(dl)))
  #
  # Create a data frame with the best pseudotimes
  best.tau.m <- reshape2::melt(dl$vb.tau, varnames=c('c', 'fit'),
                               value.name='pseudotime')
  #
  # If not plotting all fits then filter data frames
  if (! is.null(fits)) {
    best.o.m <- filter(best.o.m, fit %in% fits)
    best.tau.m <- filter(best.tau.m, fit %in% fits)
  }
  #
  # Combine the data frames and join other data
  df. <- best.o.m %>%
    left_join(best.tau.m) %>%
    dplyr::mutate(
      ordering.label=factor('ordering', c('ordering', 'pseudotime')),
      pseudotime.label=factor('pseudotime', c('ordering', 'pseudotime'))) %>%
    dplyr::left_join(dl$cell.map)
  ggplot2::ggplot(
      df.,
      aes(x=ordering.label, xend=pseudotime.label,
          y=ordering, yend=pseudotime,
          color=capture)) +
    ggplot2::geom_segment(alpha=.3) + facet_wrap(~ fit)
}
