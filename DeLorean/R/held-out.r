#' Calculate posterior covariance and estimate parameters for
#' held out genes given pseudotimes estimated by DeLorean model.
#'
#' @param dl de.lorean object
#' @param held.out Held out gene expression levels
#' @param posterior.sample Posterior sample to use
#'
#' @export
#'
held.out.posterior <- function(
    dl,
    held.out,
    posterior.sample=dl$best.sample)
{
    # Zero mean the genes
    gene.means <- (
        held.out
        %>% group_by(gene)
        %>% dplyr::summarise(mean=mean(x))
    )
    held.out <- (
        held.out
        %>% left_join(gene.means)
        %>% mutate(x.adj = x - mean)
    )
    with(dl, {
        tau <- (
            samples.l$tau
            %>% filter(iter == posterior.sample)
            %>% dplyr::select(cell, tau)
        )
        held.out <- (
            held.out
            %>% left_join(tau)
        )
        cov.fn <- function(K.tau, psi, omega) {
            psi * K.tau + omega * diag(nrow(K.tau))
        }
        calc.K <- function(.data, psi, omega) {
            K.tau <- cov.matern.32(
                cov.calc.dists(.data$tau, period=opts$period),
                opts$length.scale)
            # stopifnot(is.positive.definite(K.tau))
            K <- cov.fn(K.tau, psi, omega)
            # stopifnot(is.positive.definite(K))
            K
        }
        optimise <- function(.data) {
            likelihood <- function(par) {
                psi <- exp(par[1])
                omega <- exp(par[2])
                K <- calc.K(.data, psi, omega)
                gp.ll <- gp.log.marg.like(.data$x.adj, K=K)
                return(
                    gp.ll[1,1]
                    + dlnorm(psi,
                            meanlog=hyper$mu_psi,
                            sdlog=sqrt(hyper$sigma_psi),
                            log=TRUE)
                    + dlnorm(omega,
                            meanlog=hyper$mu_omega,
                            sdlog=sqrt(hyper$sigma_omega),
                            log=TRUE)
                )
            }
            # Find the optimal values of psi and omega
            optimum <- optim(
                c(hyper$mu_psi, hyper$mu_omega),
                likelihood,
                control=list(fnscale=-1))
            psi <- exp(optimum$par[1])
            omega <- exp(optimum$par[2])
            return(data.frame(psi=psi, omega=omega))
        }
        calc.posterior <- function(params) {
            .data <- held.out %>% filter(gene == params$gene)
            # Make posterior predictions on the test inputs
            K <- calc.K(.data, params$psi, params$omega)
            Kstar <- params$psi * cov.matern.32(
                cov.calc.dists(.data$tau, test.input, period=opts$period),
                opts$length.scale)
            Kstarstar <- params$psi * cov.matern.32(
                cov.calc.dists(test.input, period=opts$period),
                opts$length.scale)
            posterior <- gp.predict(.data$x.adj, K, Kstar, Kstarstar)
            posterior.df <- (
                gp.predictions.df(posterior)
                %>% dplyr::mutate(gene=params$gene)
                %>% left_join(gene.means)
                %>% dplyr::mutate(mean=mu + mean)
                %>% dplyr::rename(var=Sigma)
            )
            posterior.df$x <- test.input
            posterior.df
        }
        params <- (
            held.out
            %>% group_by(gene)
            %>% do(optimise(.))
            %>% ungroup()
        )
        posterior <- (
            params
            %>% group_by(gene)
            %>% do(calc.posterior(.))
            %>% ungroup()
        )
        return(list(held.out=held.out, params=params, posterior=posterior))
    })
}


#' Select held out genes by those with highest variance
#'
#' @param dl de.lorean object
#' @param expr Expression matrix of all genes
#' @param num.held.out Number to select
#'
#' @export
#'
held.out.select.genes <- function(dl, expr, num.held.out) {
    cells.fit <- as.character(dl$cell.map$cell)
    genes.fit <- as.character(dl$gene.map$gene)
    held.out.var <- apply(
        expr[! rownames(expr) %in% genes.fit, cells.fit],
        1,
        var)
    names(tail(sort(held.out.var), num.held.out))
}


#' Melt held out genes
#'
#' @param dl de.lorean object
#' @param expr Expression matrix of all genes
#' @param held.out.genes Genes to hold out
#'
#' @export
#'
held.out.melt <- function(dl, expr, held.out.genes) {
    cells.fit <- as.character(dl$cell.map$cell)
    melt.expr(dl, expr[held.out.genes, cells.fit])
}


#' Filter the genes
#'
#' @param posterior The posterior of some held out genes
#' @param genes Genes to filter
#'
held.out.posterior.filter <- function(posterior, genes) {
  return(list(
    held.out = posterior$held.out %>% filter(gene %in% genes),
    params = posterior$params %>% filter(gene %in% genes),
    posterior = posterior$posterior %>% filter(gene %in% genes)))
}


#' Order the genes by the variation of their posterior mean
#'
#' @param posterior The posterior of some held out genes
#'
held.out.posterior.by.variation <- function(posterior) {
  return((
    posterior$posterior %>%
    group_by(gene) %>%
    dplyr::summarise(dynamic.var=var(mu)) %>%
    arrange(-dynamic.var))$gene)
}


#' Join with another data frame. Useful for adding gene names etc..
#'
#' @param posterior The posterior of some held out genes
#' @param .df Data frame to join with
#'
held.out.posterior.join <- function(posterior, .df) {
  return(list(
      held.out = posterior$held.out %>% left_join(.df),
      params = posterior$params %>% left_join(.df),
      posterior = posterior$posterior %>% left_join(.df)))
}



#' Plot the posterior of held out genes
#'
#' @param dl de.lorean object
#' @param posterior The posterior of some held out genes
#' @param facets Variables to wrap facets on
#'
plot.held.out.posterior <- function(dl, posterior, facets=~ gene) (
    plot.add.mean.and.variance(
        ggplot(posterior$posterior %>% left_join(dl$gene.meta))) +
    geom_point(
        data=posterior$held.out %>% left_join(dl$cell.meta),
        aes(x=tau, y=x, color=capture)) +
    facet_wrap(facets)
)
