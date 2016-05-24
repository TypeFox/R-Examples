#' Plot marginal effects from two-way interactions in linear regressions
#'
#' @param obj fitted model object from \code{lm}.
#' @param term1 character string of the first constitutive term of the
#' interaction's variable name.
#' @param term2 character string of the second constitutive term of the
#' interaction's variable name.
#' @param fitted2 numeric vector of fitted values of \code{term2} to plot for.
#' If unspecified, then all unique observed values are used.
#' @param ci numeric. confidence interval level, expressed on the ]0, 100[]
#' interval.
#' @param plot boolean. return plot if TRUE; return data.frame of marginal
#' effects estimates if FALSE
#'
#' @return a \code{gg} class ggplot2 object
#'
#' @examples
#' # Estimate model
#' states <- as.data.frame(state.x77)
#' m1 <- lm(Murder ~ Income * Population, data = states)
#'
#' # Plot marginal effect of Income across the observed range of Population
#' # on the Murder rate
#' plot_me(m1, 'Income', 'Population', ci = 95)
#'
#' # Return marginal effects as a data frame
#' plot_me(m1, 'Income', 'Population', plot = FALSE)
#'
#'
#'
#' @source Inspired by:
#' \url{http://www.statsblogs.com/2013/08/27/creating-marginal-effect-plots-for-linear-regression-models-in-r/}
#'
#' @importFrom stats coef qnorm vcov
#' @import ggplot2
#'
#' @export

plot_me <- function(obj, term1, term2, fitted2, ci = 90, plot = TRUE) {

    if (class(obj) != 'lm') stop('Only lm model objects can be used.',
            call. = FALSE)

    beta_hat <- coef(obj)
    cov <- vcov(obj)

    int_term12 <- sprintf('%s:%s', term1, term2)
    int_term21 <- sprintf('%s:%s', term2, term1)

    if (!(term1 %in% names(beta_hat))) stop(sprintf('%s not found.', term1),
                                            call. = FALSE)
    if (all(!(c(int_term12, int_term21) %in% names(beta_hat)))) {
        stop('Interaction term not found.', call. = FALSE)
    }
    if (int_term12 %in% names(beta_hat)) int_term <- int_term12
    else int_term <- int_term21

    term2_dist <- obj$model[, term2]
    term2_dist <- data.frame(fake_y = 0, real_x = term2_dist)
    names(term2_dist) <- c('term1', 'term2')

    if (missing(fitted2)) fitted2 <- sort(unique(term2_dist$term2))

    # Estimated marginal effect
    dy_dx <- beta_hat[term1] + beta_hat[int_term] * fitted2

    # Standard error
    se_dy_dx <- sqrt(cov[term1, term1] + fitted2 ^ 2 * cov[int_term, int_term] +
                         2 * fitted2 *cov[term1, int_term])

    # Confidence interval
    z <- qnorm(ci / 100)
    upper <- dy_dx + z * se_dy_dx
    lower <- dy_dx - z * se_dy_dx

    parts <- data.frame(cbind(fitted2, dy_dx, lower, upper))

    if (plot) {
        if (length(parts$fitted2) > 5) {
            ggplot(parts, aes(fitted2, dy_dx)) +
                geom_rug(data = term2_dist, aes(x = term2, y = term1), sides = 'b',
                         alpha = 0.1) +
                geom_hline(yintercept = 0, linetype = 'dotted') +
                geom_line() +
                geom_ribbon(ymin = lower, ymax = upper, alpha = 0.1) +
                xlab(sprintf('\n%s', term2)) +
                ylab(sprintf('Marginal effect of\n%s\n', term1)) +
                theme_bw()
        }
        else if (length(parts$fitted2) <= 5) {
            ggplot(parts, aes(fitted2, dy_dx, ymin = lower, ymax = upper)) +
              #  geom_rug(data = term2_dist, aes(x = term2, y = term1), sides = 'b',
              #           alpha = 0.1) +
                geom_hline(yintercept = 0, linetype = 'dotted') +
                geom_pointrange() +
                scale_x_continuous(breaks = parts$fitted) +
                xlab(sprintf('\n%s', term2)) +
                ylab(sprintf('Marginal effect of\n%s\n', term1)) +
                theme_bw()
        }

    } else {
        return(parts)
    }
}
