#' @title  Fit the PRO model
#'
#' @description
#' Fit the PRO model to data. 
#' Reference: X Luo, S Gee, V Sohal, D Small (In Press). A Point-process Response Model for Optogenetics Experiments on Neural Circuits. _Statistics in Medicine_.
#' @export
#'
#' @importFrom stats as.formula binomial glm
#' 
#' @param spike A binary vector represents spiking (1) or no spiking (0).
#' @param flash A binary vector of the same length of \code{spike}, 1 for flashing and 0 for non-flashing.
#' @param ... Additional parameters, see \code{\link{model.pro}}.
#' @return  a \code{\link{glm}} object of the fitted PRO coefficients.
#' @examples
#' n <- 500
#' set.seed(100)
#' re <- sim.lif(n, rbinom(n, 1, 0.14), 7, 3)
#' fit.pro <- pro(re$sbin, re$I)
#' summary(fit.pro)
pro <- function(spike, flash, ...) {
    model.f <- "spike~logpf+logcumf+loglogcumsqpf+logcumf:loglogcumsqpf"
    dd <-  model.pro(spike, flash, ...)
    fit <- glm(as.formula(model.f), family=binomial, data=dd)
    return(fit) 
}



