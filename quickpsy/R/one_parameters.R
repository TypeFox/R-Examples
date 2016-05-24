#' Obtains the parameters for one condition
#' \code{one_parameters} obtains the parameters for one condition
#' @keywords internal
#' @export
one_parameters <- function(d, x, k, n, psyfunguesslapses, funname, parini,
                           pariniset, guess, lapses, optimization, groups) {
  nllfun <- create_nll(d, x, k, n, psyfunguesslapses)

  if (optimization == 'DE') {
    if (is.data.frame(parini) || is.atomic(parini))
      stop('parini should be specified as a list of the type list(c(par1min, par1max), c(par2min, par2max),...', call. = F)
    else if (is.list(parini)) {
      parini <- matrix(unlist(parini), ncol = 2, byrow = T)
      mod <- DEoptim(nllfun, lower = parini[,1], upper = parini[,2])$optim
      para <- mod$bestmem
    }
    else
      stop('parini should be specified as a list of the type list(c(par1min, par1max), c(par2min, par2max),...', call. = F)

  }
  if (optimization== 'optim') {

    if (pariniset) {
      if (is.atomic(parini))
        para <- optim(parini, nllfun)$par
      if (is.list(parini)){
        parini <- matrix(unlist(parini), ncol = 2, byrow = T)
        para <- optim(.5 * (parini[,1] + parini[,2]),
                      nllfun, method = 'L-BFGS-B',
                      lower = parini[,1],
                      upper = parini[,2])$par
      }
    }
    else {
      if (length(groups) == 0) parini <- parini$par
      else parini <- semi_join(parini, d, by = groups)$par

      if (funname == 'weibull_fun') {
        if (parini[1] < 0) parini[1] <- .Machine$double.eps
        if (parini[2] < 0) parini[2] <- .Machine$double.eps
      }
      if (funname == 'cum_normal_fun') {
        if (parini[2] < 0) parini[2] <- 0
      }

      para <- optim(parini, nllfun)$par
    }
  }
  if (optimization != 'DE' && optimization != 'optim')
    stop('optimization should be \'optim \' or \'DE\'.', call. = F)
  data.frame(parn = paste0('p', seq(1, length(para))), par = para)
}


