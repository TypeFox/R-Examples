#' Performs the bootstrap for one condition
#' \code{one_bootstrap} performs the bootstrap for one condition
#' @import boot DEoptim
#' @importFrom tidyr gather
#' @keywords internal
#' @export
parn <- 'No te quejes'
one_bootstrap <- function(d, x, k, n, psyfunguesslapses, funname,
                           guess, lapses, parini, pariniset, optimization,
                          bootstrap, B,
                          groups, ypred) {

  if (length(groups) != 0) ypred <- semi_join(ypred, d, by = groups)

  if (bootstrap == 'parametric') ypred <- ypred$ypred
  if (bootstrap == 'nonparametric') ypred <- d[[k]] / d[[n]]

  calculate_par <- function(f)
    parameters(f, x, k, n, psyfunguesslapses, funname,
               parini, pariniset, guess, lapses, optimization, groups)$par

  create_fake_data <- function(f, mle){
    kfake <- rbinom(length(f[[x]]), f[[n]], mle)
    f[[k]] <- kfake
    f$y <- kfake / f[[n]]
    f
  }

  b <- boot(d, calculate_par, R = B, sim = 'parametric',
            ran.gen = create_fake_data, mle = ypred)
  fake_par <- b$t
  colnames(fake_par) <- paste0('p',1:length(fake_par[1,]))
  long <- data.frame(fake_par, sample = 1:length(fake_par[,1]))
  long %>% gather(parn, par, -sample) %>% arrange(sample)

}

