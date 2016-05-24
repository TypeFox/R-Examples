#' Pair comparisons of the parameters using bootstrap (for each parameter)
#' \code{one_parcomparisons} Calculates the bootstrap confidence intervals for the
#' difference in the parameters for two groups for all possible pairs
#' of groups (for each parameter)
#' @keywords internal
#' @export
one_parcomparisons <- function(d, para, groups, ci) {
  V1 <- NULL
  V2 <- NULL
  diffe <- NULL
  difinf <- NULL
  difsup <- NULL
  condlevels <- d %>% filter(sample==1)
  if (length(groups) != 0) para <- semi_join(para, condlevels, by ='parn')
  condlevels <- condlevels %>% select(match(groups,names(condlevels)))

  combinations <- as.data.frame(t(combn(nrow(condlevels),2)))

  putnames <- function(f) {
    cond1 <- condlevels[f$V1,]
    cond2 <- condlevels[f$V2,]
    names(cond2) <- paste0(names(cond2), '2')
    cbind(cond1, cond2)
  }
  pairnames <- combinations %>% group_by(V1,V2) %>%
    do(putnames(.))

  one_difcom <- function(f) {
    cond1 <- semi_join(para, condlevels[f$V1,], by=groups)
    cond2 <- semi_join(para, condlevels[f$V2,], by=groups)
    data.frame(dif= cond1$par -cond2$par)
  }
  dif <- pairnames %>% group_by_(.dots = names(pairnames)) %>%
    do(one_difcom(.)) %>% ungroup() %>% select(-V1, -V2)


  one_bootcom <- function(f) {
    cond1 <- semi_join(d, condlevels[f$V1,], by=groups)
    cond2 <- semi_join(d, condlevels[f$V2,], by=groups)

    data.frame(par1 = cond1$par, par2= cond2$par,
               diffe = cond1$par -cond2$par) %>%
      summarise(difinf = quantile(diffe, .5*(1 - ci))[[1]],
                difsup = quantile(diffe, 1 - .5*(1 - ci))[[1]],
                signif = ifelse(difinf * difsup < 0,' ','*'))
  }
  ci <- pairnames %>% group_by_(.dots = names(pairnames)) %>%
    do(one_bootcom(.)) %>% ungroup() %>% select(-V1, -V2)

 merge(dif,ci)

}

