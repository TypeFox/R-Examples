# SIMULATIONS - droplet, array ------------------------------
#summary workhorse + plot, should be not called directly by user

#dube's simulation
sim_dpcr_dube <- function(m, n, times, pos_sums, n_panels) {
  mt <- m*times
  nt <- n*times
  distr_molec <- sample.int(nt, mt, replace = TRUE)
  counts <- rep(0, nt)
  for (i in distr_molec)
    counts[i] <- counts[i] + 1
  counts <- sample(counts) #is additional shuffle needed?
  #extract panels and calculate the number of positive partion (if required)
  if (pos_sums == TRUE) {
    times_v <- 1L:n_panels
    gr_id2 <- times_v*n
    gr_id1 <- gr_id2 - n + 1
    matrix(vapply(times_v, function(i) sum(counts[gr_id1[i]:gr_id2[i]] > 0), 0), 
           ncol = n_panels)
  } else {
    times_v <- 1L:n_panels
    gr_id2 <- times_v*n
    gr_id1 <- gr_id2 - n + 1
    sapply(times_v, function(i) counts[gr_id1[i]:gr_id2[i]])
  }
}

#our simple simulation
#assumption: the number of molecules in each panel is determined by the normal distribution
#with mean  equal to m
sim_dpcr_multi <- function(m, n, times, pos_sums, n_panels) {
  total_m <- m*times
  ms <- rnorm(times, m, 0.05*m)
  ms <- round(ms + (total_m - sum(ms))/times, 0)
  delta <- total_m - sum(ms)
  indices <- sample(1:times, abs(delta))
  ms[indices] <- ms[indices] + sign(delta)
  #extract panels and calculate the number of positive partion (if required)
  if (pos_sums) {
    matrix(vapply(ms[1L:n_panels], function(x) 
      sum(table(factor(sample.int(n, x, replace = TRUE), levels = 1L:n)) > 0), 0), 
      ncol = n_panels)
  } else {
    vapply(ms[1L:n_panels], function(x) 
      as.vector(table(factor(sample.int(n, x, replace = TRUE), levels = 1L:n))), 
      rep(0, n))
  }
}


sim_dpcr <- function(m, n, times, dube, pos_sums, n_panels) {
  n <- num2int(n)
  if (n_panels > times) 
    stop("The 'n_panels' argument cannot have larger value than the 'times'
         argument", call. = TRUE, domain = NA)
  if (m < 0) 
    stop("'m' must be a non-negative integer.", call. = TRUE, domain = NA)
  if (n < 1)  
    stop("'n' must be an integer bigger than 1.", call. = TRUE, domain = NA)
  if (dube) {
    res <- sim_dpcr_dube(m, n, times, pos_sums, n_panels)
  } else {
    res <- sim_dpcr_multi(m, n, times, pos_sums, n_panels)
  }
  res
}