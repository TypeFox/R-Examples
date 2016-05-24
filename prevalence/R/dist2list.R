## BUGS-dist to list

dist2list <-
function(d, type) {
  ## check if d is numeric value
  if (!is.na(suppressWarnings(as.numeric(d)))) {
    x <- list(dist = "fixed", par = as.numeric(d))
    out <- checkSeSp(x, type)

  } else {
    ## extract distribution and parameters
    dst <- as.character(parse(text = d)[[1]])[1]
    par <- as.numeric(as.character(parse(text = d)[[1]])[-1])

    ## check distribution
    if (!any(c("fixed", "dunif", "dbeta", "dpert") == dst))
      stop(paste("Distribution must be either",
                 "'fixed', 'dunif', 'dbeta' or 'dpert'"))

    ## check length of par
    dst_nr <- which(c("fixed", "dunif", "dbeta", "dpert") == dst)
    len <- c(1, 2, 2, 3)[dst_nr]
    if (length(par) != len)
      stop(paste("Distribution", dst, "requires", len, "parameters"))

    ## check dist and pars
    x <-
      switch(
        dst[1],
        "fixed" = list(dist = "fixed",
                       par = par[1]),
        "dunif" = list(dist = "uniform",
                       min = par[1], max = par[2]),
        "dbeta" = list(dist = "beta",
                       alpha = par[1], beta = par[2]),
        "dpert" = list(dist = "pert",
                       a = par[1], m = par[2], b = par[3]))
    out <- checkSeSp(x, type)
  } 

  ## return distribution in list format
  return(out)
}