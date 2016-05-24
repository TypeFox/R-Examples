`MCepi` <-
  function (init, params, vac, costs, T = 40, MCreps = 1000,
            quant = c(0.025, 0.975), midepi = FALSE, start = 7, ...) 
{
  SEQ <- SimulateEpidemicQuantiles
  out <- SEQ(init, T, params, vac$frac, vac$stop, costs, MCreps,
             quant[1], quant[2], midepi, start, ...)
  out$quant <- quant
  out$call <- match.call()
  class(out) <- "MCepi"
  return(out)
}

`MCmanage` <-
  function (init, epistep, vacgrid, costs, pinit = list(b = 0.1, k = 0.02, nu = 0.2, mu = 0.1),
            hyper = list(bh = c(1, 3), kh = c(1, 3), nuh = c(1, 1), muh = c(1, 1)),
            vac0 = list(frac=0, stop=0), T = 40, MCreps = 30,
            MCvits = 50, MCMCpits = 1000, bkrate=1, vacsamps = 50,
            quant = c(0.025, 0.975), start = 7, ...) 
{
  out <- SimulateManagementQuantiles(epistep = epistep, T, 
                                     pinit = pinit, init = init, hyper = hyper,
                                     vac0=vac0, costs = costs, vacgrid = vacgrid, MCvits = MCvits, 
                                     MCMCpits = MCMCpits, bkrate=bkrate, vacsamps = vacsamps,
                                     start = start + 1, nreps = MCreps, lowerq = quant[1],
                                     upperq = quant[2], ...)
  out$quant <- quant
  out$call <- match.call()
  class(out) <- "MCepi"
  return(out)
}

`optvac` <-
  function (init, params, vacgrid, costs, T = 40, MCvits = 100, 
            midepi = FALSE, start = 7) 
{
  VSTP <- VarStopTimePolicy
  C <- VSTP(init$S0, init$I0, T, params$b, params$k, params$nu, 
            params$mu, costs$vac, costs$death, costs$infect, MCvits, 
            vacgrid$fracs, vacgrid$stops, midepi, start)
  out <- list(C=C, vacgrid=vacgrid)
  out$call <- match.call()
  class(out) <- "optvac"
  return(out)
}

`getcost` <-
  function (obj) 
{
  if(class(obj) == "epiman") {
    T <- nrow(obj$soln)
    return(obj$soln$C[T])
  } else {
    T <- length(obj$Median$C)
    q1 <- paste("q", signif(obj$quant[1], 3), sep="")
    q3 <- paste("q", signif(obj$quant[2], 3), sep="")
    df <- data.frame(obj$Q1$C[T], obj$Mean$C[T], obj$Median$C[T], obj$Q3$C[T])
    names(df) <- c(q1, "mean", "median", q3)
    return(df)
  }
}


`getvac` <-
  function (obj) 
{
  if(class(obj) == "epiman") {
    T <- nrow(obj$soln)
    return(obj$soln$V[T])
  } else {
    T <- length(obj$Median$C)
    q1 <- paste("q", signif(obj$quant[1], 3), sep="")
    q3 <- paste("q", signif(obj$quant[2], 3), sep="")
    df <- data.frame(obj$Q1$V[T], obj$Mean$V[T], obj$Median$V[T], obj$Q3$V[T])
    names(df) <- c(q1, "mean", "median", q3)
    return(df)
  }
}


`getpolicy` <-
  function (obj, which = c("best", "worst")) 
{
  which <- match.arg(which)
  if (which == "best") {
    pol <- which(obj$C == min(obj$C), arr.ind = TRUE)[1,]
    frac <- obj$vacgrid$fracs[pol[1]]
    stop <- obj$vacgrid$stops[pol[2]]
    name <- "best"
  }
  else {
    pol <- which(obj$C == max(obj$C), arr.ind = TRUE)[1,]
    frac <- obj$vacgrid$fracs[pol[1]]
    stop <- obj$vacgrid$stops[pol[2]]
    name <- "worst"
  }
  cost <- ceiling(obj$C[pol[1], pol[2]])
  df <- data.frame(row = pol[1], col = pol[2], frac = frac, 
                   stop = stop, cost = cost)
  row.names(df) <- name
  return(df)
}
