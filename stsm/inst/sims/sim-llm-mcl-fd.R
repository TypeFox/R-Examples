#
# Simulation experiments:
# the series in the data set 'llm' are fit by 
# maximizing the spectral log-likelihood function (ML-FD)
# where one of the parameters is concentrated out of
# the likelihood function
#
# the following procedures are considered:
# procedure 0: StructTS(), benchmark results
# procedure 1: ML-FD Netwon-Raphson 
# procedure 2: ML-FD scoring algorithm
# procedure 3: ML-FD BFGS algorithm from 'optim()'
# procedure 4: ML-FD L-BFGS-B algorithm from 'optim()'
#

library("stsm")

# load data

##NOTE
# requires data set generated in file "datagen-llm.R" in the 
# same folder as this file

iter <- ncol(llm)

# initial parameter values

initpars <- c(var2 = 1)
initcpar <- c(var1 = 1)
nop <- NULL

# scoring algorithm parameters

tol <- 0.001
maxiter <- 100 #250
step <- NULL # (ignored for 'maxlik.fd.optim()')

# line search parameters 
# they are ignored if 'step' is fixed to a numeric, e.g. 'step <- 1'

ls <- list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1)
#ls <- list(type = "brent.fmin", tol = .Machine$double.eps^0.25, cap = 1)
#ls <- list(type = "wolfe", ftol = 0.0001, gtol = 0.5, cap = 1)

# parameterization

transP <- NULL #NULL #"StructTS" #"square"

# barrier term in log-likelihood function (no barrier with mu = 0)

bar <- list(type = "1", mu = 0)

# storage list

Mres <- vector("list", 6)
names(Mres) <- paste("proc", 0:(length(Mres) - 1), sep = "")
tmp1 <- array(dim = c(iter, length(initpars) + 1), 
  dimnames = list(NULL, names(c(initcpar, initpars))))
tmp2 <- rep(NA, iter)
tmp3 <- array(dim = c(maxiter + 1, length(initpars) + 1, iter))
tmp4 <- array(dim = c(iter, 2), dimnames = list(NULL, c("fcn", "grd")))
Mres <- lapply(Mres, function(x) list(M = tmp1, iter = tmp2,
  time = tmp2, paths = tmp3, lscounts = tmp4))

rm(tmp1, tmp2, tmp3, tmp4)

# begin main loop

for (i in seq_len(iter))
{
  m <- stsm.model(model = "local-level", y = llm[,i], 
    pars = initpars, cpar = initcpar, nopars = nop, transPars = transP)

  # procedure 0: ML-TD StructTS(), benchmark results

  res0 <- StructTS(m@y, type = "level")$coef[c(2,1)]

  Mres$proc0$M[i,] <- res0

  # procedure 1: ML-FD Netwon-Raphson 
  # with optimized step length if 'step' is NULL

  res1 <- maxclik.fd.scoring(m = m, step = NULL, 
    information = "observed", ls = ls, barrier = bar, 
    control = list(maxit = maxiter, tol = tol, trace = TRUE))

  Mres$proc1$M[i,] <- c(res1$model@cpar, res1$model@cpar * res1$par)
  Mres$proc1$iter[i] <- res1$iter
  Mres$proc1$paths[seq(res1$iter+1),,i] <- res1$Mpars
  if (length(res1$ls.counts) == 2) # 'optimize()' does not return this
    Mres$proc1$lscounts[i,] <- res1$ls.counts

#i=1
#> res1$par
#      var2 
#0.02742823 
#> res1$iter
#[1] 7
#> get.cpar(res1$model, rescale = TRUE)
#    var1 
#1722.995 
#> get.pars(res1$model, rescale = TRUE)
#    var2 
#47.25869 
#> c(res1$model@cpar, res1$model@cpar * res1$par)
#      var1       var1 
#1722.99473   47.25869 
#> coef(res1)
#      var1       var2 
#1722.99473   47.25869 

  # procedure 2: ML-FD scoring algorithm
  # with optimized step length if 'step' is NULL

  res2 <- maxclik.fd.scoring(m = m, step = NULL, 
    information = "expected", ls = ls, barrier = bar, 
    control = list(maxit = maxiter, tol = tol, trace = TRUE))

  Mres$proc2$M[i,] <- c(res2$model@cpar, res2$model@cpar * res2$par)
  Mres$proc2$iter[i] <- res2$iter
  Mres$proc2$paths[seq(res2$iter+1),,i] <- res2$Mpars
  if (length(res2$ls.counts) == 2)
    Mres$proc2$lscounts[i,] <- res2$ls.counts

#i=1
#> res2$par
#      var2 
#0.02742822 
#> res2$iter
#[1] 2
#> get.cpar(res2$model, rescale = TRUE)
#    var1 
#1722.995 
#> get.pars(res2$model, rescale = TRUE)
#    var2 
#47.25868 
#> c(res2$model@cpar, res2$model@cpar * res2$par)
#      var1       var1 
#1722.99476   47.25868 

  # procedure 3: ML-FD BFGS algorithm from 'optim()'
  # NOTE 'count' is recorded (not 'iter')

  res3 <- try(maxlik.fd.optim(m, inf = 99999, method = "BFGS", 
    gr = "analytical"), silent = TRUE)
  #try(res3 <- maxlik.fd.optim(m, inf = 99999, method = "BFGS", 
  #  gr = "numerical", hessian = FALSE), silent = TRUE)

  if (!inherits(res3, "try-error"))
  {
    Mres$proc3$M[i,] <- c(res3$model@cpar, res3$model@cpar * res3$par)
    Mres$proc3$iter[i] <- res3$iter[1]
  }

#i=1
#> res3$par
#      var2 
#0.02742814 
#> res3$iter
#function gradient 
#      45        9 
#> get.cpar(res3$model, rescale = TRUE)
#    var1 
#1722.995 
#>  get.pars(res3$model, rescale = TRUE)
#    var2 
#47.25856 
#> c(res3$model@cpar, res3$model@cpar * res3$par)
#      var1       var1 
#1722.99516   47.25856 

  # procedure 4: ML-FD L-BFGS-B algorithm from 'optim()'

  #try(res4 <- maxlik.fd.optim(m, barrier = bar, inf = 99999, 
  #  method = "L-BFGS-B", gr = "analytical", hessian = FALSE), silent = TRUE)
  res4 <- try(maxlik.fd.optim(m, barrier = bar, inf = 99999, 
    method = "L-BFGS-B", gr = "numerical"), silent = TRUE)

  # record the path followed by the optimization method
  # using the same stopping criterion as in 'maxlik.fd.scoring()'
  if (!inherits(res4, "try-error"))
  {
    Mpars <- c(get.cpar(m, TRUE), get.pars(m, TRUE))
    conv <- FALSE
    r <- 0
    while (!conv)
    {
      optout <- maxlik.fd.optim(m = m, barrier = bar, inf = 99999, 
        #method = "L-BFGS-B", gr = "analytical", hessian = FALSE,
        method = "L-BFGS-B", gr = "numerical", hessian = FALSE,
        optim.control = list(trace = FALSE, maxit = r))
      #Mpars <- rbind(Mpars, optout$par)
      Mpars <- rbind(Mpars, 
        c(get.cpar(optout$model, TRUE), get.pars(optout$model, TRUE)))
      #conv <- optout$convergence == 0
      if (sqrt(sum((Mpars[r+1,-1] - Mpars[r+2,-1])^2)) < tol || r > maxiter)
        conv <- TRUE
      r <- r + 1
    }

    Mres$proc4$M[i,] <- Mpars[r+1,]
    Mres$proc4$iter[i] <- r
    Mres$proc4$paths[seq(r+1),,i] <- Mpars
  }

#i=1
#> res4$par
#      var2 
#0.02744035 
#> r
#[1] 6
#> get.cpar(res4$model, rescale = TRUE)
#    var1 
#1722.932 
#> get.pars(res4$model, rescale = TRUE)
#    var2 
#47.27785 
#> c(res4$model@cpar, res4$model@cpar * res4$par)
#      var1       var1 
#1722.93189   47.27785 

  # trace simulation

  #print(i)
  trace.iter <- 100 * (i / iter)
  if (trace.iter %% 10 == 0)
  {
    cat(paste(trace.iter, "% complete.\n", sep = ""))
    #print(colMeans(Mres$proc1$M, na.rm = TRUE))
    #print(colMeans(Mres$proc2$M, na.rm = TRUE))
    #print(colMeans(Mres$proc3$M, na.rm = TRUE))
    #print(colMeans(Mres$proc4$M, na.rm = TRUE))
  }
}

# store results

#setwd("")
#save(Mres, file = "sim-llm-mcl-fd.rda")

# summary of results

colMeans(Mres$proc0$M)
colMeans(Mres$proc1$M)
colMeans(Mres$proc2$M)
colMeans(Mres$proc3$M)
colMeans(Mres$proc4$M)

> colMeans(Mres$proc2$M)
    var1     var2 
1597.116  123.611 
> colMeans(Mres$proc4$M)
     var1      var2 
1596.9995  123.6166 

summary(Mres$proc1$M)
summary(Mres$proc2$M)

mean(Mres$proc0$iter)
mean(Mres$proc1$iter)
mean(Mres$proc2$iter)
mean(Mres$proc3$iter)
mean(Mres$proc4$iter)

> mean(Mres$proc2$iter)
[1] 2.068
> mean(Mres$proc4$iter)
[1] 6.036

summary(Mres$proc1$iter)
summary(Mres$proc2$iter)
summary(Mres$proc4$iter)
