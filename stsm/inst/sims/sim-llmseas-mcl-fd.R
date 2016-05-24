#
# Simulation experiments:
# the series in the data set 'llmseas' are fit by 
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
# requires data set generated in file "datagen-llmseas.R" in the 
# same folder as this file

iter <- ncol(llmseas)

# initial parameter values

initpars <- c(var2 = 1, var3 = 1)
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

# graphical parameters to trace the simulation

plotid <- c(1, 2, 4) + 1
#plab <- names(Mres)[plotid]
plab <- c("StructTS", "Newton-Raphson", "Scoring", "BFGS", "L-BFGS-B")[plotid]
colors <- c("red", "green", "blue")
arrow.angle <- c(20, 30, 40)
legtext <- NULL

# begin main loop

for (i in seq_len(iter))
{
  m <- stsm.model(model = "llm+seas", y = llmseas[,i], 
    pars = initpars, cpar = initcpar, nopars = nop, transPars = transP)

  # procedure 0: ML-TD StructTS(), benchmark results

  res0 <- StructTS(m@y, type = "BSM", 
    fixed = c(NA, 0, NA, NA))$coef[c(4,1,3)]

  Mres$proc0$M[i,] <- res0

#i=1
#  epsilon     level      seas 
#172.29956  35.51876 176.74453 

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
#      var2       var3 
#0.05535256 0.82006341 
#> res1$iter
#[1] 8
#> get.cpar(res1$model, rescale = TRUE)
#    var1 
#200.3328 
#> get.pars(res1$model, rescale = TRUE)
#     var2      var3 
# 11.08893 164.28562 
#> c(res1$model@cpar, res1$model@cpar * res1$par)
#     var1      var2      var3 
#200.33282  11.08893 164.28562 
#> coef(res1)
#     var1      var2      var3 
#200.33282  11.08893 164.28562 

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
#      var2       var3 
#0.05535222 0.82002091 
#> res2$iter
#[1] 6
#> get.cpar(res2$model, rescale = TRUE)
#    var1 
#200.3375 
#> get.pars(res2$model, rescale = TRUE)
#     var2      var3 
# 11.08912 164.28091 
#> c(res2$model@cpar, res2$model@cpar * res2$par)
#     var1      var2      var3 
#200.33746  11.08912 164.28091 

  # procedure 3: ML-FD BFGS algorithm from 'optim()'
  # NOTE 'count' is recorded (not 'iter')

  res3 <- try(maxlik.fd.optim(m, inf = 99999, method = "BFGS", 
    gr = "analytical"), silent = TRUE)
  #res3 <- try(maxlik.fd.optim(m, inf = 99999, method = "BFGS", 
  #  gr = "numerical"), silent = TRUE)

  if (!inherits(res3, "try-error"))
  {
    Mres$proc3$M[i,] <- c(res3$model@cpar, res3$model@cpar * res3$par)
    Mres$proc3$iter[i] <- res3$iter[1]
  }

#i=1
#> res3$par
#      var2       var3 
#0.05535264 0.82006470 
#> res3$iter
#function gradient 
#      52       13 
#> get.cpar(res3$model, rescale = TRUE)
#    var1 
#200.3327 
#> get.pars(res3$model, rescale = TRUE)
#     var2      var3 
# 11.08894 164.28574 
#> c(res3$model@cpar, res3$model@cpar * res3$par)
#     var1      var2      var3 
#200.33265  11.08894 164.28574 

  # procedure 4: ML-FD L-BFGS-B algorithm from 'optim()'

  #res4 <- try(maxlik.fd.optim(m, barrier = bar, inf = 99999, 
  #  method = "L-BFGS-B", gr = "analytical"), silent = TRUE)
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
        method = "L-BFGS-B", gr = "numerical", 
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
#      var2       var3 
#0.05538285 0.82012291 
#> r
#[1] 13
#> get.cpar(res4$model, rescale = TRUE)
#   var1 
#200.315 
#> get.pars(res4$model, rescale = TRUE)
#     var2      var3 
# 11.09402 164.28293 
#> c(res4$model@cpar, res4$model@cpar * res4$par)
#     var1      var2      var3 
#200.31502  11.09402 164.28293 

  # trace simulation

  # '[,-1,i]' the first column is removed, 
  # it is the value inferred from the maximized likelihood
  # the parameter concentrated out of the likelihood function
  tmp <- rbind(na.omit(Mres$proc1$paths[,-1,i]), na.omit(Mres$proc2$paths[,-1,i]),
    na.omit(Mres$proc3$paths[,-1,i]), na.omit(Mres$proc4$paths[,-1,i]))
  rg <- apply(tmp, 2, range)

  options(warn = -1)
  plot(rg[,1], rg[,2], type = "n", xlab = "var2", ylab = "var3", 
    main = paste("series", i))

  for (k in seq_along(plotid))
  {
    idk <- plotid[k]
    mp <- na.omit(Mres[[idk]]$paths[,,i])
    if (length(mp) > 0) 
    {
      for (j in seq(2, nrow(mp)))
      {
        arrows(mp[j-1,2], mp[j-1,3], mp[j,2], mp[j,3], 
          lty = 1, col = colors[k], length = 0.15, angle = arrow.angle[k])
      }
    }

    legtext <- c(legtext, 
      paste(plab[k], ": ", Mres[[idk]]$iter[i], " iter.", sep = ""))
  }

  legend("topright", legtext, cex = 1, adj = 0,
    col = colors, lty = 1, xjust = 0, yjust = 2.0, bty = "n", horiz = FALSE)
  legtext <- NULL
  options(warn = 0)

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
#save(Mres, file = "sim-llmseas-mcl-fd.rda")

# summary of results

colMeans(Mres$proc0$M)
colMeans(Mres$proc1$M)
colMeans(Mres$proc2$M)
colMeans(Mres$proc3$M)
colMeans(Mres$proc4$M)

> colMeans(Mres$proc2$M)
     var1      var2      var3 
291.40347  11.22479 107.83065 
> colMeans(Mres$proc4$M)
     var1      var2      var3 
290.79962  11.22067 108.38712 

summary(Mres$proc1$M)

mean(Mres$proc0$iter)
mean(Mres$proc1$iter)
mean(Mres$proc2$iter)
mean(Mres$proc3$iter)
mean(Mres$proc4$iter)

> mean(Mres$proc2$iter)
[1] 4.56
> mean(Mres$proc4$iter)
[1] 11.832

summary(Mres$proc2$iter)
