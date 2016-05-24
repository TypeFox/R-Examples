#
# Simulation experiments:
# the series in the data set 'llmseas' are fit by 
# maximizing the spectral log-likelihood function (ML-FD)
#
# the following procedures are considered:
# procedure 0: StructTS(), benchmark results
# procedure 1: ML-FD Netwon-Raphson 
# procedure 2: ML-FD scoring algorithm
# procedure 3: ML-FD Newton-Raphson experimental version
# procedure 4: ML-FD BFGS algorithm from 'optim()'
# procedure 5: ML-FD L-BFGS-B algorithm from 'optim()'
#
# some of the results obtained here are reported in Table 1
# in the vignette of the 'stsm' package
#

##NOTE added in optim() convergence if (... || r > maxiter)

library("stsm")

# load data

##NOTE
# requires data set generated in file "datagen-llmseas.R" in the 
# same folder as this file

iter <- ncol(llmseas)

# initial parameter values

initpars <- c(var1 = 1, var2 = 1, var3 = 1)
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
tmp1 <- array(dim = c(iter, length(initpars)), dimnames = list(NULL, names(initpars)))
tmp2 <- rep(NA, iter)
tmp3 <- array(dim = c(maxiter + 1, length(initpars), iter))
tmp4 <- array(dim = c(iter, 2), dimnames = list(NULL, c("fcn", "grd")))
Mres <- lapply(Mres, function(x) list(M = tmp1, iter = tmp2,
  time = tmp2, paths = tmp3, lscounts = tmp4))

rm(tmp1, tmp2, tmp3, tmp4)

# graphical parameters to trace the simulation

plotid <- c(1, 2, 3, 5) + 1
#plab <- names(Mres)[plotid]
plab <- c("StructTS", "Newton-Raphson", "Scoring", 
  "NR-mod.", "BFGS", "L-BFGS-B")[plotid]
colors <- c("red", "green", "blue", "orange")
arrow.angle <- c(20, 30, 40, 50)
legtext <- NULL

# begin main loop

for (i in seq_len(iter))
{
  m <- stsm.model(model = "llm+seas", y = llmseas[,i], 
    pars = initpars, nopars = nop, transPars = transP)

  # procedure 0: ML-TD StructTS(), benchmark results

  res0 <- StructTS(m@y, type = "BSM", 
    fixed = c(NA, 0, NA, NA))$coef[c(4,1,3)]

  Mres$proc0$M[i,] <- res0

#i=1
#  epsilon     level      seas 
#172.29956  35.51876 176.74453 

  # procedure 1: ML-FD Netwon-Raphson 
  # with optimized step length if 'step' is NULL

  res1 <- maxlik.fd.scoring(m = m, step = step, 
    information = "observed", ls = ls, barrier = bar, 
    control = list(maxit = maxiter, tol = tol, trace = TRUE))

  Mres$proc1$M[i,] <- res1$par
  Mres$proc1$iter[i] <- res1$iter
  Mres$proc1$paths[seq(res1$iter+1),,i] <- res1$Mpars
  if (length(res1$ls.counts) == 2) # 'optimize()' does not return this
    Mres$proc1$lscounts[i,] <- res1$ls.counts

#i=1
#> res1$par
#     var1      var2      var3 
#200.33269  11.08894 164.28572 
#> res1$iter
#[1] 18

  # procedure 2: ML-FD scoring algorithm
  # with optimized step length if 'step' is NULL

  res2 <- maxlik.fd.scoring(m = m, step = step, 
    information = "expected", ls = ls, barrier = bar, 
    control = list(maxit = maxiter, tol = tol, trace = TRUE))

  Mres$proc2$M[i,] <- res2$par
  Mres$proc2$iter[i] <- res2$iter
  Mres$proc2$paths[seq(res2$iter+1),,i] <- res2$Mpars
  if (length(res2$ls.counts) == 2)
    Mres$proc2$lscounts[i,] <- res2$ls.counts

#i=1
#> res2$par
#     var1      var2      var3 
#200.33274  11.08894 164.28566 
#> res2$iter
#[1] 8

  # procedure 3: ML-FD Newton-Raphson experimental version
  # with optimized step length if 'step' is NULL

  res3 <- maxlik.fd.scoring(m = m, step = step, 
    information = "mix", ls = ls, barrier = bar,
    control = list(maxit = maxiter, tol = tol, trace = TRUE))

  Mres$proc3$M[i,] <- res3$par
  Mres$proc3$iter[i] <- res3$iter
  Mres$proc3$paths[seq(res3$iter+1),,i] <- res3$Mpars
  if (length(res3$ls.counts) == 2)
    Mres$proc3$lscounts[i,] <- res3$ls.counts

#i=1
#> res3$par
#     var1      var2      var3 
#200.33274  11.08894 164.28566 
#> res3$iter
#[1] 8

  # procedure 4: ML-FD BFGS algorithm from 'optim()'
  # NOTE 'count' is recorded (not 'iter')

  res4 <- try(maxlik.fd.optim(m, barrier = bar, inf = 99999, 
    method = "BFGS", gr = "analytical"), silent = TRUE)

  if (!inherits(res4, "try-error"))
  {
    Mres$proc4$M[i,] <- res4$par
    Mres$proc4$iter[i] <- res4$iter[1]
  }

#i=1
#> res4$par
#    var1     var2     var3 
#2100.735 1389.820 3680.320 
#> res4$iter
#function gradient 
#     100      100 

  # procedure 5: ML-FD L-BFGS-B algorithm from 'optim()'

  res5 <- try(maxlik.fd.optim(m, barrier = bar, inf = 99999, 
    method = "L-BFGS-B", gr = "analytical"), silent = TRUE)

  # record the path followed by the optimization method
  # using the same stopping criterion as in 'maxlik.fd.scoring()'
  if (!inherits(res5, "try-error"))
  {
    Mpars <- m@pars
    conv <- FALSE
    r <- 0
    while (!conv)
    {
      optout <- maxlik.fd.optim(m = m, barrier = bar, inf = 99999, 
        method = "L-BFGS-B", gr = "analytical", 
        optim.control = list(trace = FALSE, maxit = r))
      Mpars <- rbind(Mpars, optout$par)
      #conv <- optout$convergence == 0
      if (sqrt(sum((Mpars[r+1,] - Mpars[r+2,])^2)) < tol || r > maxiter)
        conv <- TRUE
      r <- r + 1
    }

    Mres$proc5$M[i,] <- Mpars[r+1,]
    Mres$proc5$iter[i] <- r
    Mres$proc5$paths[seq(r+1),,i] <- Mpars
  }

#i=1
#> Mpars[r+1,]
#     var1      var2      var3 
#200.33488  11.08881 164.28346 
#> r
#[1] 31

  # trace simulation

  tmp <- rbind(na.omit(Mres$proc1$paths[,,i]), na.omit(Mres$proc2$paths[,,i]),
    na.omit(Mres$proc3$paths[,,i]), na.omit(Mres$proc5$paths[,,i]))
  rg <- apply(tmp, 2, range)

  options(warn = -1)
  tmp <- persp(rg[,1], rg[,2], cbind(rep(Inf, 2), rep(Inf, 2)),
    zlim = rg[,3], box = TRUE, axes = TRUE, border = NA,
    xlab = "var1", ylab = "var2", zlab = "var3",  
    expand = 0.8, ticktype = "detailed",
    theta = 30, phi = 25, main = paste("series", i))

  for (k in seq_along(plotid))
  {
    idk <- plotid[k]
    mp <- na.omit(Mres[[idk]]$paths[,,i])
    if (length(mp) > 0) 
    {
      a <- trans3d(mp[,1], mp[,2], mp[,3], tmp)
      for (j in seq(2, nrow(mp)))
      {
        arrows(a$x[j-1], a$y[j-1], a$x[j], a$y[j], 
          lty = 1, col = colors[k], length = 0.15, angle = arrow.angle[k])
      }
    }

    legtext <- c(legtext, 
      paste(plab[k], ": ", Mres[[idk]]$iter[i], " iter.", sep = ""))
  }

  legend(par("usr")[1], par("usr")[4]*1.5, legtext, cex = 1, adj = 0,
    col = colors, lty = 1, xjust = 0, yjust = 2.0, bty = "n", horiz = FALSE)
  legtext <- NULL
  options(warn = 0)

  #print(i)
  trace.iter <- 100 * (i / iter)
  if (trace.iter %% 10 == 0)
  {
    cat(paste(trace.iter, "% complete.\n", sep = ""))
    #print(colMeans(Mres$proc1$M, na.rm = TRUE))
    print(colMeans(Mres$proc2$M, na.rm = TRUE))
    #print(colMeans(Mres$proc3$M, na.rm = TRUE))
    #print(colMeans(Mres$proc4$M, na.rm = TRUE))
    #print(colMeans(Mres$proc5$M, na.rm = TRUE))
  }
}

# store results

#setwd("")
#save(Mres, file = "sim-llmseas-ml-fd.rda")

# summary of results

colMeans(Mres$proc0$M)
colMeans(Mres$proc1$M)
colMeans(Mres$proc2$M)
colMeans(Mres$proc3$M)
colMeans(Mres$proc4$M)
colMeans(Mres$proc5$M)

summary(Mres$proc0$M)
summary(Mres$proc1$M)
