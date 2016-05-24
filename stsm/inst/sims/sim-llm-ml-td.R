#
# Simulation experiments:
# the series in the data set 'llm' are fit by 
# maximizing the time domain log-likelihood function (ML-TD)
#
# the following procedures are considered:
# procedure 0: 'StructTS()'
# procedure 1: ML-TD scoring algorithm
# procedure 2: ML-TD BFGS algorithm from 'optim()'
# procedure 3: ML-TD L-BFGS-B algorithm from 'optim()'
# procedure 5: ML-TD L-BFGS-B, 'var1' is concentrated out of the likelihood function
#
# some of the results obtained here are reported in Table 1
# in the vignette of the 'stsm' package
#

##NOTE added in optim() convergence if (... || r > maxiter)

library("stsm")

# load data

##NOTE
# requires data set generated in file "datagen-llm.R" in the 
# same folder as this file

iter <- ncol(llm)

# initial parameter values

initpars <- c(var1 = 1, var2 = 1)
nop <- NULL

# scoring algorithm parameters

tol <- 0.001
maxiter <- 100 #250
step <- NULL # (ignored for 'maxlik.fd.optim()')

# line search parameters 
# they are ignored if 'step' is fixed to a numeric, e.g. 'step <- 1'

ls <- list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1)
#ls <- list(type = "brent.fmin", tol = .Machine$double.eps^0.25, cap = 1)

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

plotid <- c(3, 1) + 1
#plab <- names(Mres)[plotid]
plab <- c("StructTS", "Scoring", "BFGS", "L-BFGS-B", "", "L-BFGS-B clik")[plotid]
colors <- c("red", "blue", "green")
arrow.angle <- c(30, 30, 30)
legtext <- NULL

# begin main loop

for (i in seq_len(iter))
{
  m <- stsm.model(model = "local-level", y = llm[,i], 
    pars = initpars, nopars = nop, transPars = transP)

  # procedure 0: ML-TD StructTS()
 
  res0 <- StructTS(m@y, type = "level")$coef[c(2, 1)]

  Mres$proc0$M[i,] <- res0

#i=1
#   epsilon      level 
#1742.59713   29.99584 

  # procedure 1: ML-TD scoring algorithm
  # with optimized step length if 'step' is NULL

  res1 <- maxlik.td.scoring(m = m, step = step, 
    KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
    ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
    control = list(maxit = maxiter, tol = tol, trace = TRUE))

  Mres$proc1$M[i,] <- res1$par
  Mres$proc1$iter[i] <- res1$iter
  Mres$proc1$paths[seq(res1$iter+1),,i] <- res1$Mpars
  if (length(res1$ls.counts) == 2) # 'optimize()' does not return this
    Mres$proc1$lscounts[i,] <- res1$ls.counts

#i=1
#> res1$par
#    var1     var2 
#1742.597   29.996 
#> res1$iter
#[1] 10

  # procedure 2: ML-TD BFGS
  # NOTE 'count' is recorded (not 'iter')

  res2 <- maxlik.td.optim(m, KF.version = "KFKSDS", 
    KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
    barrier = bar, inf = 99999, method = "BFGS", gr = "numerical")

  Mres$proc2$M[i,] <- res2$par
  Mres$proc2$iter[i] <- res2$iter[1]
  #Mres$proc2$paths[seq(res2$iter+1),,i] <- res2$Mpars

#i=1
#> res2$par
#    var1     var2 
#7767.246 3806.246 
#> res2$iter
#function gradient 
#     101      100 
#the problem may be that it starts from values close to 0
#and as bounds are not anyhow imposed unfeasible values are tried 
#at the first iterations
#adding a barrier term (trial and error values) did not make any improvements

  # procedure 3: ML-TD L-BFGS-B algorithm from 'optim()'

  res3 <- try(maxlik.td.optim(m, KF.version = "KFKSDS", 
    KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
    barrier = bar, inf = 99999, method = "L-BFGS-B", gr = "numerical"), 
    silent = TRUE)

  # record the path followed by the optimization method
  # using the same stopping criterion as in 'maxlik.fd.scoring()'
  if (!inherits(res3, "try-error"))
  {
    Mpars <- m@pars
    conv <- FALSE
    r <- 0
    while (!conv)
    {
      optout <- maxlik.td.optim(m, KF.version = "KFKSDS", 
        KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
        barrier = bar, inf = 99999, method = "L-BFGS-B", gr = "numerical", 
        optim.control = list(trace = FALSE, maxit = r))
      Mpars <- rbind(Mpars, optout$par)
      #conv <- optout$convergence == 0
      if (sqrt(sum((Mpars[r+1,] - Mpars[r+2,])^2)) < tol || r > maxiter)
        conv <- TRUE
      r <- r + 1
    }

    Mres$proc3$M[i,] <- Mpars[r+1,]
    Mres$proc3$iter[i] <- r
    Mres$proc3$paths[seq(r+1),,i] <- Mpars
  }

#i=1
#      var1       var2 
#1742.59717   29.99584 
#> res3$iter
#function gradient 
#      37       37 
#> r
#[1] 34

  # procedure 4: procedure in package 'KFAS' (currently not supported,
  # the last version of package has changed the interface)
#i=1
#> res4$opt$par
#      var1       var2 
#1742.59712   29.99586 
#> res4$opt$counts
#function gradient 
#      37       37 

  # procedure 5: ML-TD L-BFGS-B 
  # 'var1' is concentrated out of the likelihood function

  mc <- stsm.model(model = "local-level", y = llm[,i], 
    #pars = initpars[-2], cpar = initpars[2], nopars = nop, transPars = transP)
    pars = initpars[-1], cpar = initpars[1], nopars = nop, transPars = transP)

  res5 <- try(maxlik.td.optim(mc, KF.version = "KFKSDS", 
    KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
    barrier = bar, inf = 99999, method = "L-BFGS-B", gr = "numerical"), 
    silent = TRUE)
  
  if (!inherits(res5, "try-error"))
    Mres$proc5$M[i,] <- coef(res5)

  # the code below is very time consuming
  # res5$iter[1] can be used instead although it does not
  # use the same stopping rule set in the other procedures
  #
  # record the path followed by the optimization method
  # using the same stopping criterion as in 'maxlik.fd.scoring()'
  if (!inherits(res5, "try-error"))
  {
    Mpars <- c(mc@pars, mc@cpar)
    conv <- FALSE
    r <- 0
    while (!conv)
    {
      optout <- maxlik.td.optim(mc, KF.version = "KFKSDS", 
        KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
        barrier = bar, inf = 99999, method = "L-BFGS-B", gr = "numerical", 
        optim.control = list(trace = FALSE, maxit = r))
      Mpars <- rbind(Mpars, 
        c(get.cpar(optout$model, TRUE), get.pars(optout$model, TRUE)))
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
#pars = initpars[-1], cpar = initpars[1]
#> res5$par
#      var2 
#0.01782615 
#> Mpars[r+1,]
#      var1       var2 
#1742.59717   29.99584 
#> c(res5$model@cpar, res5$model@cpar * res5$par)
#      var1       var1 
#1724.32544   30.73809 
#> coef(res5)
#      var1       var2 
#1724.32544   30.73809 
#> r
#[1] 34
#> res5$iter
#function gradient 
#      30       30 

  # trace simulation

  tmp <- rbind(na.omit(Mres$proc1$paths[,,i]), na.omit(Mres$proc3$paths[,,i]), 
    na.omit(Mres$proc5$paths[,,i]))
  rg <- apply(tmp, 2, range)

  options(warn = -1)
  plot(rg[,1], rg[,2], type = "n", xlab = "var1", ylab = "var2", 
    main = paste("series", i))

  for (k in seq_along(plotid))
  {
    idk <- plotid[k]
    mp <- na.omit(Mres[[idk]]$paths[,,i])
    if (length(mp) > 0) 
    {
      for (j in seq(2, nrow(mp)))
      {
        arrows(mp[j-1,1], mp[j-1,2], mp[j,1], mp[j,2], 
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
    #print(colMeans(Mres$proc5$M, na.rm = TRUE))
  }
}

# store results

#setwd("")
#save(Mres, file = "sim-llm-ml-td.rda")

# summary of results

colMeans(Mres$proc0$M)
colMeans(Mres$proc1$M)
colMeans(Mres$proc2$M)
colMeans(Mres$proc3$M)
colMeans(Mres$proc4$M)
colMeans(Mres$proc5$M)

> colMeans(Mres$proc0$M)
      var1       var2 
1609.61466   99.09333 
> colMeans(Mres$proc1$M)
      var1       var2 
1609.02088   99.15593 
> colMeans(Mres$proc3$M)
      var1       var2 
1609.32534   99.14996 
> colMeans(Mres$proc5$M)
     var1      var2 
1591.8766  100.4652 

mean(Mres$proc0$iter)
mean(Mres$proc1$iter)
mean(Mres$proc2$iter)
mean(Mres$proc3$iter)
mean(Mres$proc4$iter)
mean(Mres$proc5$iter)

> mean(Mres$proc1$iter)
[1] 8.507
> mean(Mres$proc3$iter)
[1] 33.624
> mean(Mres$proc4$iter)
