## ----initalization, echo = FALSE, results = 'hide', message = FALSE----
library("knitr")
opts_chunk$set(concordance = TRUE, tidy = FALSE,
               strip.white = TRUE, dev = 'pdf',
               prompt = TRUE)

# Change command prompt appearance
#render_sweave() # Render output as Sweave, no syntax highlight etc.
options(replace.assign = TRUE, width = 70,
        continue = "+  ", prompt = "> ",
        useFancyQuotes = FALSE)

read_chunk('GMCM-Standalone.R')  # Read chunks in the script to input here

## ----chunk_1, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
#
# General initial parameters, packages, and auxiliary functions used throughout
#

rm(list = ls()) # Remove any objects in workspace (ensure a clean run)

# Set working directory if nessesary
# setwd("<the directory of this file>")

# Should all saved data be recomputed?
recompute <- FALSE
# Set recompute to TRUE if saved data is to be recomputed and not just loaded.
# WARNING! Setting will increase the running time to ~ 5 days on a regular
# laptop (depending on the use of parallel computations).
# Alternatively, "saved.RData" can be deleted or particular objects can
# be removed using rm() to force the recomputation

# We load all saved objects that are computationally time consuming
if (file.exists("saved.RData")) {
  load("saved.RData")
}


# Global plotting parameters used in some plots
axiscex <- 1.7

#
# Loading needed libraries
#

#install.packages(c("GMCM", "idr", "Hmisc", "RColorBrewer", "jpeg"))
library("GMCM")
library("idr")
library("Hmisc")
library("RColorBrewer")
library("jpeg")

# To enable multicore processing, change the %do% to %dopar% and
# uncomment the relevant packages below.
library("foreach")
# library("doMC") #library("doParallel") # Use this package on windows
# registerDoMC(detectCores())

# Defining auxiliary functions used
# Pretty big number formatting
prettyN <- function (x) {
  prettyNum(x, big.mark = "{,}", scientific = FALSE)
}

# Function for appending objects to an existing .RData file
# I used modified version of a function created by user "flodel":
# http://stackoverflow.com/questions/11813096/updating-an-existing-rdata-file
resave <- function(..., list = character(), file) {
  if (!file.exists(file)) { # If file does not exists resave functions as save
    save(..., list = list, file = file)
  }
  previous  <- load(file)  # Returns the loaded object names
  var.names <- c(list, as.character(substitute(list(...)))[-1L])
  for (var in var.names) {
    assign(var, get(var, envir = parent.frame()))
  }
  save(list = unique(c(previous, var.names)), file = file)
}


################################################################################
# Create Figure 1 of Section "2. Gaussian mixture copula models"
################################################################################

# Setting a seed for the simulation of data
set.seed(120312)

# Choosing some arbitrary marignal distributions for the general GMCM
inv.M1 <- function(u) {
  M1 <- function(u1) qbeta(u1, shape1 = 0.5, shape2 = 0.9)
  M2 <- function(u2) qchisq(u2, df = 5)
  u[,1] <- M1(u[,1])
  u[,2] <- M2(u[,2])
  return(u)
}

# Choosing some arbitrary marignal distributions for the special GMCM
inv.M2 <- function(u) {
  M1 <- function(u1) qbeta(u1, shape1 = 0.5, shape2 = 0.9)
  M2 <- function(u2) qbeta(u2, shape1 = 1.3, shape2 = 1.1)
  u[,1] <- M1(u[,1])
  u[,2] <- M2(u[,2])
  return(u)
}

# Simulate general GMCM data
theta1 <- list(m = 3, d = 2,
               pie = c(1/6,3/6,2/6),
               mu = list(comp1 = c(0,0),
                         comp2 = c(2,-1)*2,
                         comp3 = c(-1,-0.5)*2),
               sigma = list(comp1 = cbind(c(1,0),c(0,1)),
                            comp2 = cbind(c(2,1.78),c(1.78,2)),
                            comp3 = cbind(c(1,-0.53),c(-0.53,1))))
data1 <- SimulateGMCMData(n = 10000, theta = theta1)
z1 <- data1$z
u1 <- data1$u
x1 <- inv.M1(u1)

# Simulate special GMCM data
par <- c(0.8, mu = 2, sigma = 1, rho = 0.8)
data2 <- SimulateGMCMData(n = 10000, d = 2, par = par)
theta2 <- data2$theta
z2 <- data2$z
u2 <- data2$u
x2 <- inv.M2(u2)

#
# Create Figure 1
#

jpeg("Figure1.jpg", height = 2*7*0.5, width = 3*7*0.5, units = "in", res = 100)
{
  lab.cex <- 1
  # Setting plotting parameters
  par(mgp = c(2.3,0.8,0),
      oma = c(0,2,0,0)+0.1,
      mar = c(3, 3.3, 2, 0.1),
      xaxs = "i", yaxs = "i",
      mfrow = c(2,3),
      xpd = TRUE)

  col1 <- c("white", "grey60")
  col2 <- c("black", "orange", "steelblue")
  pchs <- c(1,3,4)
  pcex <- 0.6
  labs <- pretty(seq(0,1))

  # Panel A: Observed
  plot(x1, axes = FALSE, type = "n", xlab = "", ylab = "")
  points(x1, col = col2[data1$K], pch = pchs[data1$K], cex = pcex)
  mtext(expression(paste(x[1])), 1, line = 2, cex = lab.cex)
  mtext(expression(paste(x[2])), 2, line = 2, cex = lab.cex)
  axis(1, labels = FALSE, at = axTicks(1), lwd = axiscex)
  axis(2, labels = FALSE, at = axTicks(2), lwd = axiscex)
  mtext("A", line = 0.9, adj = -0.05, cex = 1, font = 2)
  mtext("Observed process", line = 0.5, cex = 1, font = 2)
  mtext("General GMCM", line = 3.5, side = 2, cex = 1, font = 2)

  # Panel B: Copula density
  plot(1, type = "n", axes = FALSE, xlim = c(0,1), xlab = "",
       ylim = c(0,1), ylab = "")#, asp = 1)
  points(u1, col = col2[data1$K], pch = pchs[data1$K], cex = pcex)
  mtext(expression(paste(u[1])), 1, line = 2, cex = lab.cex)
  mtext(expression(paste(u[2])), 2, line = 2, cex = lab.cex)
  axis(1, at = labs, lwd = axiscex)
  axis(2, at = labs, lwd = axiscex)
  mtext("B", line = 0.9, adj = -0.05, cex = 1, font = 2)
  mtext("Copula process", line = 0.5, cex = 1, font = 2)

  # Panel C: Latent
  plot(1, type = "n", axes = FALSE, xlim = range(z1[,1]),
       xlab = "", ylim = range(z1[,2]), ylab = "")#, asp = 1)
  mtext(expression(paste(z[1])), 1, line = 2, cex = lab.cex)
  mtext(expression(paste(z[2])), 2, line = 2, cex = lab.cex)
  points(z1, col = col2[data1$K], pch = pchs[data1$K], cex = pcex)
  axis(1, lwd = axiscex)
  axis(2, lwd = axiscex)
  mtext("C", line = 0.9, adj = -0.05, cex = 1, font = 2)
  mtext("Latent process", line = 0.5, cex = 1, font = 2)


  col2 <- c("steelblue", "black", "orange")

  # Panel D: Observed process
  plot(x2, axes = FALSE, type = "n", xlab = "", ylab = "")
  points(x2, col = col2[data2$K], pch = pchs[data2$K], cex = pcex)
  mtext(expression(paste(x[1])), 1, line = 2, cex = lab.cex)
  mtext(expression(paste(x[2])), 2, line = 2, cex = lab.cex)
  axis(1, labels=FALSE, at = axTicks(1), lwd = axiscex)
  axis(2, labels=FALSE, at = axTicks(2), lwd = axiscex)
  mtext("D", line = 0.9, adj = -0.05, cex = 1, font = 2)
  mtext("Special GMCM", line = 3.5, side = 2, cex = 1, font = 2)

  # Panel E: Copula density
  plot(1, type = "n", axes = FALSE, xlim = c(0,1), xlab = "",
       ylim = c(0,1), ylab = "")#, asp = 1)
  points(u2, col = col2[data2$K], pch = pchs[data2$K], cex = pcex)
  mtext(expression(paste(u[1])), 1, line = 2, cex = lab.cex)
  mtext(expression(paste(u[2])), 2, line = 2, cex = lab.cex)
  axis(1, at = labs, lwd = axiscex)
  axis(2, at = labs, lwd = axiscex)
  mtext("E", line = 0.9, adj = -0.05, cex = 1, font = 2)

  # Panel F: Latent process
  plot(1, type = "n", axes = FALSE, xlim = range(z2[,1]), xlab = "",
       ylim = range(z2[,2]), ylab = "")#, asp = 1)
  points(z2, col = col2[data2$K], pch = pchs[data2$K], cex = pcex)
  mtext(expression(paste(z[1])), 1, line = 2, cex = lab.cex)
  mtext(expression(paste(z[2])), 2, line = 2, cex = lab.cex)
  axis(1, lwd = axiscex)
  axis(2, lwd = axiscex)
  mtext("F", line = 0.9, adj = -0.05, cex = 1, font = 2)

}
dev.off()


################################################################################
# The following performs the speed tests of section
# "4.2. Runtime and technical comparison"
################################################################################

# Set seed for the simulation of data
set.seed(1)

# Setting parameters for the simulation
# Start values:
alpha <- 0.5
mu    <- 2.5
sigma <- 0.5
rho   <- 0.8

speed.eps      <- 0.001  # Convergence epsilon
max.ite        <- 1000   # Maximum iterations
speed.par      <- c(0.7, 2, 1, 0.9)  # True parameters
n.observations <- c(1000, 10000, 100000) # Simulated observations


# We only do the computation is the saved speed.res is not present
# or if recompute is set to TRUE.
if (!exists("speed.res") | recompute) {

  speed.res <- foreach(n = n.observations, .combine = "rbind",
                       .packages = c("GMCM", "idr")) %do% {

    # Simulating data
    x <- SimulateGMCMData(n = n, par = speed.par)

    # Ranking and scaling
    u <- Uhat(x$u)

    # Li et al. (2011) PEM algorithm and results
    idr.time <- system.time({
      idr.out <- est.IDR(u, mu, sigma, rho, alpha,
                         eps = speed.eps, max.ite = max.ite)
    })

    if (length(idr.out$loglik.trace) == 2) {
      stop("est.IDR only took two interations")
    }

    li.res <- unlist(idr.out$para[c("p", "mu", "sigma", "rho")])

    # GMCM package PEM timing
    gmcm.time.pem <- system.time({
      my.res <- fit.meta.GMCM(u,  init.par = c(alpha, mu, sigma, rho),
                              method = "PEM", trace.theta = TRUE,
                              max.ite = max.ite, verbose = FALSE)
    })

    # GMCM package NM timing
    gmcm.time.nm <- system.time({
      my.res2 <- fit.meta.GMCM(u, init.par = c(alpha, mu, sigma, rho),
                               method = "NM", trace.theta = TRUE,
                               max.ite = max.ite, verbose = FALSE)
    })

    # Getting convergence information
    aa <- c(it = length(idr.out$loglik.trace),   idr.time[3])
    bb <- c(it = ncol(my.res[[2]]$loglik.tr),    gmcm.time.pem[3])
    cc <- c(it = unname(my.res2[[2]]$counts[1]), gmcm.time.nm[3])

    res <- rbind("idr-pkg PEM" = aa, "gmcm-pkg PEM" = bb, "gmcm-pkg NM" = cc)
    res <- cbind(res, "$s/n$" = res[,2]/res[,1])
    res <- cbind(res, "Rel.\\ speed" = res[,3]/res[3,3])
    res <- cbind(rbind(unlist(idr.out$para)[c(1,3,4,2)],
                       my.res[[1]], my.res2[[1]]), res)
    res <- cbind(n.obs = n, res)
    # speed.res <- rbind(speed.res, res)
    cat("Speed test with", n, "observations finished\n"); flush.console()

    return(res)
  }

  resave(speed.res, file = "saved.RData")
}

# Print the results
#print(speed.res)

################################################################################
# The following performs the comparison of the different fitting procedures
# in section "4.2. Runtime and technical comparison"
# NOTE: Some of the computations are VERY time consuming. I.e., the computation
# time is measured in days for some as warned below.
################################################################################

# The structure of this section is as follows:
# 1. Some initalization of parameters and aux. functions
# 2. Simulate the "n.sims" datasets and starting parameters
# 3. Fit the special model using NM to the generated datasets
# 4. Fit the special model using SANN ...
# 5. Fit the special model using L-BFGS ...
# 6. Fit the special model using L-BFGS-B ...
# 7. Fit the special model using PEM ...
# 8. Fit the special model using PEM from the idr-package ...
# 9. Format the results and create Figure 4

# Setting seed to ensure reproducibiity of results
set.seed(1)

# Setting parameters of simulations
n.sims  <- simulation.n.sims  <- 1000  # Number of simulated datasets
n.obs   <- simulation.n.obs   <- 10000 # Number of observations in each dataset
max.ite <- simulation.max.ite <- 2000  # Maximum number of interations

d <- 2 # Dimension of data

# Define the true parameters
true.par  <- simulation.true.par <-
  c(alpha1 = 0.90, mu = 3, sigma = 2, rho = 0.5)

# List the available fitting methods
methods   <- c("NM", "SANN", "L-BFGS", "L-BFGS-B", "PEM (GMCM)", "PEM (idr)")
n.methods <- length(methods)

# Define an error and warning handler
# The function is based off
# http://svn.r-project.org/R/trunk/src/library/base/demo/error.catching.R
tryCatch.W.E <- function(expr) {
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  e.handler <- function(e){
    e
  }
  value <-
    withCallingHandlers(tryCatch(expr, error = e.handler), warning = w.handler)
  if ("error" %in% class(value)) {
    E <- value
    value <- NULL
  } else {
    E <- NULL
  }
  return(list("value" = value, "error" = E, "warning" = W))
}

#
# Simulate GMCM data and starting parameters
# Approximate runtime: 13 seconds
#

cat("Simulating", n.sims, "datasets and start parameters... "); flush.console()

sim.tmp <- foreach(i = seq_len(n.sims), .packages = "GMCM") %do% {
  # Draw random starting values and simulate from the GMCM
  return(list(startval = c(rbeta(1,1.5,1.5), abs(rnorm(1,3,1)),
                           rchisq(1,1), rbeta(1,1.5,1.5)),
              simdata = SimulateGMCMData(n = n.obs, par = true.par, d = d)))
}

simulation.start.par <- t(sapply(sim.tmp, "[[", "startval"))
colnames(simulation.start.par) <- names(true.par)

simulation.data <- lapply(sim.tmp, "[[", "simdata")

rm(sim.tmp)

# To save space in the RData file. Data is not saved. Here.
cat("Done.\n")


#
# Fit using Nelder-Mead
#

# Only do the computation if stored data is not available or if forced by
# recompute
# Approximate runtime: 15 minutes (6 mins in parallel)
if (!exists("simulation.res.NM") | recompute) {
  st.tot <- proc.time()
  simulation.res.NM <-
    foreach(i = seq_len(n.sims), .packages = "GMCM", .inorder = FALSE) %do% {

      x <- Uhat(simulation.data[[i]]$u)
      K <- simulation.data[[i]]$K
      st <- proc.time()
      tmp <-
        tryCatch.W.E(fit.meta.GMCM(u = x, init.par = simulation.start.par[i, ],
                                   method = "NM", trace.theta = TRUE,
                                   positive.rho = TRUE, max.ite = max.ite,
                                   verbose = FALSE, reltol = 1e-4))
      time <- (proc.time()-st)[3]
      par  <- tmp$value[[1]]
      ite  <- tmp$value[[2]]$counts[[1]]
      Khat <- (get.IDR(x, par)$idr <= 0.5) + 1
      acc  <- sum(diag(table(K, Khat)))/length(K)

      return(list("par" = par, "error" = tmp$error, "warning" = tmp$warning,
                  "ite" = ite, "acc" = acc, "time" = time))
    }

  resave(simulation.res.NM, file = "saved.RData")
  cat("\n## Nelder-Mead done ##\n\n")
  cat("Time:", (proc.time()-st.tot)[3] %/% 60 ,"min ellapsed\n")
}

#
# Fit using SANN
#

set.seed(65) # Since the SANN procedure is stochastic

# Approximate runtime: 309 minutes ~ 5 hours (WARNING !!!) (2.2 hrs with par)
#
if (!exists("simulation.res.SANN") | recompute) {
  st.tot <- proc.time()
  simulation.res.SANN <-
    foreach(i = seq_len(n.sims), .packages = "GMCM", .inorder = FALSE) %do% {

      x <- Uhat(simulation.data[[i]]$u)
      K <- simulation.data[[i]]$K
      st <- proc.time()
      tmp <-
        tryCatch.W.E(fit.meta.GMCM(u = x, init.par = simulation.start.par[i, ],
                                   method = "SANN", trace.theta = TRUE,
                                   positive.rho = TRUE, max.ite = 3000,
                                   verbose = FALSE))
      time <- (proc.time()-st)[3]
      par <- tmp$value[[1]]
      ite <- tmp$value[[2]]$counts[[1]]
      if(!is.null(par)) {
        Khat <- (get.IDR(x, par)$idr <= 0.5)+1
        acc <- sum(diag(table(K, Khat)))/length(K)
      } else {
        Khat <- NULL
        acc <- NULL
      }

      return(list("par" = par,
                  "error" = tmp$error,
                  "warning" = tmp$warning,
                  "ite" = ite,
                  "acc" = acc, "time" = time))
    }

  resave(simulation.res.SANN, file = "saved.RData")
  cat("\n## SANN done ##\n\n")
  cat("Time:", (proc.time()-st.tot)[3] %/% 60 ,"min ellapsed\n")
}

#
# Fit using L-BFGS
#

# Approximate runtime: 85 minutes ~ 1.5 hours (!!!)

if (!exists("simulation.res.LBFGS") | recompute) {
  st.tot <- proc.time()
  simulation.res.LBFGS <-
    foreach(i = seq_len(n.sims), .packages = "GMCM", .inorder = FALSE) %do% {

      x <- Uhat(simulation.data[[i]]$u)
      K <- simulation.data[[i]]$K
      st <- proc.time()
      tmp <-
        tryCatch.W.E(fit.meta.GMCM(u = x, init.par = simulation.start.par[i, ],
                                   method = "L-BFGS", trace.theta = TRUE,
                                   positive.rho = TRUE, max.ite = max.ite,
                                   verbose = FALSE, factr = 1e-4))
      time <- (proc.time()-st)[3]
      par <- tmp$value[[1]]
      ite <- tmp$value[[2]]$counts[[1]]
      if(!is.null(par)) {
        Khat <- (get.IDR(x, par)$idr<=0.5)+1
        acc <- sum(diag(table(K, Khat)))/length(K)
      } else {
        Khat <- NULL
        acc <- NULL
      }
      return(list("par" = par, "error" = tmp$error,
                  "warning" = tmp$warning, "ite" = ite,
                  "acc" = acc, "time" = time))
    }
  cat("\n## L-BFGS done ##\n\n")
  resave(simulation.res.LBFGS, file = "saved.RData")
}

#
# Fit using L-BFGS-B
#

# Approximate runtime: 96 minutes ~ 1.5 hours (half if parallel)

if (!exists("simulation.res.LBFGSB") | recompute) {
  st.tot <- proc.time()
  simulation.res.LBFGSB <-
    foreach(i = seq_len(n.sims), .packages = "GMCM", .inorder = FALSE) %do% {

      x <- Uhat(simulation.data[[i]]$u)
      K <- simulation.data[[i]]$K
      st <- proc.time()
      tmp <-
        tryCatch.W.E(fit.meta.GMCM(u = x, init.par = simulation.start.par[i, ],
                                   method = "L-BFGS-B", trace.theta = TRUE,
                                   positive.rho = TRUE, max.ite = max.ite,
                                   verbose = FALSE, factr = 1e-4))
      time <- (proc.time()-st)[3]
      par <- tmp$value[[1]]
      ite <- tmp$value[[2]]$counts[[1]]
      if(!is.null(par)) {
        Khat <- (get.IDR(x, par)$idr<=0.5)+1
        acc <- sum(diag(table(K, Khat)))/length(K)
      } else {
        Khat <- NULL
        acc <- NULL
      }

      return(list("par" = par, "error" = tmp$error, "warning" = tmp$warning,
                  "ite" = ite, "acc" = acc, "time" = time))
    }
  cat("\n\n## L-BFGS-B done ##\n\n\n")
  resave(simulation.res.LBFGSB, file = "saved.RData")
}


#
# Fit using PEM in GMCM
#

# Approximate runtime: 175 minutes ~ 3 hours (!!!)

if (!exists("simulation.res.PEM") | recompute) {
  st.tot <- proc.time()
  simulation.res.PEM <-
    foreach(i = seq_len(n.sims), .packages = "GMCM", .inorder = FALSE) %do% {
      x <- Uhat(simulation.data[[i]]$u)
      K <- simulation.data[[i]]$K
      st <- proc.time()
      tmp <-
        tryCatch.W.E(fit.meta.GMCM(u = x, init.par = simulation.start.par[i, ],
                                   method = "PEM", trace.theta = TRUE,
                                   positive.rho = TRUE, max.ite = max.ite,
                                   verbose = FALSE, eps = 1e-4))
      time <- (proc.time()-st)[3]
      par <- tmp$value[[1]]
      ite <- length(tmp$value[[2]]$theta.tr)
      if(!is.null(par)) {
        Khat <- (get.IDR(x, par)$idr<=0.5)+1
        acc <- sum(diag(table(K, Khat)))/length(K)
      } else {
        Khat <- NULL
        acc <- NULL
      }
      return(list("par" = par, "error" = tmp$error, "warning" = tmp$warning,
                  "ite" = ite, "acc" = acc, "time" = time))
    }
  cat("\n\n## PEM in GMCM done ##\n\n\n")
  resave(simulation.res.PEM, file = "saved.RData")
}

#
# Fit using PEM in idr
#

# WARNING: if the file is missing or recompute is TRUE,
# then running the following if condition will take a long time!

# Approximate runtime: 5810 minutes ~ 4 days (!!! WARNING !!!)

if (!exists("simulation.res.PEMidr") | recompute) {
  st.tot <- proc.time()
  simulation.res.PEMidr <-
    foreach(i = seq_len(n.sims), .packages = "idr", .inorder = FALSE) %do% {
      x <- Uhat(simulation.data[[i]]$u)
      K <- simulation.data[[i]]$K
      start.par <- simulation.start.par[i, ]
      st <- proc.time()
      tmp <-
        tryCatch.W.E(est.IDR(x = x, mu = start.par[2], sigma = start.par[3],
                             rho = start.par[4], p = 1-start.par[1],
                             eps = 1e-4, max.ite = max.ite))
      time <- (proc.time()-st)[3]
      par <- unlist(tmp$value$para)[c("p","mu","sigma","rho")]
      par[1] <- 1 - par[1]
      ite <- length(tmp$value$loglik.trace)
      if(!is.null(par)) {
        Khat <- (get.IDR(x, par)$idr<=0.5)+1
        acc <- sum(diag(table(K, Khat)))/length(K)
      } else {
        Khat <- NULL
        acc <- NULL
      }
      return(list("par" = par, "error" = tmp$error, "warning" = tmp$warning,
                  "ite" = ite, "acc" = acc, "time" = time))
    }
  cat("\n## PEM in idr done ##\n\n")
  resave(simulation.res.PEMidr, file = "saved.RData")
}



#
# Extract the information from the fitting procedures and
# create Figure 4
#

# Function to get data from the lists
Get <- function(x, y) {
  x <- lapply(x, "[[", y)
  x <- x[!sapply(x, is.null)]
  if (identical(x, list())) {
    integer(0)
  } else {
    simplify2array(x)
  }
}

# Function to handle some formatting
formatForFig <- function (x, digits = 2) {
  m1 <- mean(x, na.rm = TRUE)
  m2 <- sd(x, na.rm = TRUE)
  m1 <- formatC(m1, format = "f", digits = digits)
  m2 <- formatC(m2, format = "f", digits = digits)
  paste(m1, "\n(", m2, ")", sep = "")
}

# Initialize object to store information about the number of iterations to
# reach convergence, warning messages, and error messages.
simulation.ites  <- vector("list", length(methods))
warning.messages <- vector("list", length(methods))
error.messages   <- vector("list", length(methods))
names(simulation.ites)  <- rev(methods) # Naming the lists
names(warning.messages) <- rev(methods)
names(error.messages)   <- rev(methods)

# Parameter estimates
jpeg("Figure4.jpg", height = 7, width = 2*7, units = "in", res = 100)
{
  # Setting plotting parameters
  par(mfrow = c(1,6),
      mar = c(4.1, 0, 3.1, 0),
      oma = c(0,6,0,2),
      cex = 1.2,
      mgp = c(2.3,0.8,0))

  cols <- brewer.pal(n = length(methods), "Dark2") # Defining colours

  # Defining some helpful objects
  estimate  <- c("alpha[1]", "mu", "sigma", "rho")
  estimate2 <- c("hat(alpha)[1]", "hat(mu)", "hat(sigma)", "hat(rho)")
  ylims <- rbind(c(0, -10, 0, -0.5, 0), c(1, 10, 6, 1, 1))

  # Plotting the six "panels"
  for (i in 1:6) {
    if (i <= 4) {
      plot(1, ylim = c(0.5,6.5) ,
           type = "n", xlim = ylims[,i], axes = FALSE,
           ylab = "",
           xlab = parse(text = estimate2[i]))
      header <- parse(text=paste(estimate,"==",true.par,sep=""))[i]
      title(main = header)
      if (i == 1)
        axis(2, at = length(methods):1, lab = methods, las = 2, lwd = axiscex)
      abline(v = true.par[i], col = "dark gray")
      box(lwd = axiscex)
    } else {
      if (i == 6) {
        plot(1, xlim = c(0,2000), ylim = c(0.5,6.5), type = "n", axes = FALSE,
             xlab = "Iterations", ylab = "", main = "")
        axis(3, lwd = axiscex)
        box(lwd = axiscex)
      }
      if (i == 5) {
        plot(1, xlim = c(0,1), ylim = c(0.5,6.5), type = "n", axes = FALSE,
             xlab = "", ylab = "", xaxt = "n")
      }

    }

    if (i != 5) {
      if (i %in% 2:3) {
        lab <- c("", axTicks(1)[-c(1,length(axTicks(1)))], "")
      } else {
        lab <- axTicks(1)
      }
      axis(1, at = axTicks(1), labels = lab, lwd = axiscex)
    }

    for (j in length(methods):1) {

      data <- switch(j,
                     simulation.res.PEMidr,
                     simulation.res.PEM,
                     simulation.res.LBFGSB,
                     simulation.res.LBFGS,
                     simulation.res.SANN,
                     simulation.res.NM)


      not.failed <- sapply(data, function(x) is.null(x[["error"]]))
      pars <- simplify2array(lapply(data, "[[", "par")[not.failed])

      error.messages[[j]] <-
        sapply(data[!not.failed], function(x) x$error$message)

      warning.messages[[j]] <-
        Get(data,"warning")

      simulation.ites[[j]] <-
        n.ites <- Get(data, "ite")
      acc <- formatForFig(Get(data, "acc"))
      tim <- round(sum(Get(data, "time"), na.rm = TRUE) / 60)

      war <- length(Get(data,"warning"))/2
      err <- sum(!not.failed)

      if (i <= 4) {
        y <- pars[i, ]
        x <- rep(j, length(y))
        points(y,
               x + runif(length(x), -1, 1)*0.18,
               #x + rnorm(length(x), sd = 0.1),
               cex = 0.55, pch = 21,
               bg = paste(cols[j], "40", sep = ""),
               col = "#00000060")

        segments(median(y, na.rm = TRUE), y0 = j-0.3, y1 = j+0.3, lwd = 2)
      } else {
        if (i == 6)
          segments(n.ites, y0 = j-0.3, y1 = j+0.3,
                   col = paste(cols[j], "80", sep = ""))
        segments(median(n.ites, na.rm = TRUE), y0 = j-0.4, y1 = j+0.4,
                 lwd = 2)
        if (j == 5) text(j, 1000, "(3000)")
        if (i == 5) {
          labels <- c("Accuracy", "Time", "Warnings", "Errors")
          at <- seq(0.15,0.9,l = length(labels))
          axis(3, at = at, las = 2, labels = labels, tick = FALSE, hadj = 0.5,
               lwd = axiscex)
          axis(1, at = at, las = 2, labels = labels, tick = FALSE, hadj = 0.5,
               lwd = axiscex)

          text(at, j, c(acc, tim, war, err), cex = 0.8)
        }

      }
    }
  }
}
dev.off()

#
# Extract the number of times the maximum number interations was hit
#

simulation.max.ite.hit <-
  sapply(simulation.ites, function(x) sum(x %in% (simulation.max.ite+c(0,1))))

#print(simulation.max.ite.hit)


# Check errors and warnings
lapply(error.messages,
       function(x) table(if (length(x) == 0) ("Non") else (x)))

lapply(warning.messages,
       function(x) table(if (length(x) == 0) ("Non") else (unlist(x[1, ]))))


################################################################################
# The following performs the reproducibility analysis of u133VsExon
# ranked data using GMCM
# The results appear in section 5.1. Reproducibility of microarray results
################################################################################

# Load ranked data from GMCM package
data(u133VsExon) # The object is named u133VsExon

# Rename the data and Benjamini-Hochberg adjust the p-values
p.vals <- u133VsExon
adj.p.vals <- data.frame(u133 = p.adjust(p.vals$u133, method = "BH"),
                         exon = p.adjust(p.vals$exon, method = "BH"))

# Ranking and scaling using Uhat
u <- Uhat(1 - p.vals)  # Which is the same as Uhat(-log(p.vals))

# Define a list of three different intial parameters
init.par <- list(c(0.5, 1, 1, 0.5),
                 c(0.999, 5, 1, 0.5),
                 c(0.5, 1, 5, 0.999))

# Fitting GMCMs for each of the initial parameters
# We only do the computation if it is not present as and .RData object.
if (!exists("bcell.ite") | recompute) {
  log.lik <- -Inf
  for (i in 1:length(init.par)) {
    tmp.bcell.ite <- fit.meta.GMCM(u = u,
                                   init.par = init.par[[i]],
                                   method = "NM",
                                   trace.theta = TRUE)
    if (tmp.bcell.ite[[2]]$value > log.lik) {
      cat(i, "\n\n");
      bcell.ite <- tmp.bcell.ite
    }
  }
  resave(bcell.ite, file = "saved.RData")
}

# Extract the parameters and number of iterations
bcell.par <- bcell.ite[[1]]
bcell.n.ite <- unname(bcell.ite[[2]]$counts[1])

# Extractng the idr and IDR values
bcell.res <- get.IDR(x = u, par = bcell.par)

# Get marginally significant genes
exon.sig <- adj.p.vals$exon <= 0.05
u133.sig <- adj.p.vals$u133 <= 0.05

# Plotting the results and creating Figure 5
jpeg("Figure5.jpg", height = 7*0.5, width = 3*7*0.5, units = "in", res = 100)
{
  bcell.cols <- c("grey", "steelblue")#, "black")

  par(mgp = c(2.3,0.8,0),
      oma = c(0,0,0,0)+0.1,
      mar = c(3, 3.3, 2, 0.1),
      xaxs = "i", yaxs = "i",
      #mfrow = c(1,2),
      xpd = TRUE)
  layout(rbind(c(3,1,1,2,2,3)))

  # Determining colouring
  col.num <- (bcell.res$IDR <= 0.05) + 1

  labs <- pretty(seq(0,1))

  plot(u, main = "Ranked values",#"GC cells vs B-cells\nGMCM process",
       cex = 0.4, pch = 16, axes = FALSE, xlab = "", ylab = "",
       col = bcell.cols[col.num])
  mtext(expression(hat(z)[list(g, scriptstyle(U133))]), 1, line = 2, cex = 0.7)
  mtext(expression(hat(z)[list(g, scriptstyle(Exon))]), 2, line = 2, cex = 0.7)
  axis(1, at = labs, lwd = axiscex)
  axis(2, at = labs, lwd = axiscex)
  mtext("A", line = 0.9, adj = -0.05, cex = 1, font = 2)

  plot(GMCM:::qgmm.marginal(u = u, theta = meta2full(bcell.par, d = 2)),
       main = "Pseudo observations",#"GC cells vs B-cells\nEstimated GMM process",
       cex = 0.3, pch = 16, axes = FALSE, xlab = "", ylab = "",
       col = bcell.cols[col.num])
  mtext(expression(hat(z)[list(g, scriptstyle(U133))]), 1, line = 2, cex = 0.7)
  mtext(expression(hat(z)[list(g, scriptstyle(Exon))]), 2, line = 2, cex = 0.7)
  axis(1, lwd = axiscex)
  axis(2, lwd = axiscex)

  legend("topleft", pch = c(16, 16), col = bcell.cols, text.font = 2,
         legend = paste(c("Irreproducible genes\nn =",
                          "Reproducible genes\nn ="), table(bcell.res$Khat)),
         bty = "n", inset = 0.01, y.intersp = 2)
  par <- round(bcell.par,2)
  A <- bquote(alpha[1] == .(par[1]))
  B <- bquote(mu == .(par[2]))
  C <- bquote(sigma == .(par[3]))
  D <- bquote(rho == .(par[4]))
  legend("bottomright", bty = "n", inset = 0.01, y.intersp = 1,
         legend = c(A, B, C, D, expression()))
  mtext("B", line = 0.9, adj = -0.05, cex = 1, font = 2)
}
dev.off()



################################################################################
# The following performs the image segmentation in
# in section 5.3. Image segmentation using the general GMCM
################################################################################

# Define the number of colours to segment the data into
n.cols <- 10

fig7.file <- "./Figure7.jpg"


if (!file.exists(fig7.file)) {
  # The image was originally downloaded from:
  # http://totallyfreeimages.com/3167/STS-27,-Orbiter-Atlantis,-Liftoff
  # download.file(paste0("http://tfi.s3.amazonaws.com/previews/standard/d/6/",
  #                      "5fecf808486b157ea6971048bb29d6f4e58bd1d6.jpg"),
  #               destfile = file)
  #http://upload.wikimedia.org/wikipedia/commons/d/d3/Atlantis_taking_off_on_STS-27.jpg
  #   link <-
  #     paste0("http://i.space.com/images/i/000/010/634/original/",
  #            "shuttle-atlantis-lifts-off.jpg?1309297321")
  url <- "http://people.math.aau.dk/~abilgrau/GMCM/STS-27.jpg"
  download.file(url, destfile = fig7.file, method = "internal", mode = "wb")
}


# Read the image
pic <- readJPEG(fig7.file)

nn <- dim(pic)[1]
mm <- dim(pic)[2]

seg.gmcm <- seg.km <- pic  # To hold the different segmented pics

# Format the picture in a matrix with 3 columns
pic.rgbmat <- cbind(red   = c(pic[,,1]),
                    green = c(pic[,,2]),
                    blue  = c(pic[,,3]))

#
# Segmentation using GMCM
#

gmcm.file <- gsub(".jpg$", ".gmcm.RData", fig7.file)
# Approximate runtime: 50 minutes
if (!exists("seg.res.gmcm") | recompute) {
  best.loglik <- -Inf
  for (i in 1:10) {
    cat("i =", i, "\n"); flush.console()

    # Set a seed
    set.seed(i)

    # Choose starting parameters
    start.theta <- choose.theta(pic.rgbmat, m = n.cols,
                                iter.max = 10)

    # Fit the full GMCM
    seg.res.gmcm <-
      fit.full.GMCM(u = pic.rgbmat,
                    theta = start.theta,
                    method = "PEM",
                    max.ite = 100,
                    verbose = TRUE,
                    eps = 1e-4,
                    convergence.criterion = "GMCM")

    # Compute the log likelihood of the maximizing parameters
    loglik <- GMCM:::dgmcm.loglik(seg.res.gmcm, u = pic.rgbmat)

    # Update if a better estiamte is found
    if (best.loglik < loglik) {
      cat("Loglik updated (", loglik, ")\n", sep = ""); flush.console()
      best.i <- i
      best.loglik <- loglik
      best.seg.res.gmcm <- seg.res.gmcm
    }
  }

  # Save results
  seg.res.gmcm <- best.seg.res.gmcm
  resave(seg.res.gmcm, file = "saved.RData")
}

# Get the most likely component
gmcm.class <- apply(get.prob(pic.rgbmat, theta = seg.res.gmcm), 1, which.max)

# Extracting the colour for the estimated component
mus <- do.call(rbind, seg.res.gmcm$mu)
ranked.rgb.cols <- GMCM:::pgmm.marginal(mus, seg.res.gmcm)
rgb.cols <-
  sapply(1:3, function(i) quantile(x = pic.rgbmat[,i],
                                   prob = ranked.rgb.cols[,i],
                                   type = 1))

# Arrange into 3 dimensional array
gmcm.rgbmat <- rgb.cols[gmcm.class, ]
seg.gmcm[,,1] <- matrix(gmcm.rgbmat[,1], nn, mm)
seg.gmcm[,,2] <- matrix(gmcm.rgbmat[,2], nn, mm)
seg.gmcm[,,3] <- matrix(gmcm.rgbmat[,3], nn, mm)


#
# Segmentation using K-means clustering
#

km.file <- gsub(".jpg$", ".km.RData", fig7.file)
if (!exists("seg.res.km") | recompute) {
  system.time(
    seg.res.km <-
      kmeans(x = pic.rgbmat, centers = n.cols, iter.max = 100)
  )
  resave(seg.res.km, file = "saved.RData")
}

km.rgbmat <- seg.res.km$centers[seg.res.km$cluster, ]
seg.km[,,1] <- matrix(km.rgbmat[,1], nn, mm)
seg.km[,,2] <- matrix(km.rgbmat[,2], nn, mm)
seg.km[,,3] <- matrix(km.rgbmat[,3], nn, mm)

#
# Create segmented pictures
#

writeJPEG(seg.gmcm, target = gsub(".RData$", ".jpg", gmcm.file), quality = 0.8)
writeJPEG(seg.km,   target = gsub(".RData$", ".jpg", km.file),   quality = 0.8)



################################################################################
# The following performs the Fresh vs Frozen reproducibility analysis
# in section 5.2. Effects of cryopreservation on reproducibility
################################################################################

set.seed(1)

# Load ranked data from GMCM package
data(freshVsFrozen) # The object is named freshVsFrozen

freshfroz.tstat <- freshVsFrozen[, c(1, 3)]
freshVsFrozen$Fresh.adj.pval <-
  p.adjust(freshVsFrozen$PreVsPost.Fresh.pval,  method="BH")
freshVsFrozen$Frozen.adj.pval <-
  p.adjust(freshVsFrozen$PreVsPost.Frozen.pval, method="BH")
freshfroz.pval  <- freshVsFrozen[, c(2, 4)]
freshVsFrozen$adj.FreshVsFrozen.pval <-
  p.adjust(freshVsFrozen$FreshVsFrozen.pval, method = "BH")
# Absolute t-score and ranking
freshfroz <- Uhat(abs(freshfroz.tstat))
freshfroz.n.fits <- 40

if (!exists("best.par.freshfroz") | recompute) {
  best.ll <- -Inf
  best.par <- NULL
  best.init.par <- NULL
  for (i in 1:freshfroz.n.fits) {
    init.par <- c(runif(1, 0.05, 0.95), rchisq(2, df = 1), runif(1, 0.05, 0.95))

    # Fit model
    gmcm.par <-
      fit.meta.GMCM(u = freshfroz, init.par = init.par, method = "NM",
                    max.ite = 1000, reltol = 1e-5, verbose = FALSE)
    # Compute loglikelihood
    ll <- GMCM:::dgmcm.loglik(u = freshfroz, theta = meta2full(gmcm.par, d = 2))

    cat("Fit", i, "done.\tll =", sprintf("%.5f", ll),
        "\tbest.ll =", sprintf("%.5f", best.ll)); flush.console()
    if (ll > best.ll) {
      best.ll <- ll
      best.par <- gmcm.par
      best.init.par <- init.par
      cat("\t(best.ll updated!)"); flush.console()
    }
    cat("\n")
  }
  best.par.freshfroz <- best.par
  resave(best.par.freshfroz, file = "saved.RData")
}

# Get idr values
freshfroz.idr <- get.IDR(x = freshfroz, par = best.par.freshfroz)
freshfroz.group <- (freshfroz.idr$idr < 0.5) + (freshfroz.idr$IDR < 0.05) + 1

# Tests for dependency between significane of FreshVsFrozen compared to
# idr
suppressWarnings({
freshfroz.ctest <- cor.test(freshfroz.idr$idr, freshVsFrozen$FreshVsFrozen.pval,
                            method = "spearman")
})

freshfroz.table <- table(irreproducible = freshfroz.group == 1,
                     sig.diff = freshVsFrozen$adj.FreshVsFrozen.pval <= 0.05)
freshfroz.htest <- fisher.test(freshfroz.table)


#
# Overlap with differentially expressed genes across Fresh and Frozen
#

# Plotting and reporting
is.sig <- freshVsFrozen$adj.FreshVsFrozen.pval <= 0.05
ng <- c(table(freshfroz.group), sum(is.sig))
legend <- paste0(c("Irreproducible\nn = ",
                   "Reproducible\nn = ",
                   "Highly reproducible\nn = ",
                   "Adj. p-value < 0.05\n(Fresh vs Frozen)\nn = "), ng)

jpeg("Figure6.jpg", height = 7*0.5, width = 3*7*0.5, units = "in", res = 100)
{
  freshfroz.cols <- c("grey", "steelblue", "black", "red")
  freshfroz.cols2 <- freshfroz.cols[ifelse(is.sig, 4, freshfroz.group)]


  par(mfrow = c(1, 3),
      mar = c(3, 3.3, 2, 0.1),
      oma = c(0.1, 0.1, 0.1, 0.1),
      mgp = c(2.1, 1, 0),
      xaxs = "i", yaxs = "i",
      xpd = TRUE)
  labs <- pretty(seq(0,1))

  plot(1 - freshfroz.pval, #asp = 1,
       pch = 16, cex = 0.4, xlab = "", ylab = "", axes = FALSE,
       main = "1 - p-value", col = freshfroz.cols2)
  points(1 - freshfroz.pval[is.sig, ], pch = 15, cex = 0.6,
        col = freshfroz.cols2[is.sig])

  mtext(side = 1, expression(1 - p[list(g, scriptstyle(fresh))]),
        line = 2, cex = 0.6)
  mtext(side = 2, expression(1 - p[list(g, scriptstyle(frozen))]),
        line = 2, cex = 0.6)
  axis(1, at = labs, lwd = axiscex)
  axis(2, at = labs, lwd = axiscex)
  mtext("A", line = 0.9, adj = -0.05, cex = 1, font = 2)

  plot(freshfroz, xlab = "", ylab = "", pch = 16, cex = 0.4, #asp = 1,
       axes = FALSE, main = "Ranked values",
       col = freshfroz.cols2, xpd = TRUE)
  points(freshfroz[is.sig, ], pch = 15, cex = 0.6,
         col = freshfroz.cols2[is.sig])
  mtext(expression(hat(u)[list(g, scriptstyle(fresh))]), 1,
        line = 2, cex = 0.7, xpd = TRUE)
  mtext(expression(hat(u)[list(g, scriptstyle(frozen))]), 2,
        line = 2, cex = 0.7, xpd = TRUE)
  axis(1, at = labs, lwd = axiscex)
  axis(2, at = labs, lwd = axiscex)
  mtext("B", line = 0.9, adj = -0.05, cex = 1, font = 2)

  freshfroz.pseudo <-
    GMCM:::qgmm.marginal(u = freshfroz,
                         theta = meta2full(best.par.freshfroz, d = 2))
  plot(freshfroz.pseudo, xlab = "", ylab = "",  pch = 16, cex = 0.4, #asp = 1,
       axes = FALSE, main = "Pseudo observations",
       col = freshfroz.cols2, xpd = TRUE)
  points(freshfroz.pseudo[is.sig, ], pch = 15, cex = 0.6,
         col = freshfroz.cols2[is.sig])
  mtext(expression(hat(z)[list(g,scriptstyle(fresh))]), 1,
        line = 2, cex = 0.7,xpd = TRUE)
  mtext(expression(hat(z)[list(g,scriptstyle(frozen))]), 2,
        line = 2, cex = 0.7, xpd = TRUE)
  axis(1); axis(2)
  mtext("C", line = 0.9, adj = -0.05, cex = 1, font = 2)

  legend("topleft", pch = c(16,16,16), col = freshfroz.cols[-4], text.font = 2,
         pt.cex = rep(0.5,3)*2.5,
         inset = 0.01, legend = legend[-4],
         bg = "#FFFFFF00", bty = "n", y.intersp = 2)
  legend("bottomright", pch = 15, col = freshfroz.cols[4], text.font = 2,
         pt.cex = 0.6*2.5,
         inset = 0.01, legend = legend[4],
         bg = "#FFFFFF00", bty = "n", y.intersp = 2)
}
dev.off()

## ----start_example, echo = FALSE, eval = TRUE-----------------------
library("GMCM")
set.seed(100)
n <- 10000
sim <- SimulateGMCMData(n = n, theta = rtheta(m = 3, d = 2))

## ----start_example, echo = TRUE, eval = FALSE-----------------------
#  library("GMCM")
#  set.seed(100)
#  n <- 10000
#  sim <- SimulateGMCMData(n = n, theta = rtheta(m = 3, d = 2))

## ----make_simulation_example_png, results = 'hide', echo = FALSE----
jpeg("Figure2.jpg", width = 3*7*0.5, height = 7*0.5, units = "in", res = 100)
{
  par(mar = c(3, 3.3, 2, 0.1),
      oma = c(0, 0, 0, 0) + .1,
      mgp = c(2.1, 1, 0),
      xaxs = "i", yaxs = "i",
      xpd = TRUE)
  layout(rbind(c(3,1,2,3)), widths = c(1/6, 1/3, 1/3, 1/6))

  pchs <- c(1,3,4)
  eks.col <- c("orange", "black", "steelblue")
  labs <- pretty(seq(0,1))

  # Panel A
  plot(sim$z, col = eks.col[sim$K], main = "Latent GMM process", axes = FALSE,
       xlab = "", ylab = "", cex = 0.9, pch = pchs[sim$K])
  mtext(expression(z[1]), 1, line = 2, cex = 0.6)
  mtext(expression(z[2]), 2, line = 2, cex = 0.6)
  axis(1, lwd = axiscex)
  axis(2, lwd = axiscex)
  mtext("A", line = 0.9, adj = -0.05, cex = 1, font = 2)

  # Panel B
  plot(sim$u, col = eks.col[sim$K], main = "GMCM process", axes = FALSE,
       xlab = "", ylab = "", cex = 0.9, pch = pchs[sim$K])
  mtext(expression(u[1]), 1, line = 2, cex = 0.6)
  mtext(expression(u[2]), 2, line = 2, cex = 0.6)
  axis(1, at = labs, lwd = axiscex)
  axis(2, at = labs, lwd = axiscex)
  mtext("B", line = 0.9, adj = -0.05, cex = 1, font = 2)
}
dev.off()

## ----run_example, results = 'hide', eval = FALSE--------------------
#  ranked.data <- Uhat(sim$u)
#  start.theta <- choose.theta(ranked.data, m = 3)
#  mle.theta <- fit.full.GMCM(u = ranked.data, theta = start.theta,
#                             method = "NM", max.ite = 10000, reltol = 1e-4)
#  kappa <- get.prob(ranked.data, theta = mle.theta)
#  Khat <- apply(kappa, 1, which.max)

## ----time_run_example, results = 'hide', echo = FALSE, warning=FALSE----
st <- proc.time()
info <- capture.output({
ranked.data <- Uhat(sim$u)
start.theta <- choose.theta(ranked.data, m = 3)
mle.theta <- fit.full.GMCM(u = ranked.data, theta = start.theta,
                           method = "NM", max.ite = 10000, reltol = 1e-4)
kappa <- get.prob(ranked.data, theta = mle.theta)
Khat <- apply(kappa, 1, which.max)
})
elapsed <- round((proc.time() - st)[3], 1)
nm.ite <- as.numeric(sub(" *([0-9]+) [a-z ]+", "\\1", info[length(info)]))

## ----simulation_example_results, results = 'hide', echo = FALSE-----
fit <- mle.theta
Khat.tmp <- rep(NA, length(Khat))
Khat.tmp[Khat==1] <- 1
Khat.tmp[Khat==2] <- 2
Khat.tmp[Khat==3] <- 3
Khat <- Khat.tmp

# Plot the results
jpeg("Figure3.jpg", width = 3*7*0.5, height = 7*0.5, units = "in", res = 100)
{

  par(mfrow = c(1, 3),
      mar = c(3, 3.3, 2, 0.1),
      oma = c(0, 0, 0, 0) + 0.1,
      mgp = c(2.1, 1, 0),
      xaxs = "i", yaxs = "i",
      xpd = TRUE)
  labs <- pretty(seq(0,1))

  plot(sim$u, col = eks.col[Khat], main = "Clustering", axes = FALSE,
       xlab = "", ylab = "", cex = 0.7,
       pch = pchs[Khat])
  mtext(expression(u[1]), 1, line = 2, cex = 0.6)
  mtext(expression(u[2]), 2, line = 2, cex = 0.6)
  axis(1, at = labs, lwd = axiscex)
  axis(2, at = labs, lwd = axiscex)
  mtext("A", line = 0.9, adj = -0.05, cex = 1, font = 2)

  simfit <- SimulateGMCMData(n, theta = fit)
  plot(simfit$z, col = eks.col[simfit$K],
       main = "Model check 1", axes = FALSE,
       xlab = "", ylab = "", cex = 0.7,
       pch = pchs[simfit$K])
  mtext(expression(z[1]), 1, line = 2, cex = 0.6)
  mtext(expression(z[2]), 2, line = 2, cex = 0.6)
  axis(1, lwd = axiscex)
  axis(2, lwd = axiscex)
  mtext("B", line = 0.9, adj = -0.05, cex = 1, font = 2)

  plot(simfit$u, col = eks.col[simfit$K],
       main = "Model check 2", axes = FALSE,
       xlab = "", ylab = "", cex = 0.7, pch = pchs[simfit$K])
  mtext(expression(u[1]), 1, line = 2, cex = 0.6)
  mtext(expression(u[2]), 2, line = 2, cex = 0.6)
  axis(1, at = labs, lwd = axiscex)
  axis(2, at = labs, lwd = axiscex)
  mtext("C", line = 0.9, adj = -0.05, cex = 1, font = 2)
}
dev.off()

## ----confusion_table, results = 'asis', echo = FALSE----------------
confusion.tab <- table(K = sim$K, Khat)
acc <- round(sum(diag(confusion.tab))/n,3)*100
km.tmp  <- km <- kmeans(ranked.data, centers = 3)$cluster

# Rename clusters to get correct colours and
# make the confusion matrix easy to read
km.tmp[km == 1] <- 2
km.tmp[km == 2] <- 3
km.tmp[km == 3] <- 1
km <- km.tmp
confusion.tab.km <- table(K = sim$K, Khat = km)
acc.km <- round(sum(apply(confusion.tab.km, 1, max))/n, 3)*100

latex(cbind(confusion.tab, confusion.tab.km),
      cgroup = c("$\\hat{H}$ (GMCM)", "$\\hat{H}$ ($k$-means)"),
      n.cgroup = c(3,3),
      rgroup = "$H$",
      n.rgroup = 3,
      title = "",
      file = "",
      label = "confusion.mat",
      caption.loc = "bottom",
      caption = "Confusion matrices of GMCM and $k$-means clustering results.")

## ----table_speed_res, tidy = FALSE, echo = FALSE, results='asis'----
# Formatting the results stored in the speed.res object
speed.res <-
  cbind(Package = toupper(gsub("-pkg (PEM|NM)", "", rownames(speed.res))),
        Algorithm = gsub("(gmcm|idr)-pkg ", "", rownames(speed.res)),
        as.data.frame(speed.res, row.names = 1:nrow(speed.res)))
speed.res <- as.data.frame(speed.res)
speed.res$elapsed <- round(speed.res$elapsed, 2)
speed.res$"Rel.\\ speed" <- round(speed.res$"Rel.\\ speed", 1)
colnames(speed.res) <-  # Some extra formatting of the table
  gsub("elapsed", "Runtime ($s$)",
       gsub("^it$", "Iterations ($n$)",
            colnames(speed.res)))
speed.res[, c(4:7, 10)] <-
  round(speed.res[, c(4:7, 10)],3)
latex(speed.res[,-c(1,3,4:7)][,c(1,3,2,4,5)],
      title = "$p$ / Package",
      file = "",
      caption.loc = "bottom",
      caption = paste("Runtime comparisons of the \\pkg{idr} and",
                      "\\pkg{GMCM} packages with increasing number",
                      "of observations $p$. The benchmarked",
                      "optimization procedures are the pseudo EM algorithm ",
                      "(PEM) and the Nelder-Mead (NM) method.",
                      "The runtime is given in seconds.",
                      "The last column shows the relative speed per",
                      "iteration compared to the fastest procedure."),
      label = "speed.tab",
      rowname = paste("\\pkg{", gsub("IDR", "idr", speed.res$Package), "}"),
      rgroup = unique(format(speed.res$n.obs,
                             scientific = FALSE, big.mark = ",")))

## ----table_equivalent_optima, results='asis',echo=FALSE-------------
tab <- matrix("$\\cdot$", 3, 4)
colnames(tab) <- c("\\alpha_1", "\\mu", "\\sigma", "\\rho")
colnames(tab) <- paste("$", colnames(tab), "$", sep = "") #"\\approx$"
rownames(tab) <- 1:3
tab[1,1] <- "$1$"
tab[2,c(1,4)] <- c("$0$", "$0$")
tab[3,-1] <- c("$0$", "$1$", "$0$")
caption <- paste("Equivalent optima in pure noise. A dot ($\\cdot$) denotes",
                 "an arbitrary value. The given values need only to be",
                 "approximate.")
latex(tab, title = "Situation",
      caption = caption,
      label = "BadEst",
      file = "",
      caption.loc = "bottom")

## ----freshfroz_group_pct, echo=FALSE, results='hide'----------------
tmp <- table(freshfroz.group)
freshfroz.group.pct <-
  paste0("$", prettyN(tmp), "$ $(", 100*round(tmp/sum(tmp), 3), "\\%)$")
freshfroz.group.pct2 <-
  paste0("$", prettyN(sum(tmp[2:3])), "$ $(",
         100*round(sum(tmp[2:3])/sum(tmp), 3), "\\%)$")

## ----sessionInfo, echo=FALSE, results='asis'------------------------
toLatex(sessionInfo())

