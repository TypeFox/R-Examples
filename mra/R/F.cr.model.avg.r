
F.cr.model.avg <- function(fits = ls(pattern = "^fit"), what = "survival",
                              fit.stat = "qaicc") 
{

# F.cr.model.avg.JB.R
# Perform model averaging on a list of capture-recapture models.
# Code modified by Jeff Bromaghin on 7 August 2012 to reduce memory overhead.



# Check value of what.
if(substring(what, 1, 1) == "s")
  {
  x.name <- "s.hat"
  x.se <- "se.s.hat"
  want.n.hat <- FALSE
  } else if(substring(what, 1, 1) == "c")
  {
  x.name <- "p.hat"
  x.se <- "se.p.hat"
  want.n.hat <- FALSE
  } else if(substring(what, 1, 1) == "n")
  {
  x.name <- "n.hat"
  x.se <- "se.n.hat"
  want.n.hat <- TRUE
  } else
  {
  stop(paste("Invalid option. Cannot model average '", what, "'.", sep = ""))
  }



# Check value of fit.stat
if( !(fit.stat == "aicc") & !(fit.stat == "qaicc"))
  {
  stop(paste("Invalid option. Cannot model average '", fit.stat, "'.",
             sep = ""))
  }



# Go through the model list once to determine how many models will be used in
# model averaging.  Also record nan and ns from each model to make sure they are
# all equal.
n.mod <- length(fits)
use.mod <- vector("logical", n.mod)
nan.mod <- vector("numeric", n.mod)
ns.mod <- vector("numeric", n.mod)
for(li1 in 1:n.mod)
  {
  
  # Retrieve model from memory.
  fit <- get(fits[li1], pos = .GlobalEnv)
  
  # Should this model be used?
  if("cr" %in% class(fit))
    {
    if((fit$exit.code == 1) & (fit$cov.code == 0) & (fit$df > 0))
      {
      use.mod[li1] <- TRUE
      nan.mod[li1] <- fit$aux$nan
      ns.mod[li1] <- fit$aux$ns
      } else
        {
        use.mod[li1] <- FALSE
        nan.mod[li1] <- NA
        ns.mod[li1] <- NA
        }
    } else
    {
    warning(paste("Object", fits[li1], "in fits is not a CR object and has been ignored."))
    use.mod[li1] <- FALSE
    nan.mod[li1] <- NA
    ns.mod[li1] <- NA
    }
  }



# Check that nan and ns are all equal for all model objects.
dum1 <- min(nan.mod, na.rm=TRUE)
dum2 <- max(nan.mod, na.rm=TRUE)
if(dum1 == dum2)
  {
  nan <- dum1
  } else
  {
  stop(paste("Number of individuals differ among models. Cannot model average."))
  }

dum1 <- min(ns.mod, na.rm=TRUE)
dum2 <- max(ns.mod, na.rm=TRUE)
if(dum1 == dum2)
  {
  ns <- dum1
  } else
  {
  stop(paste("Number of occasions differ among models. Cannot model average."))
  }
rm(dum1, dum2, nan.mod, ns.mod)



# Store dimension names
dim.nms <- dimnames(fit$histories)



# Determine number of models to use in model averaging and allocate required
# memory.  
n.mod.good <- sum(use.mod)
if(n.mod.good < 1)
  {
  stop(paste("Number of good models equals 0."))
  }
if(want.n.hat)
  {
  nan <- 1
  }
n.stats <- nan*ns
fits <- fits[use.mod]
stats <- matrix(0, n.mod.good, n.stats)
se.stats <- stats
all.fit.stat <- vector("numeric", n.mod.good)



# Go through model list a second time and pull out requested statistics.
for(li1 in 1:n.mod.good)
  {

  # Retrieve model from memory.
  fit <- get(fits[li1], pos = .GlobalEnv)
  
  # Load statistics
  stats[li1, 1:n.stats] <- unlist(fit[x.name])
  se.stats[li1, 1:n.stats] <- unlist(fit[x.se])
  all.fit.stat[li1] <- unlist(fit[fit.stat])
  }


# Compute model weights
delta.AIC <- all.fit.stat - min(all.fit.stat)
dum1 <- exp(-0.5*delta.AIC)
dum2 <- sum(dum1)
wi.array <- dum1/dum2
rm(dum1, dum2)



# Compute model-averaged real parameters.
a1 <- stats*wi.array
theta.average <- apply(a1, 2, sum)



# Compute model-averaged standard errors.
var.theta <- se.stats^2
dum <- matrix(theta.average, nrow=n.mod.good, ncol=n.stats, byrow=TRUE)
a2 <- sqrt(var.theta + (stats - dum)^2)*wi.array
rm(dum)
se.theta.average <- apply(a2, 2, sum)



# Compute conditional standard errors.
a2 <- se.stats*wi.array
se.conditional.theta.average <- apply(a2, 2, sum)



# Compute proportion of variance due to model selection uncertainty.
mod.selection.proportion <- (se.theta.average - se.conditional.theta.average)/
                             se.theta.average


                             
# Load model-averaged estimates into nan by ns matrices.
hat <- matrix(theta.average, nrow=nan, ncol=ns, byrow=FALSE)
se.hat <- matrix(se.theta.average, nrow=nan, ncol=ns, byrow=FALSE)
se.hat.conditional <- matrix(se.conditional.theta.average, nrow=nan, ncol=ns,
                             byrow=FALSE)
mod.selection.proportion <- matrix(mod.selection.proportion, nrow=nan, ncol=ns,
                                   byrow=FALSE)



# Construct model fit summary table.
AIC.table <- data.frame(fits, all.fit.stat, delta.AIC, wi.array,
                        stringsAsFactors=FALSE)
AIC.table <- AIC.table[order(AIC.table[,3]),]
names(AIC.table) <- c("model", fit.stat, paste("delta.", fit.stat, sep=""),
                       paste(fit.stat, ".weight", sep=""))



# Construct results list.
if(substring(what, 1, 1) == "s")
  {
  dimnames(hat) <- dim.nms                                   
  dimnames(se.hat) <- dim.nms                                   
  dimnames(se.hat.conditional) <- dim.nms                                   
  dimnames(mod.selection.proportion) <- dim.nms                                   
  results <- list(fit.table=AIC.table, s.hat=hat, se.s.hat=se.hat,
                  se.s.hat.conditional=se.hat.conditional,
                  mod.selection.proportion=mod.selection.proportion)
  } else if(substring(what, 1, 1) == "c")
  {
  dimnames(hat) <- dim.nms                                   
  dimnames(se.hat) <- dim.nms                                   
  dimnames(se.hat.conditional) <- dim.nms                                   
  dimnames(mod.selection.proportion) <- dim.nms                                   
  results <- list(fit.table=AIC.table, p.hat=hat, se.p.hat=se.hat,
                  se.p.hat.conditional=se.hat.conditional,
                  mod.selection.proportion=mod.selection.proportion)
  } else
  {
  names(hat) <- dim.nms[[2]]
  names(se.hat) <- dim.nms[[2]]
  names(se.hat.conditional) <- dim.nms[[2]]
  names(mod.selection.proportion) <- dim.nms[[2]]
  results <- list(fit.table=AIC.table, n.hat=hat, se.n.hat=se.hat,
                  se.n.hat.conditional=se.hat.conditional,
                  mod.selection.proportion=mod.selection.proportion)
  results$n.hat.lower <- results$n.hat - 1.96*results$se.n.hat  
  results$n.hat.upper <- results$n.hat + 1.96*results$se.n.hat  
  results$nhat.v.meth <- fit$nhat.v.meth + 3
  results$n.hat.conf <- 0.95
  results$intervals <- fit$intervals
  class(results) <- c("nhat", "cr")
  }

results
}
