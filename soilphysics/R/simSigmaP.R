simSigmaP <- 
function(voidratio, stress,
	what.out = c("sigmaP", "CI", "SI"),
	method = c("casagrande", "VCLzero", "reg1", "reg2", 
	"reg3", "reg4", "pacheco"),
	n4VCL = 3, nsim = 100)
{
    what.out <- match.arg(what.out)

    # mean and vcov of the simulation from p-normal
    xy <- data.frame(x = log10(stress), y = voidratio)
    fit <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4), data = xy)
    beta <- coef(fit)
    mcov <- vcov(fit)

    # parameter simulation
    sim <- mvrnorm(nsim, beta, mcov)

    # prediction (voidratio) auxiliar function
    fpred <- function(b, x = log10(stress))
    {
       pred <- b[1] + b[2]*x + b[3]*x^2 + b[4]*x^3 + b[5]*x^4
       return(pred)
    }

    # calculation of simulated sigmaP
    pb <- tkProgressBar(title = "Simulating Preconsolidation Stress", 
       label = "SIMULATION PROGRESS",
       min = 0, max = nsim * length(method), width = 400L)

    msigmaP <- matrix(nrow = nsim, ncol = length(method))
    colnames(msigmaP) <- method
    for(j in 1:length(method)) {
       i = 1
       repeat{
          pred <- fpred(sim[i, ])
          msigmaP[i, j] <- sigmaP(pred, stress, method = method[j],
             n4VCL = n4VCL, graph = FALSE)[[what.out]]
          k <- i + nsim*(j - 1)
          setTkProgressBar(pb, k, 
             label = sprintf("SIMULATION PROGRESS (%.0f%%)", 
             100 * k / (nsim * length(method))))
          i = i + 1
          if (i > nsim) break()  
       }
    }
    Sys.sleep(0.5)
    close(pb)

    # output
    class(msigmaP) <- "simSigmaP"
    return(msigmaP)
}
