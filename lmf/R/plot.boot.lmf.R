plot.boot.lmf <-
function(x,
                          what = c("all"),
                          ...)
{
  #Check if x is of class "boot.lmf"
  if(class(x) != "boot.lmf")
    stop("'x' does not contain bootstrap replicates from a 'lmf' model")
  #Plot the parameter distribution from the bootstrap replicates for "leslie"
  if(!is.null(x$lboot) & (what == "projection" | what == "all"))
  {
    #Plot lambda
    plot(density(x$luvboot[, 1]), main = "Growth rate", xlab =
           expression(lambda), col = "blue")
    #Set the option to click to see next plot and the plot window to have
    #nage numbers of plots
    par(ask = TRUE, mfrow = c(1, x$nage))
    #Plot stable age distribution
    for(i in 1 : x$nage)
    {
      plot(density(x$luvboot[, 1 + i]), main = "Stable age distribution",
           xlab = bquote(u[.(x$uage[i])]), col = "blue")
    }
    #Plot reproductive values
    for(i in 1 : x$nage)
    {
      plot(density(x$luvboot[, x$nage + 1 + i]), main = "Reproductive value",
           xlab = bquote(v[.(x$uage[i])]), col = "blue")
    }
    #Set the plot window to have "nage^2" number of plots (maximum 4 * 4 number of plots)
    par(mfrow = pmin(dim(x$l), c(4, 4)))
    #Plot projection matrix
    for(i in 1 : ifelse(x$nage == 1, dim(x$l)[1], dim(x$l)[1]^2))
    {
      plot(density(sapply(x$lboot, '[', i)),
           main = "Projection matrix (l)", xlab = paste("Component ", i,
                                                        sep = ""), col = "blue")
    }
  }
  #Plot the parameter distribution from the bootstrap replicates for "alpha"
  if(!is.null(x$aMboot) & (what == "alpha" | what == "all"))
  {
    #Set the plot window to have one plots
    par(mfrow = c(1, 1))
    #Plot environmenal variance
    plot(density(x$eboot), main = "Environmental variance",
         xlab = bquote(sigma[e]^2), col = "blue")
    #Set the option to click to see next plot and the plot window to have
    #nage + 1 numbers of plots
    par(ask = TRUE, mfrow = c(1, x$nage + 1))
    #Plot demographic variances
    for(i in 1 : x$nage)
    {
      plot(density(x$djboot[, i]), main = "Demographic variance",
           xlab = bquote(sigma[paste("d", .(x$uage[i]), sep = "")]^2), col = "blue")
    }
    plot(density(x$dboot), main = "Demographic variance",
         xlab = bquote(sigma[d]^2), col = "blue")
    #Set the plot window to have "npar" number of plots
    par(mfrow = c(1, x$npar))
    #Plot temporal mean alphas
    for(i in 1 : x$npar)
    {
      plot(density(sapply(x$aMboot, '[', i)),
           main = expression("Temporal"~alpha~"-estimates"),
           xlab = bquote(alpha[.(i - 1)]~"(M)"),
           col = "blue")
    }
    #Set the plot window to have "npar^2" number of plots
    #(maximum 4 * 4 number of plots)
    par(mfrow = c(min(4, x$npar), min(4, x$npar)))
    #Plot components of the temporal covariance matrix (M)
    for(i in 1 : x$npar^2)
    {
      plot(density(sapply(x$Mboot, '[', i)),
           main = "Temporal covariance matrix (M)", xlab = paste("Component ", i,
                                                                 sep = ""), col = "blue")
    }
    
    #Set the plot window to have "npar" number of plots
    par(mfrow = c(1, x$npar))
    #Plot temporal mean alphas
    for(i in 1 : x$npar)
    {
      plot(density(sapply(x$anfboot, '[', i)),
           main = expression("Temporal"~alpha~"-estimates (M = 0)"),
           xlab = bquote(alpha[.(i - 1)]~"(M = 0)"),
           col = "blue")
    }
    #Set the plot window to have "npar^2" number of plots
    #(maximum 4 * 4 number of plots)
    par(mfrow = c(min(4, x$npar), min(4, x$npar)))
    #Plot components of the temporal covariance matrix (M)
    for(i in 1 : x$npar^2)
    {
      plot(density(sapply(x$Anfboot, '[', i)),
           main = "Temporal covariance matrix (M = 0)", xlab = paste("Component ", i,
                                                                     sep = ""), col = "blue")
    }
  }
}
