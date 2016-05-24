# 2014-09-21 CJS change all Inf to NA
# 2012-01-22 CJS made X/Y axis limits the same so p-value prints properly
# 2011-06-13 CJS returned bayesian p-values

PredictivePosteriorPlot.TSPDE.WHCH <- function( discrep  ) {
#   Given the discrepancy measures, creates a set of panel plots.
#   It is assumed that the bp has 12 columns
#     (1-2)  o,s Freeman-Tukey measures for m2
#     (3-4)  o,s Freeman-Tukey measures for u2.A
#     (5-6)  o,s Freeman-Tukey measures for u2.N
#     (7-8)  o,s Freeman-Tukey for m2+u2.A+u2.N

# Change any Inf to NA
temp <- discrep == Inf | discrep == -Inf
if(sum(temp)>0){cat(sum(temp), " infinite discrepancy measures set to NA\n")}
discrep[ temp ] <- NA

#browser()
split.screen(figs=c(2,2))  # 4 rows and 2 columns
titles <- c("Freeman-Tukey for m2", 
            "Freeman-Tukey for u2.A", 
            "Freeman-Tukey for u2.N", 
            "Total Freeman-Tukey")

saved_p_values <- rep(NA, length(titles))

for(i in 1:4){
  screen(i)
  par(cex=.5)
  par(mai=c(.40,.40,.40,.40)) # margins of plot relative to plotting position

  lims <- range(discrep[,(2*i):(2*i-1)], na.rm=TRUE)

  plot(discrep[,(2*i):(2*i-1)],
       xlab="Simulated", ylab="Observed", 
       main=titles[i], cex.main=1.5, xlim=lims, ylim=lims)
  abline(a=0, b=1)

  ## Compute Bayesian p-value
  p.value <- sum(discrep[,2*i-1]<discrep[,2*i],na.rm=TRUE)/nrow(discrep)
  saved_p_values[i] <- p.value
  
  ## Add p-value to plot
  x.loc <- mean(lims, na.rm=TRUE)
  y.loc <- min(lims,  na.rm=TRUE)
  
  text(x.loc, y.loc,
       labels=paste("Bayesian GOF P:",formatC(p.value, digits=2, format="f")),
       cex=1.5, adj=c(0,0))  
}
close.screen(all.screens=TRUE)     # exit from plots  
gof <- data.frame(statistic=titles, p.value=saved_p_values)
gof
}
