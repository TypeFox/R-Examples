# 2014-09-01 CJS Change Inf to NA 
# 2012-01-22 CJS Made X/Y axis limits the same so that Bayesian p-value prints properly
# 2011-06-13 CJS returned p-values
# 2010-03-29 CJS First creation of routine

PredictivePosteriorPlot.TSPDE.WHCH2 <- function( discrep  ) {
#   Given the discrepancy measures, creates a set of panel plots.
#   It is assumed that the discrepancy measure has 16 columns for the Bayesian p-value plot
#     ( 1- 2)  o,s Freeman-Tukey measures for m2
#     ( 3- 4)  o,s Freeman-Tukey measures for u2.A.YoY
#     ( 5- 6)  o,s Freeman-Tukey measures for u2.N.YoY
#     ( 7- 8)  o,s Freeman-Tukey measures for u2.A.1
#     ( 9-10)  o,s Freeman-Tukey measures for u2.N.1
#     (11-12)  o,s Freeman-Tukey for u2.A.YoY+u2.N.YoY
#     (13-14)  o,s Freeman-Tukey for u2.A.1  +u2.N.1
#     (15-16)  o,s Freeman-Tukey for all data (m2, YoY and Age 1)`

# Change any Inf to NA
temp <- discrep == Inf | discrep == -Inf
if(sum(temp)>0){cat(sum(temp), " infinite discrepancy measures set to NA\n")}
discrep[ temp ] <- NA

#browser()
titles <- c("Freeman-Tukey for m2", 
            "Freeman-Tukey for u2.A.YoY", 
            "Freeman-Tukey for u2.N.YoY", 
            "Freeman-Tukey for u2.A.1", 
            "Freeman-Tukey for u2.N.1", 
            "Freeman-Tukey for YoY",
            "Freeman-Tukey for Age 1",
            "Total Freeman-Tukey")
saved_p_values <- rep(NA, length(titles))

for(page in 1:2){
   split.screen(figs=c(2,2))  # 4 rows and 2 columns
   for(i in 1:4){
     screen(i)
     par(cex=.5)
     par(mai=c(.40,.40,.40,.40)) # margins of plot relative to plotting position

     ## Compute plotting limits to generate symmetric plot
     lims <- range(discrep[,(2*i):(2*i-1)])
     
     ## Plot observed vs simulated discrepancies
     plot(discrep[,((page-1)*4+(2*i)):((page-1)*4+(2*i-1))],
          xlab="Simulated", ylab="Observed", 
         main=titles[(page-1)*4+i], cex.main=1.5, xlim=lims, ylim=lims)
     abline(a=0, b=1)

     ## Compute Bayesian p-value
     p.value <- sum(discrep[,(page-1)*4+2*i-1]<discrep[,(page-1)*4+2*i], na.rm=TRUE)/nrow(discrep)
     saved_p_values[i] <- p.value
     
  ## Add p-value to plot
  x.loc <- mean(lims)
  y.loc <- min(lims)
  
  text(x.loc, y.loc,
       labels=paste("Bayesian GOF P:",formatC(p.value, digits=2, format="f")),
       cex=1.5, adj=c(0,0))  
   }
   close.screen(all.screens=TRUE)     # exit from plots for this page
   gof <- data.frame(statistic=titles, p.value=saved_p_values)
   gof
 }
}
