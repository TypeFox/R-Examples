trace_plot <- function(title=" ", results=NULL, parms_to_plot=NULL, panels=c(1,1),
   mai=if(prod(panels)>1){c(.4,.4,.4,.4)} else {c(1.02,0.82,0.82,0.42)},
   cex=if(prod(panels)>1){.5} else {1}
   ) {
#
# Takes the MCMC object from the fit (could be TPSDE etc), a list of parameters and produces
# the traceplots.
#   
# title - title of the plot
# results - the MCMC object containing the necessary information
# parms_to_plot - character vector containing the names of the parms to plot
#                e.g. c("logitP[1]", "logitP[2]")
#               - this should be an exact match
# panels  - how the plotting page should be arranged
# mai      - how big to make the margins around the plots
# cex      - character expansion factor
#

varnames <- names(results$sims.array[1,1,])

index <- match(parms_to_plot, varnames) # find where these parms exist in the array
plots_per_page <- prod(panels)

for(i in seq(1,length(index),plots_per_page)){  # plot potentially multiple plots/page
   split.screen(panels)
   par(new=TRUE)
   for(j in i:min(i+plots_per_page-1,length(index))){
      screen(j-i+1)
      par(mai=mai)
      par(cex=cex)
      matplot(results$sims.array[,,index[j]], type="l", 
              main=paste(title), lty=1,
              ylab='Estimate')
      text(x=1, y=max(results$sims.array[,,index[j]]), 
          label=varnames[index[j]], adj=c(0,1), cex=if(prod(panels)>1){2} else {1})
      if(results$n.chains >1){   # Only print Rhat if number of chains is >1
         text(x=dim(results$sims.array)[1]  ,y=max(results$sims.array[,,index[j]]), 
          label=paste("Rhat=",round(results$summary[index[j],"Rhat"],1)), adj=c(1,1), 
          cex=if(prod(panels)>1){2} else {1})
      }
      if(results$n.chains == 1){
         text(x=dim(results$sims.array)[1]  ,y=max(results$sims.array[,,index[j]]), 
          label="No RHat avail", adj=c(1,1), 
          cex=if(prod(panels)>1){2} else {1})
      }
   }
   close.screen(all.screens=TRUE)
}

return(NULL) # nothing to return
} # end of function
